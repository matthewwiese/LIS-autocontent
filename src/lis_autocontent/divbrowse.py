"""Build a Divbrowse docker-compose.yml, one service per diversity collection.

A Traefik proxy routes to each service by its divbrowse.host and divbrowse.path labels.
A collection that can't supply VCF_URL, GFF3_URL and CHROM_PATTERN unambiguously fails
the build with its reasons. Hosts and combined references come from divbrowse.yml.
"""

import functools
import os
import re

import yaml

HEADER = """# Divbrowse Docker Compose
#
# Just run: docker compose up -d
#
# On first start, each Divbrowse service will automatically:
#   1. Download VCF and GFF3 from the specified URLs
#   2. Convert VCF to Zarr format
#   3. Generate configuration
#   4. Start the server
#
# Data is persisted in ./data/<name>/ directories.
#
# Routing: the `proxy` service (Traefik) listens on PROXY_PORT (default 8080)
# and routes http://<divbrowse.host>/<divbrowse.path>/ to each Divbrowse
# service, based on the two labels set on it. To add a dataset, add a service
# using the x-divbrowse template with those labels and a depends_on on the last
# service, then `docker compose up -d`.
# Per-host root redirects live in traefik/dynamic/divbrowse.yml.
#
# Startup: Divbrowse services start one at a time, each once the one before it
# is healthy, so first-time setups never run at once. `docker compose up` waits
# through them; later starts skip setup and are quick.

x-divbrowse: &divbrowse
  build:
    context: .
    dockerfile: Dockerfile
  restart: unless-stopped
  healthcheck:
    test: ["CMD", "wget", "-q", "-O", "/dev/null", "http://localhost:8080/"]
    interval: 30s
    timeout: 10s
    retries: 3
    # First-time setup (download, VCF to Zarr) runs within this period.
    start_period: 2h
    start_interval: 30s

services:
  proxy:
    image: traefik:v3
    ports:
      - "${PROXY_PORT:-8080}:80"
    volumes:
      # Rootless Docker: set DOCKER_SOCK=$XDG_RUNTIME_DIR/docker.sock
      - ${DOCKER_SOCK:-/var/run/docker.sock}:/var/run/docker.sock:ro
      - ./traefik:/etc/traefik:ro
    restart: unless-stopped

"""

SERVICE = """  {name}:
    <<: *divbrowse
{depends_on}    environment:
      - VCF_URL={vcf_url}
      - GFF3_URL={gff3_url}
      - CHROM_PATTERN={chrom_pattern}
    labels:
      divbrowse.host: {host}
      divbrowse.path: {path}
    volumes:
      - ./data/{data_dir}:/opt/divbrowse
"""

DEPENDS_ON = """    depends_on:
      {previous}:
        condition: service_healthy
"""

CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "divbrowse.yml")

# Characters a compose service name, a data directory and an unquoted YAML value can
# all carry. Every diversity identifier in the store fits.
IDENTIFIER = re.compile(r"[A-Za-z0-9._-]+")
# One dot-separated name token: a species abbreviation, strain or assembly version.
# Strains can contain hyphens (IT97K-499-35), which are literal in CHROM_PATTERN.
TOKEN = re.compile(r"[A-Za-z0-9_-]+")
# A chromosome prefix, possibly dotted (Aradu.A). Some READMEs list chromosome names
# instead ("chr,Pt,Mt"), which cannot be turned into a pattern.
PREFIX = re.compile(r"[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)*")
# A hostname both a Traefik Host() rule and an unquoted YAML value can carry.
HOSTNAME = re.compile(r"[A-Za-z0-9.-]+")


class DivbrowseError(Exception):
    """Raised when the requested collections can't all be expressed as services."""


def _strings(value):
    return isinstance(value, list) and all(isinstance(v, str) for v in value)


@functools.lru_cache(maxsize=None)
def load_config(path=CONFIG_PATH):
    """divbrowse.yml's hosts and combined_references; DivbrowseError if unusable."""
    try:
        with open(path, encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as err:
        raise DivbrowseError(f"cannot read {path}: {err}") from err
    config = config if isinstance(config, dict) else {}
    hosts = config.get("hosts")
    references = config.get("combined_references") or {}
    problems = []
    if not (
        isinstance(hosts, dict)
        and all(isinstance(k, str) and isinstance(v, str) for k, v in hosts.items())
    ):
        problems.append("hosts must map each genus to a hostname")
    if not isinstance(references, dict):
        problems.append("combined_references must be a mapping")
        references = {}
    for name, reference in references.items():
        if not (
            isinstance(reference, dict)
            and _strings(reference.get("genomes"))
            and reference["genomes"]
            and _strings(reference.get("note", []))
        ):
            problems.append(
                f"combined reference {name} needs a genomes list, and a note list"
                " if any"
            )
    if problems:
        raise DivbrowseError(f"{path}: {'; '.join(problems)}")
    return {"hosts": hosts, "combined_references": references}


def service_name(identifier):
    """Compose service name; lowercase, since compose rejects uppercase image names."""
    return f"divbrowse-{identifier.lower()}"


def data_directory(identifier):
    """Host directory for a service's data, named for its collection."""
    return identifier


def reference_name(identifier):
    """The reference an identifier starts with: Wm82.gnm4, or aradu1_araip1.gnm1."""
    return ".".join(identifier.split(".")[:2])


def reference_genomes(index, collection):
    """(genomes, problems) for the genomes a diversity collection is called against."""
    combined = load_config()["combined_references"].get(
        reference_name(collection.collection_key)
    )
    if combined:
        genomes = [index.collections.get(path) for path in combined["genomes"]]
        missing = [p for p, g in zip(combined["genomes"], genomes) if g is None]
        if missing:
            return [], [f"its combined reference lacks {', '.join(missing)}"]
        return genomes, []
    if len(collection.derived_from) == 1:
        genome = index.collections.get(
            f"{collection.genus}/{collection.species}/genomes/{collection.derived_from[0]}"
        )
        if genome is not None:
            return [genome], []
    return [], ["it is not linked to a genome assembly"]


def genome_parts(index, genome):
    """(gff3_url, escaped chromosome-name stem, problems) for one reference genome."""
    problems = []
    annotations = [
        a
        for a in index.collections.values()
        if a.collection_type == "annotations"
        and a.genus == genome.genus
        and a.species == genome.species
        and a.derived_from == [genome.collection_key]
    ]
    gff3 = []
    if len(annotations) != 1:
        listed = (
            f": {', '.join(a.collection_key for a in annotations)}"
            if annotations
            else ""
        )
        problems.append(
            f"its assembly {genome.collection_key} has {len(annotations)} annotations{listed}"
        )
    else:
        gff3 = [
            r.relative_path
            for r in annotations[0].data_files
            if r.canonical_type == "gene_models_main"
            and r.relative_path.endswith(".gff3.gz")
        ]
        if len(gff3) != 1:
            problems.append(
                f"annotation {annotations[0].collection_key} has {len(gff3)} "
                "gene_models_main GFF3 files"
            )

    abbrev = genome.metadata.get("scientific_name_abbrev")
    prefix = genome.metadata.get("chromosome_prefix")
    tokens = genome.collection_key.split(".")[:2]
    if not (isinstance(abbrev, str) and TOKEN.fullmatch(abbrev)):
        problems.append(
            f"genome {genome.collection_key} has no usable scientific_name_abbrev"
        )
    if len(tokens) != 2 or not all(TOKEN.fullmatch(t) for t in tokens):
        problems.append(
            f"genome {genome.collection_key} has no <strain>.<gnm> in its identifier"
        )
    if not (isinstance(prefix, str) and PREFIX.fullmatch(prefix)):
        problems.append(
            f"genome {genome.collection_key} has chromosome_prefix {prefix!r}, "
            "which is not a single prefix"
        )
    if problems:
        return None, None, problems
    strain, gnm = tokens
    return (
        f"{index.datastore_url}/{annotations[0].path}/{gff3[0]}",
        f"{abbrev}.{strain}.{gnm}.{prefix}".replace(".", r"\."),
        [],
    )


def reference_parts(index, collection):
    """(gff3_url, chrom_pattern, problems) for a collection's reference."""
    genomes, problems = reference_genomes(index, collection)
    gff3_urls, stems = [], []
    for genome in genomes:
        gff3_url, stem, genome_problems = genome_parts(index, genome)
        problems.extend(genome_problems)
        gff3_urls.append(gff3_url)
        stems.append(stem)
    if problems:
        return None, None, problems
    if len(stems) == 1:
        return gff3_urls[0], stems[0] + "[0-9]+", []
    # divbrowse greps contig names for substrings; anchor to the chromosomes alone.
    return " ".join(gff3_urls), f"^({'|'.join(stems)})[0-9]+$", []


def select_vcf(collection):
    """(vcf, problems): the collection's one VCF, or its single .main one of several."""
    vcfs = [
        r.relative_path
        for r in collection.data_files
        if r.relative_path.endswith(".vcf.gz")
    ]
    # A single .main VCF is the full dataset; any others are subsets of it.
    main = [vcf for vcf in vcfs if vcf.endswith(".main.vcf.gz")]
    if len(vcfs) > 1 and len(main) == 1:
        vcfs = main
    if len(vcfs) != 1:
        listed = f": {', '.join(vcfs)}" if vcfs else ""
        return None, [f"it has {len(vcfs)} VCF files{listed}"]
    return vcfs[0], []


def service_environment(index, identifier, hosts):
    """(environment, problems) for one diversity collection; hosts maps genus -> host.

    environment is None whenever there are problems, so nothing is built on a guess.
    """
    matches = [
        c
        for c in index.collections.values()
        if c.collection_type == "diversity" and c.collection_key == identifier
    ]
    if not matches:
        return None, ["no diversity collection has this identifier"]
    if len(matches) > 1:
        paths = ", ".join(c.path for c in matches)
        return None, [f"several diversity collections have this identifier: {paths}"]
    (collection,) = matches

    problems = []
    if not IDENTIFIER.fullmatch(identifier):
        problems.append("its identifier has characters a compose file can't carry")
    host = hosts.get(collection.genus)
    if host is None:
        problems.append(
            f"no Divbrowse host for genus {collection.genus}; "
            f"pass --host {collection.genus}=<hostname>"
        )
    elif not HOSTNAME.fullmatch(host):
        problems.append(f"its host {host!r} is not a hostname")
    vcf, vcf_problems = select_vcf(collection)
    problems.extend(vcf_problems)

    gff3_url, chrom_pattern, reference_problems = reference_parts(index, collection)
    problems.extend(reference_problems)
    if problems:
        return None, problems
    return {
        "vcf_url": f"{index.datastore_url}/{collection.path}/{vcf}",
        "gff3_url": gff3_url,
        "chrom_pattern": chrom_pattern,
        "host": host,
    }, []


def compose_file(index, identifiers, hosts=None):
    """The docker-compose.yml text, one service per identifier in the order given.

    Raises DivbrowseError listing every problem with every collection.
    """
    problems = {}
    blocks = []
    seen = set()
    noted = set()
    previous = None
    for identifier in identifiers:
        if identifier in seen:
            problems.setdefault(identifier, []).append(
                "it was requested more than once"
            )
            continue
        seen.add(identifier)
        environment, reasons = service_environment(
            index, identifier, load_config()["hosts"] if hosts is None else hosts
        )
        if reasons:
            problems.setdefault(identifier, []).extend(reasons)
            continue
        block = SERVICE.format(
            name=service_name(identifier),
            depends_on=DEPENDS_ON.format(previous=previous) if previous else "",
            path=identifier,
            data_dir=data_directory(identifier),
            # Compose interpolates $; $$ is a literal one.
            **{key: value.replace("$", "$$") for key, value in environment.items()},
        )
        reference = reference_name(identifier)
        note = load_config()["combined_references"].get(reference, {}).get("note")
        if note and reference not in noted:
            noted.add(reference)
            block = "".join(f"  # {line}\n" for line in note) + block
        blocks.append(block)
        previous = service_name(identifier)
    if not identifiers:
        problems[""] = ["no collections were requested"]
    if problems:
        raise DivbrowseError(
            "\n".join(
                f"{identifier or 'request'}: {'; '.join(reasons)}"
                for identifier, reasons in problems.items()
            )
        )
    return HEADER + "\n".join(blocks)
