"""Build a Divbrowse docker-compose.yml, one service per diversity collection.

A Traefik proxy routes to each service by its divbrowse.host and divbrowse.path labels.
A collection that can't supply VCF_URL, GFF3_URL and CHROM_PATTERN unambiguously fails
the build with its reasons.
"""

import re

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
# Routing: the `proxy` service (Traefik) listens on PROXY_PORT (default 80)
# and routes http://<divbrowse.host>/<divbrowse.path>/ to each Divbrowse
# service, based on the two labels set on it. To add a dataset, add a service
# using the x-divbrowse template with those labels, then `docker compose up -d`.
# Per-host root redirects live in traefik/dynamic/divbrowse.yml.

x-divbrowse: &divbrowse
  build:
    context: .
    dockerfile: Dockerfile
  restart: unless-stopped

services:
  proxy:
    image: traefik:v3
    ports:
      - "${PROXY_PORT:-80}:80"
    volumes:
      # Rootless Docker: set DOCKER_SOCK=$XDG_RUNTIME_DIR/docker.sock
      - ${DOCKER_SOCK:-/var/run/docker.sock}:/var/run/docker.sock:ro
      - ./traefik:/etc/traefik:ro
    restart: unless-stopped

"""

SERVICE = """  {name}:
    <<: *divbrowse
    environment:
      - VCF_URL={vcf_url}
      - GFF3_URL={gff3_url}
      - CHROM_PATTERN={chrom_pattern}
    labels:
      divbrowse.host: {host}
      divbrowse.path: {path}
    volumes:
      - ./data/{data_dir}:/opt/divbrowse
"""

# Public hostname per genus; --host adds to or overrides these.
DEFAULT_HOSTS = {
    "Glycine": "divbrowse.soybase.org",
    "Arachis": "divbrowse.peanutbase.org",
}

# Characters a compose service name, a data directory and an unquoted YAML value can
# all carry. Every diversity identifier in the store fits.
IDENTIFIER = re.compile(r"[A-Za-z0-9._-]+")
# One dot-separated name token: a species abbreviation, strain or assembly version.
# Strains can contain hyphens (IT97K-499-35), which are literal in CHROM_PATTERN.
TOKEN = re.compile(r"[A-Za-z0-9_-]+")
# A chromosome prefix. Some READMEs list chromosome names instead ("chr,Pt,Mt"), which
# cannot be turned into a pattern.
PREFIX = re.compile(r"[A-Za-z0-9_]+")
# A hostname both a Traefik Host() rule and an unquoted YAML value can carry.
HOSTNAME = re.compile(r"[A-Za-z0-9.-]+")


class DivbrowseError(Exception):
    """Raised when the requested collections can't all be expressed as services."""


def service_name(identifier):
    """Compose service name; lowercase, since compose rejects uppercase image names."""
    return f"divbrowse-{identifier.lower()}"


def data_directory(identifier):
    """Host directory for a service's data, named for its collection."""
    return identifier


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
    vcfs = [
        r.relative_path
        for r in collection.data_files
        if r.relative_path.endswith(".vcf.gz")
    ]
    if len(vcfs) != 1:
        listed = f": {', '.join(vcfs)}" if vcfs else ""
        problems.append(f"it has {len(vcfs)} VCF files{listed}")

    genome = None
    if len(collection.derived_from) == 1:
        genome = index.collections.get(
            f"{collection.genus}/{collection.species}/genomes/{collection.derived_from[0]}"
        )
    if genome is None:
        problems.append("it is not linked to a genome assembly")
        return None, problems

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
        return None, problems

    strain, gnm = tokens
    return {
        "vcf_url": f"{index.datastore_url}/{collection.path}/{vcfs[0]}",
        "gff3_url": f"{index.datastore_url}/{annotations[0].path}/{gff3[0]}",
        "chrom_pattern": f"{abbrev}.{strain}.{gnm}.{prefix}".replace(".", r"\.")
        + "[0-9]+",
        "host": host,
    }, []


def compose_file(index, identifiers, hosts=None):
    """The docker-compose.yml text, one service per identifier in the order given.

    Raises DivbrowseError listing every problem with every collection.
    """
    problems = {}
    blocks = []
    seen = set()
    for identifier in identifiers:
        if identifier in seen:
            problems.setdefault(identifier, []).append(
                "it was requested more than once"
            )
            continue
        seen.add(identifier)
        environment, reasons = service_environment(
            index, identifier, DEFAULT_HOSTS if hosts is None else hosts
        )
        if reasons:
            problems.setdefault(identifier, []).extend(reasons)
            continue
        blocks.append(
            SERVICE.format(
                name=service_name(identifier),
                path=identifier,
                data_dir=data_directory(identifier),
                **environment,
            )
        )
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
