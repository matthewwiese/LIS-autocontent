"""Build a Divbrowse docker-compose.yml for chosen diversity collections.

legumeinfo/divbrowse ships a docker-compose.yml with one service per diversity dataset.
Each service downloads a VCF and a GFF3, keeps the chromosomes matching CHROM_PATTERN,
and serves them under BASE_URL. All four values come from datastore metadata:

    VCF_URL        the diversity collection's one VCF
    GFF3_URL       gene_models_main of the one annotation of the linked genome
    CHROM_PATTERN  <abbrev>.<strain>.<gnm>.<chromosome_prefix>[0-9]+, from that genome
    BASE_URL       <base_url>/<collection>/

The file is written in the upstream layout. A collection that can't be expressed that
way -- several VCFs or none, no linked genome, several annotations, or a
chromosome_prefix that isn't a single prefix -- stops the build with every reason
reported, rather than producing a service built on a guess.
"""

import re

HEADER = """# Divbrowse Docker Compose
#
# Just run: docker compose up -d
#
# On first start, each service will automatically:
#   1. Download VCF and GFF3 from the specified URLs
#   2. Convert VCF to Zarr format
#   3. Generate configuration
#   4. Start the server
#
# Data is persisted in ./data/<name>/ directories.

services:
"""

SERVICE = """  {name}:
    build:
      context: .
      dockerfile: Dockerfile
    environment:
      - VCF_URL={vcf_url}
      - GFF3_URL={gff3_url}
      - CHROM_PATTERN={chrom_pattern}
      - BASE_URL={base_url}
    ports:
      - "{port}:8080"
    volumes:
      - ./data/{data_dir}:/opt/divbrowse
    restart: unless-stopped
"""

# Characters a compose service name, a data directory and an unquoted YAML value can
# all carry. Every diversity identifier in the store fits.
IDENTIFIER = re.compile(r"[A-Za-z0-9._-]+")
# One dot-separated name token: a species abbreviation, strain or assembly version.
# Strains can contain hyphens (IT97K-499-35), which are literal in CHROM_PATTERN.
TOKEN = re.compile(r"[A-Za-z0-9_-]+")
# A chromosome prefix. Some READMEs list chromosome names instead ("chr,Pt,Mt"), which
# cannot be turned into a pattern.
PREFIX = re.compile(r"[A-Za-z0-9_]+")
MAX_PORT = 65535


class DivbrowseError(Exception):
    """Raised when the requested collections can't all be expressed as services."""


def service_name(identifier):
    """Compose names each service's image after the service, and image names must be
    lowercase: `docker compose build` rejects divbrowse-Wm82.gnm4.div.Song_Hyten_2015.
    """
    return f"divbrowse-{identifier.lower()}"


def data_directory(identifier):
    """Host directory for a service's data, named for its collection."""
    return identifier


def service_environment(index, identifier, base_url):
    """VCF_URL, GFF3_URL, CHROM_PATTERN and BASE_URL for one diversity collection.

    Returns (environment, problems). The environment is None whenever there are
    problems, so no service is ever written from a partial or guessed value.
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
        "base_url": f"{base_url.rstrip('/')}/{identifier}/",
    }, []


def compose_file(index, identifiers, base_url, first_port):
    """The docker-compose.yml text for these collections, as services in the order given.

    Raises DivbrowseError listing every problem with every collection when any of them
    can't be expressed; nothing is returned in that case.
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
        environment, reasons = service_environment(index, identifier, base_url)
        if reasons:
            problems.setdefault(identifier, []).extend(reasons)
            continue
        blocks.append(
            SERVICE.format(
                name=service_name(identifier),
                port=first_port + len(blocks),
                data_dir=data_directory(identifier),
                **environment,
            )
        )
    if not identifiers:
        problems[""] = ["no collections were requested"]
    last_port = first_port + len(identifiers) - 1
    if identifiers and not 1 <= first_port <= last_port <= MAX_PORT:
        problems.setdefault("", []).append(
            f"ports {first_port}-{last_port} fall outside 1-{MAX_PORT}"
        )
    if problems:
        raise DivbrowseError(
            "\n".join(
                f"{identifier or 'request'}: {'; '.join(reasons)}"
                for identifier, reasons in problems.items()
            )
        )
    return HEADER + "\n".join(blocks)
