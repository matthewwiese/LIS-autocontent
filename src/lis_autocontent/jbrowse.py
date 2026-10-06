"""Plan LIS's JBrowse 2 configuration from the datastore index, as jbrowse CLI commands.

The plan keeps every name, URL and grouping populate-jbrowse2 has always produced, so a
rebuilt config keeps its existing assembly names and tracks.
"""

import functools
import os
from dataclasses import dataclass

import yaml

from .datastore_files import SRC_CHECKSUM

CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "jbrowse.yml")

# Collection types in command order: assemblies before the tracks placed on them.
PLAN_TYPES = ("genomes", "annotations", "expression", "genome_alignments")


class JBrowseError(Exception):
    """Raised when jbrowse.yml can't be used."""


@functools.lru_cache(maxsize=None)
def load_config(path=CONFIG_PATH):
    """jbrowse.yml's instances and Jekyll link instance; JBrowseError if unusable."""
    try:
        with open(path, encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as err:
        raise JBrowseError(f"cannot read {path}: {err}") from err
    config = config if isinstance(config, dict) else {}
    instances = config.get("instances")
    if not (
        isinstance(instances, dict)
        and instances
        and all(isinstance(k, str) and isinstance(v, str) for k, v in instances.items())
    ):
        raise JBrowseError(f"{path}: instances must map each id to a URL")
    if config.get("jekyll_instance") not in instances:
        raise JBrowseError(f"{path}: jekyll_instance must name one of the instances")
    return config


def jekyll_url():
    """Base URL of the instance the Jekyll site's resource links open."""
    config = load_config()
    return config["instances"][config["jekyll_instance"]]


@dataclass
class JBrowseEntry:
    """One assembly or track. ``key`` is the name entries are deduplicated on.

    ``record`` is the index's file record for ``url``, or None when the URL is built
    by convention and the collection lists no such file.
    """

    kind: str  # assembly, annotation, expression or alignment
    key: str
    url: str
    assemblies: list
    label: str = ""
    category: str = ""
    record: object = None
    track_id: str = ""
    index_url: str = ""


def gensp(genus, species):
    """The prefix populate-jbrowse2 names files with: Glycine max -> glyma."""
    return f"{genus[:3].lower()}{species[:2].lower()}"


def lookup(collection):
    """<gensp>.<collection key without its trailing token>, e.g. glyma.Wm82.gnm4."""
    key = ".".join(collection.collection_key.split(".")[:-1])
    return f"{gensp(collection.genus, collection.species)}.{key}"


def genera(index, taxa_list=None):
    """Genera to plan: a taxon list's, in order, or every genus with a description."""
    if taxa_list and os.path.isfile(taxa_list):
        with open(taxa_list, encoding="utf-8") as handle:
            return [str(taxon["genus"]) for taxon in yaml.safe_load(handle) or []]
    return sorted(
        path.name
        for path in index.root.iterdir()
        if path.is_dir()
        and not path.name.startswith(".")
        and (
            path / "GENUS" / "about_this_collection" / f"description_{path.name}.yml"
        ).is_file()
    )


def _record(collection, filename):
    return next((r for r in collection.files if r.relative_path == filename), None)


def _csi_index(url, record):
    """``url``'s CSI index URL when the collection publishes CSI only, else ""."""
    if record and ".csi" in record.indexes and ".tbi" not in record.indexes:
        return f"{url}.csi"
    return ""


def _checksum_files(collection, suffix):
    return [
        r
        for r in collection.files
        if r.src == SRC_CHECKSUM and r.relative_path.endswith(suffix)
    ]


def collection_entries(index, collection):
    """The JBrowse entries one collection contributes."""
    name = lookup(collection)
    strain, version = name.split(".")[1], name.rsplit(".", maxsplit=1)[-1]
    title = f"{collection.genus.capitalize()} {collection.species} {strain}"
    base = f"{index.datastore_url}/{collection.path}"
    prefix = (
        f"{gensp(collection.genus, collection.species)}.{collection.collection_key}"
    )
    if collection.collection_type == "genomes":
        filename = f"{prefix}.genome_main.fna.gz"
        return [
            JBrowseEntry(
                "assembly",
                name,
                f"{base}/{filename}",
                [name],
                label=f"{title} V{version.replace('gnm', '')} Genomes",
                record=_record(collection, filename),
            )
        ]
    if collection.collection_type == "annotations":
        filename = f"{prefix}.gene_models_main.gff3.gz"
        record = _record(collection, filename)
        return [
            JBrowseEntry(
                "annotation",
                name,
                f"{base}/{filename}",
                [".".join(name.split(".")[:-1])],
                label=f"{title} V{version.replace('ann', '')} Annotations",
                record=record,
                # The CLI's default ID: the file name without its last extension.
                track_id=os.path.splitext(filename)[0],
                index_url=_csi_index(f"{base}/{filename}", record),
            )
        ]
    if collection.collection_type == "expression":
        parent = ".".join(name.split(".")[:-3])
        return [
            JBrowseEntry(
                "expression",
                r.relative_path,
                f"{base}/{r.relative_path}",
                [parent],
                label=r.relative_path.split(".")[-2],
                category=".".join(r.relative_path.split(".")[1:-2]),
                record=r,
            )
            for r in _checksum_files(collection, "bw")
        ]
    if collection.collection_type == "genome_alignments":
        entries = []
        for r in _checksum_files(collection, "paf.gz"):
            parts = r.relative_path.split(".")
            first, second = ".".join(parts[:3]), ".".join(parts[4:7])
            entries.append(
                JBrowseEntry(
                    "alignment",
                    r.relative_path,
                    f"{base}/{r.relative_path}",
                    [second, first],
                    record=r,
                )
            )
        return entries
    return []


def plan(index, genus_names):
    """Every JBrowse entry for these genera, in command order.

    Within a type, entries follow genus (as given), species (as its genus description
    lists them) and collection path. An entry whose key is already planned replaces
    the earlier one in place.
    """
    grouped = {}
    for collection in index.collections.values():
        group = (collection.genus, collection.species, collection.collection_type)
        grouped.setdefault(group, []).append(collection)
    planned = {ctype: {} for ctype in PLAN_TYPES}
    for genus in genus_names:
        for species in index.taxa.get(genus, {}).get("species_in_genus", []):
            for ctype in PLAN_TYPES:
                for collection in grouped.get((genus, species, ctype), []):
                    for entry in collection_entries(index, collection):
                        planned[ctype][entry.key] = entry
    return [entry for ctype in PLAN_TYPES for entry in planned[ctype].values()]


def command(entry, out_dir):
    """The jbrowse CLI command that adds ``entry`` to the config in ``out_dir``."""
    out = os.path.abspath(out_dir)
    if entry.kind == "assembly":
        return (
            f"jbrowse add-assembly -n {entry.key} --out {out}/ -t bgzipFasta --force"
            f' --displayName "{entry.label}" {entry.url}'
        )
    if entry.kind == "annotation":
        index = f" --indexFile {entry.index_url}" if entry.index_url else ""
        return (
            f"jbrowse add-track -a {entry.assemblies[0]} --out {out}/ --force"
            f' -n "{entry.label}" --trackId {entry.track_id}{index} {entry.url}'
        )
    if entry.kind == "expression":
        return (
            f"jbrowse add-track {entry.url} --name {entry.label}"
            f" --assemblyNames {entry.assemblies[0]}"
            f" --category expression,{entry.category} --out {out} --force"
        )
    bam_url = entry.url.replace("paf.gz", "bam")
    bam_name = os.path.basename(bam_url)
    return (
        f"jbrowse add-track --assemblyNames {','.join(entry.assemblies)}"
        f" --out {out}/ {entry.url} --force"
        f";jbrowse add-track -n {bam_name} --trackId {bam_name}"
        f" -a {entry.assemblies[1]}"
        f" --out {out}/ --indexFile {bam_url}.bai {bam_url} --force"
    )
