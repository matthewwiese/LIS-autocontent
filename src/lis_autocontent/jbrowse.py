"""Plan LIS's JBrowse 2 configuration from the datastore index, as jbrowse CLI commands.

The plan keeps every name, URL and grouping populate-jbrowse2 has always produced, so a
rebuilt config keeps its existing assembly names and tracks.
"""

import functools
import hashlib
import json
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
    for key in ("built_instance", "jekyll_instance"):
        if config.get(key) not in instances:
            raise JBrowseError(f"{path}: {key} must name one of the instances")
    return config


def jekyll_url():
    """Base URL of the instance the Jekyll site's resource links open."""
    config = load_config()
    return config["instances"][config["jekyll_instance"]]


@dataclass
class JBrowseEntry:  # pylint: disable=too-many-instance-attributes
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


def _named_entry(index, collection):
    """The assembly or annotation track a genome or annotation collection names."""
    name = lookup(collection)
    strain, version = name.split(".")[1], name.rsplit(".", maxsplit=1)[-1]
    title = f"{collection.genus.capitalize()} {collection.species} {strain}"
    prefix = (
        f"{gensp(collection.genus, collection.species)}.{collection.collection_key}"
    )
    base = f"{index.datastore_url}/{collection.path}"
    if collection.collection_type == "genomes":
        filename = f"{prefix}.genome_main.fna.gz"
        return JBrowseEntry(
            "assembly",
            name,
            f"{base}/{filename}",
            [name],
            label=f"{title} V{version.replace('gnm', '')} Genomes",
            record=_record(collection, filename),
        )
    filename = f"{prefix}.gene_models_main.gff3.gz"
    record = _record(collection, filename)
    return JBrowseEntry(
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


def _expression_entries(index, collection):
    parent = ".".join(lookup(collection).split(".")[:-3])
    return [
        JBrowseEntry(
            "expression",
            r.relative_path,
            f"{index.datastore_url}/{collection.path}/{r.relative_path}",
            [parent],
            label=r.relative_path.split(".")[-2],
            category=".".join(r.relative_path.split(".")[1:-2]),
            record=r,
        )
        for r in _checksum_files(collection, "bw")
    ]


def _alignment_entries(index, collection):
    entries = []
    for r in _checksum_files(collection, "paf.gz"):
        parts = r.relative_path.split(".")
        entries.append(
            JBrowseEntry(
                "alignment",
                r.relative_path,
                f"{index.datastore_url}/{collection.path}/{r.relative_path}",
                [".".join(parts[4:7]), ".".join(parts[:3])],
                record=r,
            )
        )
    return entries


def collection_entries(index, collection):
    """The JBrowse entries one collection contributes."""
    ctype = collection.collection_type
    if ctype in ("genomes", "annotations"):
        return [_named_entry(index, collection)]
    if ctype == "expression":
        return _expression_entries(index, collection)
    if ctype == "genome_alignments":
        return _alignment_entries(index, collection)
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


# Index siblings a config may name, with the indexType JBrowse records for each.
INDEX_TYPES = {".tbi": "TBI", ".csi": "CSI", ".bai": "BAI", ".crai": "CRAI"}
# Siblings that never name an assembly's or track's data file.
SIBLING_SUFFIXES = (*INDEX_TYPES, ".fai", ".gzi")


@dataclass
class Deployment:  # pylint: disable=too-many-instance-attributes
    """One instance's deployed config, mapped onto the index's collections.

    ``placements`` maps collection path -> {"assemblies": set, "tracks": list}. The
    remaining lists name what could not be placed, for the build report.
    """

    instance: str
    url: str
    status: str = "unavailable"
    detail: str = ""
    sha256: str = ""
    placements: dict = None
    assemblies: dict = None  # assembly name -> collection path
    tracks: dict = None  # track id -> index type ("" when none)
    elsewhere: list = None  # ids served from relative or non-Data Store URLs
    unmatched: list = None  # (id, url) on the Data Store but in no known collection
    unlisted: list = None  # (id, url) for files their collection doesn't list
    unpublished: list = None  # (id, url) for index files their collection doesn't list


def _uris(node):
    """Every ``uri`` value in a JBrowse adapter, in document order."""
    if isinstance(node, dict):
        for key, value in node.items():
            if key == "uri" and isinstance(value, str):
                yield value
            else:
                yield from _uris(value)
    elif isinstance(node, list):
        for item in node:
            yield from _uris(item)


def locate(index, url):
    """(collection, relative path) for a Data Store URL, or (None, None)."""
    prefix = f"{index.datastore_url}/"
    parts = url[len(prefix) :].split("/") if url.startswith(prefix) else []
    if len(parts) < 5:
        return None, None
    return index.collections.get("/".join(parts[:4])), "/".join(parts[4:])


def _listed(collection, relative_path):
    """False only when the collection's CHECKSUM is known and lacks the file."""
    return collection.index_status != "known" or any(
        r.relative_path == relative_path for r in collection.files
    )


def _place(deployment, index, item_id, adapter):
    """The collection path an adapter's data file lives in, noting what doesn't fit."""
    uris = list(_uris(adapter))
    data = next((u for u in uris if not u.endswith(SIBLING_SUFFIXES)), None)
    collection, relative = locate(index, data) if data else (None, None)
    if relative is None:
        deployment.elsewhere.append(item_id)
        return None, None
    if collection is None:
        deployment.unmatched.append((item_id, data))
        return None, None
    if not _listed(collection, relative):
        deployment.unlisted.append((item_id, data))
    index_url = next((u for u in uris if u.endswith(tuple(INDEX_TYPES))), None)
    if index_url:
        index_collection, index_relative = locate(index, index_url)
        if index_collection is not None and not _listed(
            index_collection, index_relative
        ):
            deployment.unpublished.append((item_id, index_url))
    return collection.path, relative


def _placement(deployment, path):
    return deployment.placements.setdefault(path, {"assemblies": set(), "tracks": []})


def _add_assembly(deployment, index, assembly):
    name = assembly.get("name", "")
    path, _ = _place(deployment, index, name, assembly.get("sequence", {}))
    if path:
        deployment.assemblies[name] = path
        _placement(deployment, path)["assemblies"].add(name)


def _add_track(deployment, index, track):
    track_id, adapter = track.get("trackId", ""), track.get("adapter", {})
    index_spec = adapter.get("index", {}) if isinstance(adapter, dict) else {}
    index_url = next(iter(_uris(index_spec)), "")
    index_type = index_spec.get("indexType") or next(
        (t for s, t in INDEX_TYPES.items() if index_url.endswith(s)), ""
    )
    deployment.tracks[track_id] = index_type
    path, relative = _place(deployment, index, track_id, adapter)
    if not path:
        return
    placement = _placement(deployment, path)
    placement["assemblies"].update(track.get("assemblyNames") or [])
    entry = {"id": track_id, "type": track.get("type", ""), "file": relative}
    if index_type:
        entry["index"] = index_type
    placement["tracks"].append(entry)


def read_deployment(index, instance, url, config_path):
    """Read one instance's config.json; an unreadable one leaves it unavailable."""
    deployment = Deployment(instance, url)
    if not config_path:
        deployment.detail = "no config supplied"
        return deployment
    try:
        with open(config_path, "rb") as handle:
            raw = handle.read()
        config = json.loads(raw)
        assemblies, tracks = config["assemblies"], config["tracks"]
    except (OSError, ValueError, KeyError, TypeError) as err:
        deployment.detail = f"unreadable config: {err}"
        return deployment
    deployment.status, deployment.sha256 = "ok", hashlib.sha256(raw).hexdigest()
    deployment.placements, deployment.assemblies, deployment.tracks = {}, {}, {}
    deployment.elsewhere, deployment.unmatched = [], []
    deployment.unlisted, deployment.unpublished = [], []
    for assembly in assemblies:
        _add_assembly(deployment, index, assembly)
    for track in tracks:
        _add_track(deployment, index, track)
    return deployment


def read_deployments(index, config_paths):
    """A Deployment for every instance in jbrowse.yml; config_paths maps id -> path."""
    instances = load_config()["instances"]
    unknown = sorted(set(config_paths) - set(instances))
    if unknown:
        raise JBrowseError(f"not instances in jbrowse.yml: {', '.join(unknown)}")
    return [
        read_deployment(index, instance, url, config_paths.get(instance))
        for instance, url in instances.items()
    ]


def catalog_section(deployments):
    """(jbrowse_instances, {collection path: placements}) for the catalog."""
    instances, placements = {}, {}
    for deployment in deployments:
        record = {"url": f"{deployment.url}/", "status": deployment.status}
        if deployment.status == "ok":
            record["config_sha256"] = deployment.sha256
        else:
            record["detail"] = deployment.detail
        instances[deployment.instance] = record
        for path, placed in sorted((deployment.placements or {}).items()):
            entry = {
                "instance": deployment.instance,
                "assemblies": sorted(placed["assemblies"]),
            }
            if placed["tracks"]:
                entry["tracks"] = sorted(placed["tracks"], key=lambda t: t["id"])
            placements.setdefault(path, []).append(entry)
    return instances, placements


def planned_track_ids(entry):
    """The track IDs an entry's commands create, by the jbrowse CLI's rules."""
    name = os.path.basename(entry.url)
    if entry.kind == "annotation":
        return [entry.track_id]
    if entry.kind == "expression":
        return [os.path.splitext(name)[0]]
    if entry.kind == "alignment":
        return [os.path.splitext(name)[0], name.replace("paf.gz", "bam")]
    return []


def _items(rows):
    return [f"- `{item_id}`: {url}" for item_id, url in rows]


def _plan_drift(deployment, entries):
    """Report lines comparing the plan with the instance populate-jbrowse2 builds."""
    planned_assemblies = {e.key for e in entries if e.kind == "assembly"}
    planned_tracks = {
        track_id: ("CSI" if e.index_url else "TBI") if e.kind == "annotation" else ""
        for e in entries
        for track_id in planned_track_ids(e)
    }
    rows = [
        (
            "Assemblies planned, not deployed",
            planned_assemblies - set(deployment.assemblies),
        ),
        (
            "Assemblies deployed, not planned",
            set(deployment.assemblies) - planned_assemblies,
        ),
        ("Tracks planned, not deployed", set(planned_tracks) - set(deployment.tracks)),
        ("Tracks deployed, not planned", set(deployment.tracks) - set(planned_tracks)),
        (
            "Tracks whose deployed index type differs from the plan's",
            {
                t
                for t, kind in planned_tracks.items()
                if kind and deployment.tracks.get(t, kind) != kind
            },
        ),
    ]
    lines = [f"### Plan vs deployed: {deployment.instance}", ""]
    for label, names in rows:
        lines.append(f"- {label}: {len(names)}")
        lines.extend(f"  - `{name}`" for name in sorted(names))
    return lines


def _instance_table(deployments):
    lines = [
        "## JBrowse instances",
        "",
        "| Instance | Status | Assemblies | Tracks placed | Elsewhere "
        "| Unknown collection | File unlisted | Index unpublished |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for d in deployments:
        if d.status != "ok":
            lines.append(f"| {d.instance} | {d.status}: {d.detail} | | | | | | |")
            continue
        placed = sum(len(p["tracks"]) for p in d.placements.values())
        lines.append(
            f"| {d.instance} | ok | {len(d.assemblies)} | {placed} of {len(d.tracks)} "
            f"| {len(d.elsewhere)} | {len(d.unmatched)} | {len(d.unlisted)} "
            f"| {len(d.unpublished)} |"
        )
    return lines


def _problem_lines(deployment):
    lines = []
    for label, rows in (
        ("on the Data Store, in no known collection", deployment.unmatched),
        ("naming a file its collection doesn't list", deployment.unlisted),
        ("naming an index file its collection doesn't publish", deployment.unpublished),
    ):
        if rows:
            lines += [
                "",
                f"### {deployment.instance}: tracks {label}",
                "",
                *_items(rows),
            ]
    return lines


def _unserved_lines(index, deployments, ok):
    served = {path for d in ok for path in d.placements}
    unserved = sorted(
        c.path
        for c in index.collections.values()
        if c.collection_type in ("genomes", "annotations") and c.path not in served
    )
    caveat = "" if len(ok) == len(deployments) else " (some instances unavailable)"
    lines = ["", f"### Genome and annotation collections no instance serves{caveat}"]
    return lines + (
        ["", *(f"- {path}" for path in unserved)] if unserved else ["", "None."]
    )


def report(index, deployments, entries, built_instance):
    """A Markdown summary of what each instance serves and where it drifts."""
    ok = [d for d in deployments if d.status == "ok"]
    lines = _instance_table(deployments)
    for deployment in ok:
        lines += _problem_lines(deployment)
    lines += _unserved_lines(index, deployments, ok)
    built = next((d for d in ok if d.instance == built_instance), None)
    if built:
        lines += ["", *_plan_drift(built, entries)]
    return "\n".join(lines) + "\n"
