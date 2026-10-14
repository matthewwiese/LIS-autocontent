"""Plan LIS's Genome Context Viewer instances as gcv-docker-compose's gcvctl runs them.

Each instance in gcv.yml becomes instances/<name>.env and a load script in
data/<name>/docker-entrypoint-initdb.d/, which redis_loader runs while the instance's
Redis is empty. Copied into a gcv-docker-compose checkout, `./gcvctl load <name>` builds
an instance's database and `./gcvctl up <name>` serves it behind the shared proxy at
/<name>/. Any member that can't be loaded unambiguously from published files fails the
build with its reasons.
"""

import functools
import os
import re
import shutil
from dataclasses import dataclass, field

import yaml

from .datastore_files import STATUS_KNOWN

CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gcv.yml")
INPUTS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gcv_inputs.py")

TEMPLATE_URL = "https://github.com/legumeinfo/gcv-docker-compose/tree/{tag}"
# gcvctl's defaults for an instance's files.
ENV_DIR = "instances"
INITDB = "data/{name}/docker-entrypoint-initdb.d"
LOAD_SCRIPT = "10-load.sh"
INPUTS_SCRIPT = "gcv_inputs.py"
PROFILES = ("services", "all")
# .env keys populate-gcv or gcvctl sets, and RDB_URL/RDB_SHA256, which a release of
# built databases writes.
DERIVED_ENV = (
    "COMPOSE_PROJECT_NAME",
    "COMPOSE_FILE",
    "COMPOSE_PROFILES",
    "MICROSERVICES_PATH",
    "REDIS_DATA_DIR",
    "INITDB_DIR",
    "RDB_URL",
    "RDB_SHA256",
)

# Compose project names are lowercase.
INSTANCE_NAME = re.compile(r"[a-z0-9][a-z0-9_-]*")
ENV_KEY = re.compile(r"[A-Z_][A-Z0-9_]*")
ENV_VALUE = re.compile(r"[^\s#'\"]+")
# Words written into the load script unquoted.
SHELL_WORD = re.compile(r"[A-Za-z0-9._:/+=-]+")
FAMILY = re.compile(r"[A-Za-z0-9_.]+")
# The GCV web app's path: absolute, with the trailing slash its nginx config needs.
GCV_PATH = re.compile(r"/(?:[A-Za-z0-9._-]+/)+")


class GcvError(Exception):
    """Raised when the requested instances can't all be built from published files."""


@dataclass
class GcvLoad:  # pylint: disable=too-many-instance-attributes
    """One redis_loader gff call: an annotation, its genome's sequences and families.

    ``rename`` is the (old, new) assembly token rewritten in every sequence name.
    """

    path: str
    genus: str
    species: str
    strain: str
    fai_url: str
    gene_gff_url: str
    gfa_url: str
    rename: tuple = None


@dataclass
class GcvInstance:
    """One gcv-docker-compose deployment and the annotations it loads, in order."""

    name: str
    profile: str
    family: str
    env: dict
    loads: list = field(default_factory=list)


def _member(entry):
    """(path, rename_chromosomes) from a members entry, a path or a mapping."""
    if isinstance(entry, str):
        return entry, False
    if isinstance(entry, dict) and isinstance(entry.get("path"), str):
        if set(entry) <= {"path", "rename_chromosomes"} and isinstance(
            entry.get("rename_chromosomes", False), bool
        ):
            return entry["path"], entry.get("rename_chromosomes", False)
    return None, False


def _env_problems(env, where):
    problems = []
    if not isinstance(env, dict):
        return [f"{where} env must be a mapping"]
    for key, value in env.items():
        if not (isinstance(key, str) and ENV_KEY.fullmatch(key)):
            problems.append(f"{where} env key {key!r} is not a variable name")
        elif key in DERIVED_ENV:
            problems.append(f"{where} env sets {key}, which populate-gcv derives")
        if not (
            isinstance(value, (str, int))
            and not isinstance(value, bool)
            and ENV_VALUE.fullmatch(str(value))
        ):
            problems.append(f"{where} env {key}={value!r} can't be written unquoted")
    return problems


def _instance_problems(name, instance, env):
    """Problems with one instance's definition, before any lookup in the index."""
    where = f"instance {name}"
    if not isinstance(instance, dict):
        return [f"{where} must be a mapping"]
    problems = []
    if not (isinstance(name, str) and INSTANCE_NAME.fullmatch(name)):
        problems.append(f"{where}: name must be lowercase letters, digits, - or _")
    unknown = set(instance) - {"profile", "family", "env", "members"}
    if unknown:
        problems.append(f"{where} has unknown keys {', '.join(sorted(unknown))}")
    if instance.get("profile") not in PROFILES:
        problems.append(f"{where} profile must be one of {', '.join(PROFILES)}")
    if not (
        isinstance(instance.get("family"), str) and FAMILY.fullmatch(instance["family"])
    ):
        problems.append(f"{where} needs a family, such as legfed_v1_0")
    own_env = instance.get("env", {})
    problems.extend(_env_problems(own_env, where))
    merged = {**env, **own_env} if isinstance(own_env, dict) else env
    if instance.get("profile") == "all" and not GCV_PATH.fullmatch(
        str(merged.get("GCV_PATH", ""))
    ):
        problems.append(
            f"{where} runs the web app but sets no GCV_PATH like /{name}/gcv/"
        )
    members = instance.get("members")
    if not (isinstance(members, list) and members):
        problems.append(f"{where} needs a list of members")
    else:
        for entry in members:
            if _member(entry)[0] is None:
                problems.append(
                    f"{where} member {entry!r} is not a path or "
                    "{path, rename_chromosomes}"
                )
    return problems


def _path_problems(instances, env):
    """Web app paths shared by two instances, which the shared proxy can't route."""
    owners = {}
    for name, instance in instances.items():
        if instance["profile"] == "all":
            path = {**env, **instance.get("env", {})}["GCV_PATH"]
            owners.setdefault(path, []).append(name)
    return [
        f"GCV_PATH {path} is set by {' and '.join(names)}"
        for path, names in sorted(owners.items())
        if len(names) > 1
    ]


@functools.lru_cache(maxsize=None)
def load_config(path=CONFIG_PATH):
    """gcv.yml, checked for shape; GcvError listing every problem."""
    try:
        with open(path, encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
    except (OSError, yaml.YAMLError) as err:
        raise GcvError(f"cannot read {path}: {err}") from err
    config = config if isinstance(config, dict) else {}
    problems = []
    if not (
        isinstance(config.get("template"), str)
        and SHELL_WORD.fullmatch(config["template"])
    ):
        problems.append("template must name a gcv-docker-compose release tag")
    env = config.get("env", {})
    problems.extend(_env_problems(env, "top-level"))
    instances = config.get("instances")
    if not (isinstance(instances, dict) and instances):
        problems.append("instances must map each instance name to its definition")
        instances = {}
    for name, instance in instances.items():
        problems.extend(_instance_problems(name, instance, env))
    if not problems:
        problems.extend(_path_problems(instances, env))
    if problems:
        raise GcvError(f"{path}: " + "\n".join(problems))
    return config


def _single(records, what):
    """(record, problems) for exactly one of ``records``."""
    if len(records) == 1:
        return records[0], []
    listed = f": {', '.join(r.relative_path for r in records)}" if records else ""
    return None, [f"it has {len(records)} {what}{listed}"]


def genome_fai(index, annotation):
    """(genome collection, .fai record, problems) for an annotation's assembly."""
    if len(annotation.derived_from) != 1:
        return None, None, ["it is not linked to a genome assembly"]
    genome = index.collections.get(
        f"{annotation.genus}/{annotation.species}/genomes/{annotation.derived_from[0]}"
    )
    if genome is None or not genome.files:
        return None, None, [f"its genome {annotation.derived_from[0]} has no CHECKSUM"]
    main, problems = _single(
        [
            r
            for r in genome.data_files
            if r.canonical_type == "genome_main"
            and r.extensions == "fna.gz"
            and "/" not in r.relative_path
        ],
        f"genome_main FASTA files in {genome.collection_key}",
    )
    if problems:
        return genome, None, problems
    fai = index.companion(main, ".fai")
    if fai is None:
        return (
            genome,
            None,
            [f"{genome.collection_key}'s CHECKSUM lists no {main.filename}.fai"],
        )
    return genome, fai, []


def family_file(annotation, gene_models, family):
    """(record, problems) for the annotation's one top-level family file of ``family``."""
    stem = gene_models.filename[: -len(".gene_models_main.gff3.gz")]
    pattern = re.compile(
        rf"{re.escape(stem)}\.{re.escape(family)}\.[A-Za-z0-9]+\.gfa\.tsv\.gz"
    )
    return _single(
        [r for r in annotation.data_files if pattern.fullmatch(r.relative_path)],
        f"{family} family files",
    )


def plan_load(index, path, family, rename_chromosomes):
    """(GcvLoad, problems) for one member annotation; the load is None on any problem."""
    annotation = index.collections.get(path)
    if annotation is None:
        return None, ["no such collection in the checkout"]
    if annotation.collection_type != "annotations":
        return None, ["it is not an annotation"]
    if not annotation.files or annotation.index_status != STATUS_KNOWN:
        return None, ["it has no CHECKSUM"]
    problems = []
    gene_models, found = _single(
        [
            r
            for r in annotation.data_files
            if r.canonical_type == "gene_models_main"
            and r.extensions == "gff3.gz"
            and "/" not in r.relative_path
        ],
        "gene_models_main GFF3 files",
    )
    problems.extend(found)
    gfa = None
    if gene_models is not None:
        gfa, found = family_file(annotation, gene_models, family)
        problems.extend(found)
    genome, fai, found = genome_fai(index, annotation)
    problems.extend(found)
    rename = None
    if rename_chromosomes and genome is not None:
        old = genome.collection_key.split(".")[1]
        rename = (old, ".".join(annotation.collection_key.split(".")[1:3]))
    load = GcvLoad(
        path=path,
        genus=annotation.genus,
        species=annotation.species,
        strain=annotation.collection_key.split(".")[0],
        fai_url=fai.url(index.datastore_url) if fai else None,
        gene_gff_url=gene_models.url(index.datastore_url) if gene_models else None,
        gfa_url=gfa.url(index.datastore_url) if gfa else None,
        rename=rename,
    )
    problems.extend(shell_word_problems(load))
    return (None if problems else load), problems


def shell_word_problems(load):
    """Values of a load that the load script can't carry unquoted."""
    words = (load.genus, load.species, load.strain, *(load.rename or ()))
    urls = (load.fai_url, load.gene_gff_url, load.gfa_url)
    return [
        f"{word!r} can't be written into a shell script unquoted"
        for word in (*words, *(url for url in urls if url))
        if not SHELL_WORD.fullmatch(word)
    ]


def collision_problems(index, members):
    """Members sharing an assembly whose chromosome names would collide."""
    by_genome = {}
    for path, rename_chromosomes in members:
        annotation = index.collections.get(path)
        if annotation is not None and len(annotation.derived_from) == 1:
            genome = (annotation.genus, annotation.species, annotation.derived_from[0])
            by_genome.setdefault(genome, []).append((path, rename_chromosomes))
    problems = []
    for (_, _, genome), sharing in sorted(by_genome.items()):
        if sum(not rename for _, rename in sharing) > 1:
            paths = ", ".join(path for path, _ in sharing)
            problems.append(
                f"{paths} share assembly {genome}; set rename_chromosomes on all but one"
            )
    return problems


def plan(index, config, names=None):
    """GcvInstance records for ``names`` (default: all), in gcv.yml's order.

    Raises GcvError listing every problem with every member.
    """
    instances = config["instances"]
    names = list(instances) if not names else list(names)
    problems = [
        f"no instance {name} in gcv.yml" for name in names if name not in instances
    ]
    planned = []
    for name in (n for n in instances if n in names):
        definition = instances[name]
        members = [_member(entry) for entry in definition["members"]]
        instance = GcvInstance(
            name=name,
            profile=definition["profile"],
            family=definition["family"],
            env={**config.get("env", {}), **definition.get("env", {})},
        )
        seen = set()
        for path, rename_chromosomes in members:
            if path in seen:
                problems.append(f"{name}: {path}: it is listed more than once")
                continue
            seen.add(path)
            load, reasons = plan_load(
                index, path, definition["family"], rename_chromosomes
            )
            problems.extend(f"{name}: {path}: {reason}" for reason in reasons)
            if load is not None:
                instance.loads.append(load)
        problems.extend(f"{name}: {p}" for p in collision_problems(index, members))
        planned.append(instance)
    if problems:
        raise GcvError("\n".join(problems))
    return planned


def env_file(instance):
    """instances/<name>.env, which gcvctl passes to docker compose for the instance."""
    env = {
        "COMPOSE_PROFILES": instance.profile,
        "MICROSERVICES_PATH": f"/{instance.name}",
        **instance.env,
    }
    lines = "".join(f"{key}={value}\n" for key, value in env.items())
    return f"# {instance.name}: generated by lis-autocontent populate-gcv\n{lines}"


def load_command(load):
    """Shell lines that load one annotation."""
    rename = f" \\\n  --rename {' '.join(load.rename)}" if load.rename else ""
    lines = [
        f"echo {load.path}",
        f'python -u "$inputs" chromosomes \\\n  {load.fai_url} \\\n'
        f'  "$work/chromosomes.gff3"{rename}',
    ]
    gene_gff = load.gene_gff_url
    if load.rename:
        lines.append(
            f'python -u "$inputs" genes \\\n  {load.gene_gff_url} \\\n'
            f'  "$work/genes.gff3"{rename}'
        )
        gene_gff = '"$work/genes.gff3"'
    lines.append(
        "python -u -m redis_loader --load-type append gff \\\n"
        f"  --genus {load.genus} \\\n"
        f"  --species {load.species} \\\n"
        f"  --strain {load.strain} \\\n"
        '  --chromosome-gff "$work/chromosomes.gff3" \\\n'
        f"  --gene-gff {gene_gff} \\\n"
        f"  --gfa {load.gfa_url}"
    )
    return "\n".join(lines) + "\n"


def load_script(instance):
    """10-load.sh, run by redis_loader's entrypoint while Redis is empty."""
    header = f"""#!/bin/sh
# {instance.name}: {len(instance.loads)} annotations, {instance.family} families.
# Generated by lis-autocontent populate-gcv. redis_loader runs this only while Redis
# is empty, saving after each load.

set -o errexit -o pipefail

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
inputs=/docker-entrypoint-initdb.d/{INPUTS_SCRIPT}
"""
    return header + "".join(f"\n{load_command(load)}" for load in instance.loads)


def readme(instances, config):
    """README.md: how to use the tree, and what each instance loads."""
    rows = "".join(
        f"| `{i.name}` | {i.profile} | `/{i.name}/` | `{i.family}` | {len(i.loads)} |\n"
        for i in instances
    )
    tag = config["template"]
    return f"""# LIS GCV instances

Generated by `lis-autocontent populate-gcv` for gcv-docker-compose
[`{tag}`]({TEMPLATE_URL.format(tag=tag)}). Copy `instances/` and `data/` into a checkout
of it, then:

```console
./gcvctl proxy              # once: the shared proxy
./gcvctl load <instance>    # build the instance's database from the Data Store
./gcvctl up <instance>      # serve it at /<instance>/
```

To serve a database built elsewhere instead of loading one, set `RDB_URL` and
`RDB_SHA256` in the instance's `.env` before `./gcvctl up`.

| Instance | Profile | Path | Families | Annotations |
| --- | --- | --- | --- | --- |
{rows}"""


def render(instances, config, out_dir):
    """Write every instance's .env and load scripts, and README.md; returns the paths."""

    def write(relative, text, mode=0o644):
        path = os.path.join(out_dir, relative)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(text)
        os.chmod(path, mode)
        written.append(path)

    written = []
    write("README.md", readme(instances, config))
    for instance in instances:
        initdb = INITDB.format(name=instance.name)
        write(f"{ENV_DIR}/{instance.name}.env", env_file(instance))
        write(f"{initdb}/{LOAD_SCRIPT}", load_script(instance), 0o755)
        inputs = os.path.join(out_dir, initdb, INPUTS_SCRIPT)
        shutil.copyfile(INPUTS_PATH, inputs)
        written.append(inputs)
    return written
