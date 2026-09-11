"""End-to-end tests of the CLI subcommands and the artifacts they hand downstream.

Each subcommand runs through click against a small datastore-metadata clone with the
network stubbed out. That catches wiring mistakes in lis_cli.py -- an option missing
from a function's signature, arguments passed in the wrong order -- as well as changes
to what consumers receive: DSCensor node files, BLAST and JBrowse2 commands, and the
Jekyll YAML.

The expected artifacts in tests/data/cli/ were reviewed by hand. When a change to an
artifact is intended, regenerate them with `LIS_UPDATE_EXPECTED=1 pytest
tests/test_cli.py` and review the diff before committing it: regenerating without
reading the diff defeats the test.
"""

import json
import os
import traceback
from types import SimpleNamespace

import pytest
import requests
import yaml
from click.testing import CliRunner

from lis_autocontent import lis_cli
from lis_autocontent.catalog import CatalogBuilder
from lis_autocontent.datastore_files import DatastoreIndex

EXPECTED = os.path.join(os.path.dirname(__file__), "data", "cli")
UPDATING = bool(os.environ.get("LIS_UPDATE_EXPECTED"))
MD5 = "d41d8cd98f00b204e9800998ecf8427e"

GENOME = "Cicer/arietinum/genomes/CDCFrontier.gnm3.QT0P"
ANNOTATION = "Cicer/arietinum/annotations/CDCFrontier.gnm3.ann1.NPD7"
ALIGNMENTS = "Cicer/arietinum/genome_alignments/CDCFrontier.gnm3.wga.PXV3"
EXPRESSION = (
    "Cicer/arietinum/expression/CDCFrontier.gnm3.ann1.expr.ICC4958.Singh_Garg_2013"
)
QTL = "Cicer/arietinum/qtl/Demo.qtl.X_2000"

GENOME_BUSCO = {
    "results": {
        "Complete": 98.8,
        "Single copy": 97.0,
        "Multi copy": 1.8,
        "Fragmented": 0.2,
        "Missing": 1.0,
        "n_markers": 5366,
        "domain": "eukaryota",
        "Number of scaffolds": "2095",
        "Number of contigs": "2491",
        "Total length": "638661011",
        "Percent gaps": "0.006%",
        "Scaffold N50": "65305557",
    }
}

# A protein-mode summary carries completeness only, no assembly metrics.
ANNOTATION_BUSCO = {
    "results": {
        "Complete": 99.1,
        "Single copy": 96.5,
        "Multi copy": 2.6,
        "Fragmented": 0.3,
        "Missing": 0.6,
        "n_markers": 5366,
        "domain": "eukaryota",
    }
}


def _checksum(*names):
    return "\n".join(f"{MD5}  ./{name}" for name in names) + "\n"


@pytest.fixture(name="clone")
def fixture_clone(tmp_path, write):
    """One genus and species, with one collection of each type the subcommands turn
    into artifacts -- one per type, so the filesystem's listing order cannot reorder
    the output -- plus a qtl collection for the catalog's file prediction."""
    root = tmp_path / "datastore-metadata"

    def put(relative, content):
        text = content if isinstance(content, str) else yaml.safe_dump(content)
        write(os.path.join(root, relative), text)

    put(
        "Cicer/GENUS/about_this_collection/description_Cicer.yml",
        {
            "taxid": 3826,
            "genus": "Cicer",
            "commonName": "Chickpea",
            "description": "Demo genus description.",
            "species": ["arietinum"],
            "resources": [
                {"name": "Demo genus resource", "URL": "https://example.org/cicer"}
            ],
        },
    )
    put(
        "Cicer/arietinum/about_this_collection/description_Cicer_arietinum.yml",
        {
            "taxid": 3827,
            "genus": "Cicer",
            "species": "arietinum",
            "abbrev": "cicar",
            "commonName": "chickpea",
            "description": "Demo species description.",
            "strains": [
                {
                    "identifier": "CDCFrontier",
                    "name": "CDC Frontier",
                    "origin": "Canada",
                    "resources": [
                        {
                            "name": "Demo strain resource",
                            "URL": "https://example.org/cdc",
                        }
                    ],
                }
            ],
        },
    )
    put(
        f"{GENOME}/README.CDCFrontier.gnm3.QT0P.yml",
        {
            "identifier": "CDCFrontier.gnm3.QT0P",
            "synopsis": "Demo genome assembly",
            "scientific_name": "Cicer arietinum",
            "taxid": 3827,
            "scientific_name_abbrev": "cicar",
            "genotype": ["CDCFrontier"],
            "chromosome_prefix": "chr",
            "supercontig_prefix": "tig",
            "publication_doi": "10.1000/demo.genome",
            "license": "Open",
        },
    )
    put(
        f"{GENOME}/BUSCO/cicar.CDCFrontier.gnm3.QT0P.busco.fabales_odb10.short_summary.json",
        json.dumps(GENOME_BUSCO),
    )
    put(
        f"{ANNOTATION}/README.CDCFrontier.gnm3.ann1.NPD7.yml",
        {
            "identifier": "CDCFrontier.gnm3.ann1.NPD7",
            "synopsis": "Demo genome annotation",
            "scientific_name": "Cicer arietinum",
            "taxid": 3827,
            "scientific_name_abbrev": "cicar",
            "publication_doi": "10.1000/demo.annotation",
            "license": "Open",
        },
    )
    put(
        f"{ANNOTATION}/BUSCO/cicar.CDCFrontier.gnm3.ann1.NPD7.busco.fabales_odb10"
        ".short_summary.json",
        json.dumps(ANNOTATION_BUSCO),
    )
    put(
        f"{ALIGNMENTS}/README.CDCFrontier.gnm3.wga.PXV3.yml",
        {"identifier": "CDCFrontier.gnm3.wga.PXV3", "synopsis": "Demo alignments"},
    )
    pair = "cicar.CDCFrontier.gnm3.x.cicec.S2Drd065.gnm1.PXV3.minimap2"
    put(
        f"{ALIGNMENTS}/CHECKSUM.CDCFrontier.gnm3.wga.PXV3.md5",
        _checksum(f"{pair}.paf.gz", f"{pair}.bam", f"{pair}.bam.bai"),
    )
    put(
        f"{EXPRESSION}/README.CDCFrontier.gnm3.ann1.expr.ICC4958.Singh_Garg_2013.yml",
        {
            "identifier": "CDCFrontier.gnm3.ann1.expr.ICC4958.Singh_Garg_2013",
            "synopsis": "Demo expression study",
            "scientific_name": "Cicer arietinum",
            "taxid": 3827,
            "scientific_name_abbrev": "cicar",
            "publication_doi": "10.1000/demo.expression",
            "bioproject": "PRJNA000001",
            "expression_unit": "TPM",
            "genotype": ["ICC4958"],
        },
    )
    samples = "cicar.CDCFrontier.gnm3.ann1.expr.ICC4958.Singh_Garg_2013"
    put(
        f"{EXPRESSION}/CHECKSUM.CDCFrontier.gnm3.ann1.expr.ICC4958.Singh_Garg_2013.md5",
        _checksum(
            f"{samples}.Flower.bw",
            f"{samples}.values.tsv.gz",
            f"{samples}.Young_leaves.bw",
        ),
    )
    put(
        f"{QTL}/README.Demo.qtl.X_2000.yml",
        {
            "identifier": "Demo.qtl.X_2000",
            "synopsis": "Demo QTL",
            "scientific_name_abbrev": "cicar",
        },
    )
    return str(root)


@pytest.fixture(name="taxa_list")
def fixture_taxa_list(tmp_path):
    path = tmp_path / "taxon_list.yml"
    path.write_text(yaml.safe_dump([{"genus": "Cicer", "description": "Chickpea"}]))
    return str(path)


@pytest.fixture(name="network")
def fixture_network(monkeypatch):
    """Answer the only remote requests these subcommands make: a genome's .fai, to place
    JBrowse2 sessions, and HEAD checks for protein files, of which only the full set
    exists here. Any other request fails the test, so no run can quietly depend on the
    live datastore."""

    def fake_get(url, timeout=None):  # pylint: disable=unused-argument
        if url.endswith(".genome_main.fna.gz.fai"):
            text = "chr1\t48359943\t6\t60\t61\nchr2\t36633317\t49166994\t60\t61\n"
            return SimpleNamespace(status_code=200, text=text)
        raise AssertionError(f"unexpected GET {url}")

    def fake_head(url, timeout=None):  # pylint: disable=unused-argument
        if url.endswith(".protein.faa.gz"):
            return SimpleNamespace(status_code=200, text="")
        if url.endswith(".protein_primary.faa.gz"):
            return SimpleNamespace(status_code=404, text="")
        raise AssertionError(f"unexpected HEAD {url}")

    monkeypatch.setattr(requests, "get", fake_get)
    monkeypatch.setattr(requests, "head", fake_head)


def _invoke(command, args):
    result = CliRunner().invoke(command, [str(arg) for arg in args])
    detail = (
        "".join(traceback.format_exception(*result.exc_info)) if result.exc_info else ""
    )
    assert result.exit_code == 0, f"{result.output}\n{detail}"
    return result


def _check(name, actual):
    """Compare an artifact with its reviewed copy, or rewrite the copy when regenerating."""
    path = os.path.join(EXPECTED, name)
    if UPDATING:
        os.makedirs(EXPECTED, exist_ok=True)
        with open(path, "w", encoding="utf-8") as handle:
            handle.write(actual)
        return
    with open(path, encoding="utf-8") as handle:
        assert actual == handle.read(), f"{name} differs from its reviewed copy"


def _skip_if_updating():
    if UPDATING:
        pytest.skip(f"regenerated artifacts in {EXPECTED}; review the diff")


@pytest.mark.parametrize("verify", [False, True])
def test_populate_catalog_passes_its_options_to_the_builder(
    clone, tmp_path, monkeypatch, verify
):
    """`--verify` was declared but missing from the function signature, so every run
    failed with TypeError. The document itself is pinned in test_catalog.py; this checks
    that each option reaches the builder (the qtl collection is `verified` only if
    --verify does)."""
    monkeypatch.setattr(DatastoreIndex, "url_exists", lambda self, url: True)
    out = tmp_path / "out" / "catalog.json"
    args = [
        "--from_github",
        clone,
        "--catalog_out",
        out,
        "--datastore_url",
        "https://example.org/ds",
        "--log_file",
        tmp_path / "catalog.log",
    ]
    _invoke(lis_cli.populate_catalog, args + (["--verify"] if verify else []))
    written = json.loads(out.read_text())
    builder = CatalogBuilder(
        clone, datastore_url="https://example.org/ds", verify=verify
    )
    expected = json.loads(json.dumps(builder.build()))
    for document in (written, expected):
        del document["built_at"]
    assert written == expected
    status = {c["path"]: c["index_status"] for c in written["collections"]}[QTL]
    assert status == ("verified" if verify else "inferred")


@pytest.mark.usefixtures("network")
@pytest.mark.parametrize("with_taxa_list", [True, False])
def test_populate_dscensor_writes_the_reviewed_nodes(
    clone, taxa_list, tmp_path, with_taxa_list
):
    """One node per genome, annotation, protein set that exists, whole-genome alignment
    and bigwig, each carrying the README fields consumers read. The nodes are the same
    whether genera come from a taxa list or are discovered from the clone; passing the
    two in the wrong order once made every run fail."""
    out = tmp_path / "nodes"
    args = [
        "--nodes_out",
        out,
        "--from_github",
        clone,
        "--log_file",
        tmp_path / "d.log",
    ]
    if with_taxa_list:
        args += ["--taxa_list", taxa_list]
    _invoke(lis_cli.populate_dscensor, args)
    nodes = {
        name: json.loads((out / name).read_text())
        for name in sorted(os.listdir(out))
        if name.endswith(".json")
    }
    _check("dscensor_nodes.json", json.dumps(nodes, indent=2, sort_keys=True) + "\n")
    _skip_if_updating()


@pytest.mark.usefixtures("network")
def test_populate_blast_prints_the_reviewed_commands(clone, taxa_list, tmp_path):
    out = tmp_path / "blast"
    result = _invoke(
        lis_cli.populate_blast,
        [
            "--taxa_list",
            taxa_list,
            "--blast_out",
            out,
            "--from_github",
            clone,
            "--cmds_only",
            "--log_file",
            tmp_path / "blast.log",
        ],
    )
    commands = [
        line.replace(str(out), "<OUT>")
        for line in result.stdout.splitlines()
        if line.startswith("set -o pipefail")
    ]
    _check("blast_commands.txt", "\n".join(commands) + "\n")
    _skip_if_updating()


@pytest.mark.usefixtures("network")
def test_populate_jbrowse2_prints_the_reviewed_commands(clone, taxa_list, tmp_path):
    out = tmp_path / "jbrowse"
    result = _invoke(
        lis_cli.populate_jbrowse2,
        [
            "--jbrowse_url",
            "https://example.org/jbrowse2",
            "--taxa_list",
            taxa_list,
            "--jbrowse_out",
            out,
            "--from_github",
            clone,
            "--cmds_only",
            "--log_file",
            tmp_path / "jbrowse.log",
        ],
    )
    commands = [
        line.replace(str(out), "<OUT>")
        for line in result.stdout.splitlines()
        if line.startswith("jbrowse ")
    ]
    _check("jbrowse2_commands.txt", "\n".join(commands) + "\n")
    _skip_if_updating()


@pytest.mark.usefixtures("network")
def test_populate_jekyll_writes_the_reviewed_yaml(clone, taxa_list, tmp_path):
    out = tmp_path / "jekyll"
    _invoke(
        lis_cli.populate_jekyll,
        [
            "--taxa_list",
            taxa_list,
            "--collections_out",
            out,
            "--from_github",
            clone,
            "--log_file",
            tmp_path / "jekyll.log",
        ],
    )
    for name in (
        "genus_resources.yml",
        "species_collections.yml",
        "species_resources.yml",
    ):
        _check(f"jekyll_{name}", (out / "Cicer" / name).read_text())
    _skip_if_updating()
