"""Tests for the catalog document built from the datastore index.

These pin what consumers read: the document shape, the index-status and `src` labels,
index flags, README propagation, lineage and inheritance. Parsing details of the index
itself are covered in test_datastore_files.py; the fixture tree is in conftest.py.
"""

import json
import os

import pytest
import yaml

from lis_autocontent.catalog import SCHEMA_VERSION, CatalogBuilder
from lis_autocontent.datastore_files import README_FIELDS

GENOME = "Glycine/max/genomes/Wm82.gnm4.4PTR"
ANNOTATION = "Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"
QTL = "Glycine/max/qtl/Demo.qtl.X_1990"
MD5 = "d41d8cd98f00b204e9800998ecf8427e"

# One value for every README field the catalog publishes, in the types the datastore
# uses (taxid is an integer, genotype and keywords are lists).
EVERY_README_FIELD = {
    "synopsis": "Demo synopsis",
    "description": "Demo description",
    "scientific_name": "Glycine max",
    "taxid": 3847,
    "scientific_name_abbrev": "glyma",
    "genotype": ["Williams 82"],
    "publication_doi": "10.1111/tpj.14500",
    "publication_title": "Demo publication",
    "citation": "Demo et al., 2020",
    "license": "Open",
    "public_access_level": "public",
    "keywords": ["soybean", "genome"],
    "chromosome_prefix": "Gm",
    "supercontig_prefix": "scaffold",
    "bioproject": "PRJNA19861",
    "genbank_accession": "GCA_000004515.4",
    "sraproject": "SRP000001",
    "genetic_map": "GmComposite1999",
    "expression_unit": "TPM",
    "dataset_doi": "10.5281/zenodo.1",
    "related_to": ["Wm82.gnm4.4PTR"],
    "source": "https://example.org",
    "dataset_release_date": "2020-01-01",
}


@pytest.fixture(name="catalog")
def fixture_catalog(metadata_dir):
    return CatalogBuilder(metadata_dir).build()


def _by_path(catalog):
    return {c["path"]: c for c in catalog["collections"]}


# --- document shape ----------------------------------------------------------
EXPECTED_CATALOG = os.path.join(
    os.path.dirname(__file__), "data", "catalog.expected.json"
)

TRAITS = """---
scientific_name: Glycine max
gene_symbols:
  - GmNARK
  - NTS-1
gene_model_full_id: glyma.Wm82.gnm4.ann1.Glyma.12G040000
phenotype_synopsis: Demo nodulation phenotype
references:
  - citation: First demo reference
    doi: 10.1000/demo.first
  - citation: Second demo reference
    doi: 10.1000/demo.second
---
scientific_name: Glycine max
gene_symbols:
  - GmFT2a
gene_model_full_id: glyma.Wm82.gnm4.ann1.Glyma.16G150700
references:
  - citation: A reference without a DOI
"""


def _complete_the_tree(metadata_dir, write):
    """Add the parts of the document no other fixture populates: a data file in a
    subdirectory (1,127 real ones, all under BUSCO/), a synteny edge, curated gene
    symbols and taxon descriptions."""
    checksum = os.path.join(
        metadata_dir, ANNOTATION, "CHECKSUM.Wm82.gnm4.ann1.T8TQ.md5"
    )
    with open(checksum, "a", encoding="utf-8") as handle:
        handle.write(
            f"\n{MD5}  ./BUSCO/glyma.Wm82.gnm4.ann1.T8TQ.busco.fabales_odb10"
            ".short_summary.json\n"
        )
    synteny = os.path.join(metadata_dir, "Glycine/max/synteny/Wm82.gnm4.synt.PXV3")
    partner = "glyma.Wm82.gnm4.x.phavu.G19833.gnm2.PXV3.gff3.gz"
    write(
        os.path.join(synteny, "CHECKSUM.Wm82.gnm4.synt.PXV3.md5"),
        f"{MD5}  ./{partner}\n{MD5}  ./{partner}.tbi\n",
    )
    write(
        os.path.join(metadata_dir, "Glycine/max/gene_functions/glyma.traits.yml"),
        TRAITS,
    )
    write(
        os.path.join(
            metadata_dir,
            "Glycine/max/about_this_collection/description_Glycine_max.yml",
        ),
        "---\ntaxid: 3847\ngenus: Glycine\nspecies: max\nabbrev: glyma\n"
        "commonName: soybean\n",
    )
    write(
        os.path.join(
            metadata_dir, "Glycine/GENUS/about_this_collection/description_Glycine.yml"
        ),
        "---\ntaxid: 3846\ngenus: Glycine\nspecies:\n  - max\n",
    )


@pytest.mark.usefixtures("markers_collection")
def test_the_document_matches_the_reviewed_expected_output(metadata_dir, write):
    """The catalog is a published contract, so the whole document is pinned against a
    reviewed copy in tests/data. The named tests explain individual rules; this one
    catches the changes they do not anticipate -- a renamed key, a lost subdirectory in
    a file path, a dropped gene symbol, a miscounted stat, a reordered list.

    When a change to the document is intended, regenerate the copy with
    `LIS_UPDATE_EXPECTED=1 pytest tests/test_catalog.py` and review its diff before
    committing it: regenerating without reading the diff defeats the test."""
    _complete_the_tree(metadata_dir, write)
    document = json.loads(json.dumps(CatalogBuilder(metadata_dir).build()))
    for volatile in ("built_at", "source_commit"):  # clock and checkout dependent
        del document[volatile]
    if os.environ.get("LIS_UPDATE_EXPECTED"):
        os.makedirs(os.path.dirname(EXPECTED_CATALOG), exist_ok=True)
        with open(EXPECTED_CATALOG, "w", encoding="utf-8") as handle:
            json.dump(document, handle, indent=2)
            handle.write("\n")
        pytest.skip(f"regenerated {EXPECTED_CATALOG}; review the diff")
    with open(EXPECTED_CATALOG, encoding="utf-8") as handle:
        assert document == json.load(handle)


def test_document_carries_schema_and_build_stamp(catalog):
    """A consumer must be able to tell which catalog it is reasoning over."""
    assert catalog["schema"] == SCHEMA_VERSION
    assert catalog["built_at"].endswith("Z")
    assert "source_commit" in catalog  # None outside a git checkout, but present


# --- README propagation ------------------------------------------------------
@pytest.mark.parametrize(
    "ctype, key, has_checksum",
    [("genomes", "Wm82.gnm9.ABCD", True), ("qtl", "Demo.qtl.Y_2000", False)],
)
def test_every_readme_field_is_carried_verbatim(
    metadata_dir, write, ctype, key, has_checksum
):
    """The regression this exists to prevent: taxid and publication_doi were both
    parsed and then dropped before reaching a consumer. Whether a collection has a
    CHECKSUM changes its file list, never its metadata."""
    assert set(EVERY_README_FIELD) == set(README_FIELDS)  # keep this table complete
    path = f"Glycine/max/{ctype}/{key}"
    coll = os.path.join(metadata_dir, path)
    readme = dict(EVERY_README_FIELD, identifier=key)
    write(os.path.join(coll, f"README.{key}.yml"), yaml.safe_dump(readme))
    if has_checksum:
        write(
            os.path.join(coll, f"CHECKSUM.{key}.md5"),
            f"{MD5}  ./glyma.{key}.genome_main.fna.gz\n",
        )
    record = _by_path(CatalogBuilder(metadata_dir).build())[path]
    assert {field: record.get(field) for field in README_FIELDS} == EVERY_README_FIELD


# --- index status: how the file list was obtained -----------------------------
def test_index_status_and_src_record_how_each_file_list_was_obtained(catalog):
    """Consumers rely on the difference: `known`/`checksum` came from a CHECKSUM,
    `inferred`/`predicted` were constructed from the documented convention, `unknown`
    means neither was possible. Collapsing them presents a guess as a fact. (`verified`
    is covered by test_verify_keeps_only_files_that_exist.)"""
    records = _by_path(catalog)
    observed = {
        path: (
            records[path]["index_status"],
            sorted({f["src"] for f in records[path]["files"]}),
        )
        for path in (ANNOTATION, QTL, GENOME)
    }
    assert observed == {
        ANNOTATION: ("known", ["checksum"]),
        QTL: ("inferred", ["predicted"]),
        GENOME: ("unknown", []),
    }


def test_index_siblings_become_flags_not_entries(catalog):
    """.fai/.tbi/.gzi are evidence about other files, never data files themselves, and
    README/MANIFEST/CHECKSUM are metadata: only the 3 data files are listed."""
    annotation = _by_path(catalog)[ANNOTATION]
    names = [f["n"] for f in annotation["files"]]
    assert len(names) == 3  # protein, gene models, gene families
    assert not any(n.endswith((".fai", ".tbi", ".gzi")) for n in names)
    indexed = {f["n"]: f["i"] for f in annotation["files"] if "i" in f}
    assert indexed == {
        "glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz": [".fai"],
        "glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz": [".tbi"],
    }


# --- MANIFEST ----------------------------------------------------------------
def test_manifest_descriptions_and_applications_merge_onto_files(catalog):
    annotation = _by_path(catalog)[ANNOTATION]
    protein = [f for f in annotation["files"] if "protein_primary" in f["n"]][0]
    assert protein["description"] == "Protein sequences - primary only"
    assert protein["applications"] == ["blast", "mines"]


def test_manifest_placeholder_description_is_dropped(catalog):
    """ "MISSING" is the spec's placeholder, not a description."""
    annotation = _by_path(catalog)[ANNOTATION]
    families = [f for f in annotation["files"] if "gfa.tsv" in f["n"]][0]
    assert "description" not in families


# --- BUSCO and counts --------------------------------------------------------
def test_metrics_are_not_emitted_yet(catalog):
    """Deliberately deferred: see the NOTE in catalog.py. This test exists so the
    omission is a decision on record rather than something that quietly regresses
    in either direction."""
    genome = _by_path(catalog)[GENOME]
    assert "busco" not in genome
    assert "counts" not in genome


# --- lineage and inheritance -------------------------------------------------
def test_annotation_links_to_its_genome(catalog):
    annotation = _by_path(catalog)[ANNOTATION]
    assert annotation["derived_from"] == ["Wm82.gnm4.4PTR"]


def test_assembly_conventions_are_inherited_down_the_edge(catalog):
    """These live on the genome README only. Without inheritance a caller cannot turn
    `Gm12` into `glyma.Wm82.gnm4.Gm12` from an annotation."""
    annotation = _by_path(catalog)[ANNOTATION]
    assert annotation["inherited"] == [
        "chromosome_prefix",
        "supercontig_prefix",
        "bioproject",
    ]
    assert (
        annotation["chromosome_prefix"],
        annotation["supercontig_prefix"],
        annotation["bioproject"],
    ) == ("Gm", "scaffold", "PRJNA19861")


def test_inheritance_never_overwrites_a_published_value(metadata_dir):
    """An annotation that publishes its own prefix keeps it."""
    readme = os.path.join(metadata_dir, ANNOTATION, "README.Wm82.gnm4.ann1.T8TQ.yml")
    with open(readme, "a", encoding="utf-8") as handle:
        handle.write("chromosome_prefix: Chr\n")
    record = _by_path(CatalogBuilder(metadata_dir).build())[ANNOTATION]
    assert record["chromosome_prefix"] == "Chr"
    assert "chromosome_prefix" not in record.get("inherited", [])


def test_genomes_do_not_derive_from_themselves(catalog):
    genome = _by_path(catalog)[GENOME]
    assert "derived_from" not in genome


def test_pairwise_relationships_are_derived_from_file_names(metadata_dir, write):
    """A pairwise file is stored once, under its reference genome; the catalog lifts
    synteny and whole-genome alignments into one edge list sorted by genome pair,
    keeping any duplication epoch. The alignment collection sorts first on disk, so the
    expected order below only holds if the list is actually sorted."""
    synteny = "Glycine/max/synteny/Wm82.gnm4.synt.PXV3"
    partner = "glyma.Wm82.gnm4.x.phavu.G19833.gnm2.PXV3.gff3.gz"
    self_pair = "glyma.Wm82.gnm4.x.glyma.Wm82.gnm4.old_duplication.PXV3.gff3.gz"
    write(
        os.path.join(metadata_dir, synteny, "CHECKSUM.Wm82.gnm4.synt.PXV3.md5"),
        "\n".join(
            f"{MD5}  ./{name}" for name in (partner, partner + ".tbi", self_pair)
        ),
    )
    alignments = "Glycine/max/genome_alignments/Wm82.gnm4.wga.LXVF"
    stem = "glyma.Wm82.gnm4.x.vigun.IT97K-499-35.gnm1.LXVF"
    write(
        os.path.join(metadata_dir, alignments, "CHECKSUM.Wm82.gnm4.wga.LXVF.md5"),
        "\n".join(f"{MD5}  ./{stem}.{ext}" for ext in ("paf.gz", "bam", "bam.bai")),
    )
    pairs = CatalogBuilder(metadata_dir).build()["pairwise"]
    assert [
        (
            p["a"],
            p["b"],
            p["kind"],
            p["format"],
            p.get("epoch"),
            p.get("self"),
            p.get("i"),
        )
        for p in pairs
    ] == [
        (
            "glyma.Wm82.gnm4",
            "glyma.Wm82.gnm4",
            "synteny",
            "gff3.gz",
            "old_duplication",
            True,
            None,
        ),
        (
            "glyma.Wm82.gnm4",
            "phavu.G19833.gnm2",
            "synteny",
            "gff3.gz",
            None,
            None,
            [".tbi"],
        ),
        (
            "glyma.Wm82.gnm4",
            "vigun.IT97K-499-35.gnm1",
            "alignment",
            "bam",
            None,
            None,
            [".bai"],
        ),
        (
            "glyma.Wm82.gnm4",
            "vigun.IT97K-499-35.gnm1",
            "alignment",
            "paf.gz",
            None,
            None,
            None,
        ),
    ]
    assert pairs[1]["url"] == f"https://data.legumeinfo.org/{synteny}/{partner}"
    assert pairs[3]["collection"] == alignments


# --- urls and output ---------------------------------------------------------
@pytest.mark.parametrize(
    "datastore_url, expected",
    [
        (None, "https://data.legumeinfo.org"),
        ("https://example.org/ds/", "https://example.org/ds"),
    ],
)
def test_base_url_joins_the_datastore_url_and_collection_path(
    metadata_dir, datastore_url, expected
):
    catalog = CatalogBuilder(metadata_dir, datastore_url=datastore_url).build()
    assert catalog["datastore_url"] == expected
    assert _by_path(catalog)[ANNOTATION]["base_url"] == f"{expected}/{ANNOTATION}"


def test_write_round_trips(metadata_dir, tmp_path):
    out = str(tmp_path / "out" / "catalog.json")
    CatalogBuilder(metadata_dir).write(out)
    with open(out, encoding="utf-8") as handle:
        assert json.load(handle)["stats"]["collections"] == 3


def test_unreadable_readme_degrades_the_record_not_the_build(metadata_dir, write):
    """One malformed file must not sink a 1,000-collection build — and must not drop
    the collection either. A real datastore collection ships a README that is not
    valid YAML; its files are still real and still worth cataloguing. With no README
    there is no abbrev, so no filename is predicted either."""
    broken = os.path.join(metadata_dir, "Glycine", "max", "maps", "Bad.map.X")
    write(os.path.join(broken, "README.Bad.map.X.yml"), "---\n: : not yaml : :\n")
    catalog = CatalogBuilder(metadata_dir).build()
    assert catalog["stats"]["collections"] == 4
    record = _by_path(catalog)["Glycine/max/maps/Bad.map.X"]
    assert record["index_status"] == "unknown"  # no CHECKSUM, nothing predictable
    assert record["files"] == []
    assert "publication_doi" not in record  # nothing was invented


def test_a_collection_with_only_a_checksum_is_catalogued(metadata_dir, write):
    """Five genome_alignments collections publish a CHECKSUM and no README. Keying off
    the README alone dropped them, along with the indexed BAMs they hold."""
    wga = os.path.join(
        metadata_dir, "Arachis", "hypogaea", "genome_alignments", "X.wga.AAAA"
    )
    write(
        os.path.join(wga, "CHECKSUM.X.wga.AAAA.md5"),
        f"{MD5}  ./arahy.x.y.AAAA.bam\n{MD5}  ./arahy.x.y.AAAA.bam.bai",
    )
    catalog = CatalogBuilder(metadata_dir).build()
    record = _by_path(catalog)["Arachis/hypogaea/genome_alignments/X.wga.AAAA"]
    assert record["index_status"] == "known"
    assert record["files"] == [
        {"n": "arahy.x.y.AAAA.bam", "i": [".bai"], "src": "checksum"}
    ]


def test_a_bare_README_is_still_parsed(metadata_dir, write):
    """At least one collection ships `README` with no extension; it is still YAML."""
    wga = os.path.join(
        metadata_dir, "Arachis", "ipaensis", "genome_alignments", "Y.wga.BBBB"
    )
    write(
        os.path.join(wga, "README"),
        "identifier: Y.wga.BBBB\nsynopsis: Genome alignments\npublication_doi: 10.1/x\n",
    )
    catalog = CatalogBuilder(metadata_dir).build()
    record = _by_path(catalog)["Arachis/ipaensis/genome_alignments/Y.wga.BBBB"]
    assert record["synopsis"] == "Genome alignments"
    assert record["publication_doi"] == "10.1/x"


def test_taxa_descriptions_are_ingested(metadata_dir, write):
    """Common name, taxid and abbreviation live only in about_this_collection, so
    without this a consumer cannot map 'soybean' onto a taxon at all."""
    about = os.path.join(metadata_dir, "Glycine", "max", "about_this_collection")
    write(
        os.path.join(about, "description_Glycine_max.yml"),
        "---\ntaxid: 3847\ngenus: Glycine\nspecies: max\nabbrev: glyma\n"
        "commonName: soybean\ndescription: Soybean.\n"
        "resources:\n  - name: GlycineMine\n    URL: https://mines.legumeinfo.org/glycinemine\n"
        "    description: InterMine for Glycine\n",
    )
    catalog = CatalogBuilder(metadata_dir).build()
    assert catalog["taxa"]["Glycine/max"] == {
        "taxid": 3847,
        "abbrev": "glyma",
        "commonName": "soybean",
        "description": "Soybean.",
        "genus": "Glycine",
        "species": "max",
        "resources": [
            {
                "name": "GlycineMine",
                "URL": "https://mines.legumeinfo.org/glycinemine",
                "description": "InterMine for Glycine",
            }
        ],
    }


def test_genus_level_description_lists_its_species(metadata_dir, write):
    """A genus file uses `species` as a LIST; keying it like a species file would
    conflate the two."""
    about = os.path.join(metadata_dir, "Glycine", "GENUS", "about_this_collection")
    write(
        os.path.join(about, "description_Glycine.yml"),
        "---\ntaxid: 3847\ngenus: Glycine\ncommonName: soybean\n"
        "species:\n  - max\n  - soja\n",
    )
    catalog = CatalogBuilder(metadata_dir).build()
    assert catalog["taxa"]["Glycine"]["species_in_genus"] == ["max", "soja"]
    assert "species" not in catalog["taxa"]["Glycine"]


# --- file resolution for collections that publish no CHECKSUM -------------------------
def test_verify_keeps_only_files_that_exist(metadata_dir, monkeypatch):
    """--verify turns prediction into evidence: absent files are dropped and the
    survivors are relabelled, so a consumer can tell the two apart."""
    real = {"glyma.Demo.qtl.X_1990.qtl.tsv.gz", "glyma.Demo.qtl.X_1990.obo.tsv.gz"}
    builder = CatalogBuilder(metadata_dir, verify=True)
    monkeypatch.setattr(
        builder.index, "url_exists", lambda url: url.rsplit("/", 1)[-1] in real
    )
    catalog = builder.build()
    qtl = _by_path(catalog)[QTL]
    assert qtl["index_status"] == "verified"
    assert {f["n"] for f in qtl["files"]} == real
    assert {f["src"] for f in qtl["files"]} == {"verified"}
    # All four candidates were asked about, and no index siblings: qtl is unindexed.
    assert catalog["stats"]["probes"] == 4
    assert catalog["stats"]["verified_files"] == 2


def test_offline_index_status_of_an_indexable_type_is_unknown(
    metadata_dir, markers_collection
):
    """A predicted file of a type that can be indexed was never probed, so whether it is
    streamable is unknown -- not "no index"."""
    record = _by_path(CatalogBuilder(metadata_dir).build())[markers_collection]
    assert record["files"] == [
        {"n": "glyma.SoySNP50K.mrk.ABCD.gff3.gz", "src": "predicted", "i_unknown": True}
    ]


def test_verify_attaches_the_index_siblings_that_exist(
    metadata_dir, markers_collection, monkeypatch
):
    """Under --verify the unknown is resolved: each surviving file is probed for the
    siblings its extension allows, and only the ones found are recorded."""
    real = {"glyma.SoySNP50K.mrk.ABCD.gff3.gz", "glyma.SoySNP50K.mrk.ABCD.gff3.gz.tbi"}
    builder = CatalogBuilder(metadata_dir, verify=True)
    monkeypatch.setattr(
        builder.index, "url_exists", lambda url: url.rsplit("/", 1)[-1] in real
    )
    catalog = builder.build()
    assert _by_path(catalog)[markers_collection]["files"] == [
        {"n": "glyma.SoySNP50K.mrk.ABCD.gff3.gz", "src": "verified", "i": [".tbi"]}
    ]
    # qtl: 4 file probes and no index probes (indexed: false); markers: 1 + .tbi/.csi.
    assert catalog["stats"]["probes"] == 7
