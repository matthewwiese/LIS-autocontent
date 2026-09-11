"""Tests for the catalog document built from the datastore index.

These pin what consumers read: the document shape, the index-status and `src` labels,
index flags, README propagation, lineage and inheritance. Parsing details of the index
itself are covered in test_datastore_files.py; the fixture tree is in conftest.py.
"""

import json
import os

import pytest

from lis_autocontent.catalog import SCHEMA_VERSION, CatalogBuilder

GENOME = "Glycine/max/genomes/Wm82.gnm4.4PTR"
ANNOTATION = "Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"
QTL = "Glycine/max/qtl/Demo.qtl.X_1990"
MD5 = "d41d8cd98f00b204e9800998ecf8427e"


@pytest.fixture(name="catalog")
def fixture_catalog(metadata_dir):
    return CatalogBuilder(metadata_dir).build()


def _by_path(catalog):
    return {c["path"]: c for c in catalog["collections"]}


# --- document shape ----------------------------------------------------------
def test_document_carries_schema_and_build_stamp(catalog):
    """A consumer must be able to tell which catalog it is reasoning over."""
    assert catalog["schema"] == SCHEMA_VERSION
    assert catalog["built_at"].endswith("Z")
    assert "source_commit" in catalog  # None outside a git checkout, but present
    assert catalog["stats"]["collections"] == 3


def test_every_collection_is_found(catalog):
    assert sorted(_by_path(catalog)) == [ANNOTATION, GENOME, QTL]


# --- README propagation ------------------------------------------------------
def test_readme_fields_are_carried_verbatim(catalog):
    """The regression this exists to prevent: taxid and publication_doi were both
    parsed and then dropped before reaching a consumer."""
    genome = _by_path(catalog)[GENOME]
    assert genome["taxid"] == 3847
    assert genome["publication_doi"] == "10.1111/tpj.14500"
    assert genome["scientific_name_abbrev"] == "glyma"
    assert genome["license"] == "Open"
    assert genome["genotype"] == ["Williams 82"]


def test_qtl_keeps_its_metadata_without_a_checksum(catalog):
    """No CHECKSUM costs the file list, not the metadata."""
    qtl = _by_path(catalog)[QTL]
    assert qtl["publication_doi"] == "10.1007/bf00226154"
    assert qtl["genetic_map"] == "GmComposite1999"


# --- index status: how the file list was obtained -----------------------------
def test_index_status_distinguishes_how_the_file_list_was_obtained(catalog):
    """Four states, and consumers rely on the difference: `known` came from a CHECKSUM,
    `inferred`/`verified` were constructed from the documented convention, `unknown`
    means neither was possible. Collapsing them presents a guess as a fact."""
    qtl = _by_path(catalog)[QTL]
    assert qtl["index_status"] == "inferred"  # no CHECKSUM, vocabulary applied
    assert {f["src"] for f in qtl["files"]} == {"predicted"}

    annotation = _by_path(catalog)[ANNOTATION]
    assert annotation["index_status"] == "known"  # CHECKSUM published


def test_index_siblings_become_flags_not_entries(catalog):
    """.fai/.tbi/.gzi are evidence about other files, never data files themselves."""
    annotation = _by_path(catalog)[ANNOTATION]
    names = [f["n"] for f in annotation["files"]]
    assert len(names) == 3  # protein, gene models, gene families
    assert not any(n.endswith((".fai", ".tbi", ".gzi")) for n in names)
    indexed = {f["n"]: f["i"] for f in annotation["files"] if "i" in f}
    assert indexed == {
        "glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz": [".fai"],
        "glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz": [".tbi"],
    }


def test_metadata_files_are_not_listed_as_data(catalog):
    annotation = _by_path(catalog)[ANNOTATION]
    names = [f["n"] for f in annotation["files"]]
    assert not any(n.startswith(("README.", "MANIFEST.", "CHECKSUM.")) for n in names)


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


def test_busco_availability_is_still_reported_in_stats(catalog):
    """The build diagnostic must stay truthful while the metrics are withheld."""
    assert catalog["stats"]["with_busco"] == 1


# --- lineage and inheritance -------------------------------------------------
def test_annotation_links_to_its_genome(catalog):
    annotation = _by_path(catalog)[ANNOTATION]
    assert annotation["derived_from"] == ["Wm82.gnm4.4PTR"]


def test_assembly_conventions_are_inherited_down_the_edge(catalog):
    """chromosome_prefix lives on the genome README only. Without inheritance a
    caller cannot turn `Gm12` into `glyma.Wm82.gnm4.Gm12` from an annotation."""
    annotation = _by_path(catalog)[ANNOTATION]
    assert annotation["chromosome_prefix"] == "Gm"
    assert annotation["bioproject"] == "PRJNA19861"
    assert "chromosome_prefix" in annotation["inherited"]


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
    """A pairwise file is stored once, under its reference genome; the catalog lifts it
    into a top-level edge list, keeping any duplication epoch."""
    synteny = "Glycine/max/synteny/Wm82.gnm4.synt.PXV3"
    partner = "glyma.Wm82.gnm4.x.phavu.G19833.gnm2.PXV3.gff3.gz"
    self_pair = "glyma.Wm82.gnm4.x.glyma.Wm82.gnm4.old_duplication.PXV3.gff3.gz"
    write(
        os.path.join(metadata_dir, synteny, "CHECKSUM.Wm82.gnm4.synt.PXV3.md5"),
        "\n".join(
            f"{MD5}  ./{name}" for name in (partner, partner + ".tbi", self_pair)
        ),
    )
    pairs = CatalogBuilder(metadata_dir).build()["pairwise"]
    assert [(p["a"], p["b"], p.get("epoch"), p.get("self")) for p in pairs] == [
        ("glyma.Wm82.gnm4", "glyma.Wm82.gnm4", "old_duplication", True),
        ("glyma.Wm82.gnm4", "phavu.G19833.gnm2", None, None),
    ]
    assert pairs[1]["kind"] == "synteny"
    assert pairs[1]["i"] == [".tbi"]
    assert pairs[1]["url"] == f"https://data.legumeinfo.org/{synteny}/{partner}"


# --- urls and output ---------------------------------------------------------
def test_base_url_points_at_the_datastore(catalog):
    annotation = _by_path(catalog)[ANNOTATION]
    assert annotation["base_url"] == f"https://data.legumeinfo.org/{ANNOTATION}"


def test_datastore_url_is_overridable(metadata_dir):
    catalog = CatalogBuilder(
        metadata_dir, datastore_url="https://example.org/ds/"
    ).build()
    record = _by_path(catalog)[GENOME]
    assert record["base_url"].startswith("https://example.org/ds/Glycine")


def test_write_round_trips(metadata_dir, tmp_path):
    out = str(tmp_path / "out" / "catalog.json")
    CatalogBuilder(metadata_dir).write(out)
    with open(out, encoding="utf-8") as handle:
        assert json.load(handle)["stats"]["collections"] == 3


def test_unreadable_readme_degrades_the_record_not_the_build(metadata_dir, write):
    """One malformed file must not sink a 1,000-collection build — and must not drop
    the collection either. A real datastore collection ships a README that is not
    valid YAML; its files are still real and still worth cataloguing."""
    broken = os.path.join(metadata_dir, "Glycine", "max", "maps", "Bad.map.X")
    write(os.path.join(broken, "README.Bad.map.X.yml"), "---\n: : not yaml : :\n")
    catalog = CatalogBuilder(metadata_dir).build()
    assert catalog["stats"]["collections"] == 4
    record = _by_path(catalog)["Glycine/max/maps/Bad.map.X"]
    assert record["index_status"] == "unknown"  # no CHECKSUM either
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
    """Common name and abbreviation live only in about_this_collection, so without
    this a consumer cannot map 'soybean' onto a taxon at all."""
    about = os.path.join(metadata_dir, "Glycine", "max", "about_this_collection")
    write(
        os.path.join(about, "description_Glycine_max.yml"),
        "---\ntaxid: 3847\ngenus: Glycine\nspecies: max\nabbrev: glyma\n"
        "commonName: soybean\nresources:\n  - name: GlycineMine\n    URL: https://x\n",
    )
    catalog = CatalogBuilder(metadata_dir).build()
    entry = catalog["taxa"]["Glycine/max"]
    assert entry["commonName"] == "soybean"
    assert entry["abbrev"] == "glyma"
    assert entry["resources"][0]["name"] == "GlycineMine"


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
def test_offline_build_labels_files_as_predicted(catalog):
    """Without --verify nothing is probed, so the label must not claim more."""
    qtl = _by_path(catalog)[QTL]
    assert qtl["index_status"] == "inferred"
    assert qtl["files"]
    assert {f["src"] for f in qtl["files"]} == {"predicted"}


def test_checksum_files_are_labelled_authoritative(catalog):
    annotation = _by_path(catalog)[ANNOTATION]
    assert annotation["index_status"] == "known"
    assert {f["src"] for f in annotation["files"]} == {"checksum"}


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
    assert catalog["stats"]["probes"] == 4  # all four candidates were asked about
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
