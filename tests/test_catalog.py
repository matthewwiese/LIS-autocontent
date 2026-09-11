"""Tests for the offline catalog builder.

A miniature datastore-metadata tree exercises the four metadata layers, the
tri-state index status, lineage linking and field inheritance. No network is
touched, which is the property the catalog build is supposed to guarantee.
"""

import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from catalog import (  # noqa: E402  pylint: disable=wrong-import-position
    INDEX_SUFFIXES,
    NODE_README_FIELDS,
    README_FIELDS,
    SCHEMA_VERSION,
    CatalogBuilder,
)

GENOME_README = """---
identifier: Wm82.gnm4.4PTR
synopsis: Glycine max Williams 82 genome assembly v4.0
scientific_name: Glycine max
taxid: 3847
scientific_name_abbrev: glyma
bioproject: PRJNA19861
chromosome_prefix: Gm
supercontig_prefix: scaffold
genotype:
  - Williams 82
publication_doi: 10.1111/tpj.14500
publication_title: "Three reference-quality genome assemblies for soybean"
license: Open
"""

ANNOTATION_README = """---
identifier: Wm82.gnm4.ann1.T8TQ
synopsis: Glycine max Williams 82 annotation
scientific_name: Glycine max
taxid: 3847
scientific_name_abbrev: glyma
publication_doi: 10.1111/tpj.14500
license: Open
"""

QTL_README = """---
identifier: Demo.qtl.X_1990
synopsis: Demo QTL study
scientific_name: Glycine max
taxid: 3847
scientific_name_abbrev: glyma
publication_doi: 10.1007/bf00226154
genetic_map: GmComposite1999
"""

ANNOTATION_CHECKSUM = "\n".join(
    "d41d8cd98f00b204e9800998ecf8427e  ./{0}".format(name)
    for name in [
        "README.Wm82.gnm4.ann1.T8TQ.yml",
        "MANIFEST.Wm82.gnm4.ann1.T8TQ.yml",
        "CHECKSUM.Wm82.gnm4.ann1.T8TQ.md5",
        "glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz",
        "glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz.fai",
        "glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz.gzi",
        "glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz",
        "glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz.tbi",
        "glyma.Wm82.gnm4.ann1.T8TQ.legume.fam3.VLMQ.gfa.tsv.gz",
    ]
)

ANNOTATION_MANIFEST = """---
- name: glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz
  description: Protein sequences - primary only
  applications:
   - blast
   - mines
- name: glyma.Wm82.gnm4.ann1.T8TQ.legume.fam3.VLMQ.gfa.tsv.gz
  description: MISSING
"""

BUSCO = {
    "results": {
        "Complete": 99.4,
        "Single copy": 36.9,
        "Multi copy": 62.5,
        "Fragmented": 0.1,
        "Missing": 0.5,
        "n_markers": 5366,
        "domain": "eukaryota",
        "Number of scaffolds": "282",
        "Number of contigs": "9200",
        "Total length": "978386919",
        "Percent gaps": "2.648%",
        "Scaffold N50": "49893278",
        "Contigs N50": "419290",
    }
}


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


@pytest.fixture(name="metadata_dir")
def fixture_metadata_dir(tmp_path):
    """A miniature datastore-metadata checkout."""
    root = str(tmp_path)
    genome = os.path.join(root, "Glycine", "max", "genomes", "Wm82.gnm4.4PTR")
    annotation = os.path.join(
        root, "Glycine", "max", "annotations", "Wm82.gnm4.ann1.T8TQ"
    )
    qtl = os.path.join(root, "Glycine", "max", "qtl", "Demo.qtl.X_1990")

    _write(os.path.join(genome, "README.Wm82.gnm4.4PTR.yml"), GENOME_README)
    _write(
        os.path.join(genome, "BUSCO", "glyma.Wm82.gnm4.4PTR.busco.short_summary.json"),
        json.dumps(BUSCO),
    )
    _write(
        os.path.join(annotation, "README.Wm82.gnm4.ann1.T8TQ.yml"), ANNOTATION_README
    )
    _write(
        os.path.join(annotation, "CHECKSUM.Wm82.gnm4.ann1.T8TQ.md5"),
        ANNOTATION_CHECKSUM,
    )
    _write(
        os.path.join(annotation, "MANIFEST.Wm82.gnm4.ann1.T8TQ.yml"),
        ANNOTATION_MANIFEST,
    )
    # No CHECKSUM: the real shape of nearly every qtl and gwas collection.
    _write(os.path.join(qtl, "README.Demo.qtl.X_1990.yml"), QTL_README)
    return root


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
    assert sorted(_by_path(catalog)) == [
        "Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ",
        "Glycine/max/genomes/Wm82.gnm4.4PTR",
        "Glycine/max/qtl/Demo.qtl.X_1990",
    ]


# --- README propagation ------------------------------------------------------
def test_readme_fields_are_carried_verbatim(catalog):
    """The regression this exists to prevent: taxid and publication_doi were both
    parsed and then dropped before reaching a consumer."""
    genome = _by_path(catalog)["Glycine/max/genomes/Wm82.gnm4.4PTR"]
    assert genome["taxid"] == 3847
    assert genome["publication_doi"] == "10.1111/tpj.14500"
    assert genome["scientific_name_abbrev"] == "glyma"
    assert genome["license"] == "Open"
    assert genome["genotype"] == ["Williams 82"]


def test_node_fields_are_a_subset_of_catalog_fields():
    """The node list is deliberately narrower, but must never drift outside it."""
    assert set(NODE_README_FIELDS) <= set(README_FIELDS)


def test_qtl_keeps_its_metadata_without_a_checksum(catalog):
    """No CHECKSUM costs the file list, not the metadata."""
    qtl = _by_path(catalog)["Glycine/max/qtl/Demo.qtl.X_1990"]
    assert qtl["publication_doi"] == "10.1007/bf00226154"
    assert qtl["genetic_map"] == "GmComposite1999"


# --- index status: how the file list was obtained -----------------------------
def test_index_status_distinguishes_how_the_file_list_was_obtained(catalog):
    """Four states, and consumers rely on the difference: `known` came from a CHECKSUM,
    `inferred`/`verified` were constructed from the documented convention, `unknown`
    means neither was possible. Collapsing them presents a guess as a fact."""
    qtl = _by_path(catalog)["Glycine/max/qtl/Demo.qtl.X_1990"]
    assert qtl["index_status"] == "inferred"  # no CHECKSUM, vocabulary applied
    assert {f["src"] for f in qtl["files"]} == {"predicted"}

    annotation = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    assert annotation["index_status"] == "known"  # CHECKSUM published


def test_index_siblings_become_flags_not_entries(catalog):
    """.fai/.tbi/.gzi are evidence about other files, never data files themselves."""
    annotation = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    names = [f["n"] for f in annotation["files"]]
    assert len(names) == 3  # protein, gene models, gene families
    assert not any(n.endswith((".fai", ".tbi", ".gzi")) for n in names)
    indexed = {f["n"]: f["i"] for f in annotation["files"] if "i" in f}
    assert indexed == {
        "glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz": [".fai"],
        "glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz": [".tbi"],
    }


def test_metadata_files_are_not_listed_as_data(catalog):
    annotation = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    names = [f["n"] for f in annotation["files"]]
    assert not any(n.startswith(("README.", "MANIFEST.", "CHECKSUM.")) for n in names)


# --- MANIFEST ----------------------------------------------------------------
def test_manifest_descriptions_and_applications_merge_onto_files(catalog):
    annotation = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    protein = [f for f in annotation["files"] if "protein_primary" in f["n"]][0]
    assert protein["description"] == "Protein sequences - primary only"
    assert protein["applications"] == ["blast", "mines"]


def test_manifest_placeholder_description_is_dropped(catalog):
    """ "MISSING" is the spec's placeholder, not a description."""
    annotation = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    families = [f for f in annotation["files"] if "gfa.tsv" in f["n"]][0]
    assert "description" not in families


# --- BUSCO and counts, offline ----------------------------------------------
def test_metrics_are_parsed_from_the_committed_summary(metadata_dir):
    """The parser stays live even though emission is deferred.

    The short_summary JSON carries assembly metrics as well as completeness, so
    no .fai fetch is needed and the build stays offline. Re-enabling the metrics
    is uncommenting their inclusion, not rewriting this.
    """
    builder = CatalogBuilder(metadata_dir)
    genome_dir = os.path.join(
        metadata_dir, "Glycine", "max", "genomes", "Wm82.gnm4.4PTR"
    )
    busco, counts = builder.busco_summary(genome_dir)
    assert busco["complete_pct"] == 99.4
    assert busco["total"] == 5366
    assert counts["contig_n50"] == 419290
    assert counts["scaffold_n50"] == 49893278
    assert counts["percent_gaps"] == 2.648


def test_metrics_are_not_emitted_yet(catalog):
    """Deliberately deferred: see the NOTE in catalog.py. This test exists so the
    omission is a decision on record rather than something that quietly regresses
    in either direction."""
    genome = _by_path(catalog)["Glycine/max/genomes/Wm82.gnm4.4PTR"]
    assert "busco" not in genome
    assert "counts" not in genome


def test_busco_availability_is_still_reported_in_stats(catalog):
    """The build diagnostic must stay truthful while the metrics are withheld."""
    assert catalog["stats"]["with_busco"] == 1


def test_busco_summary_absent_is_not_an_error(metadata_dir):
    builder = CatalogBuilder(metadata_dir)
    annotation = os.path.join(
        metadata_dir, "Glycine", "max", "annotations", "Wm82.gnm4.ann1.T8TQ"
    )
    assert builder.busco_summary(annotation) == (None, None)


# --- lineage and inheritance -------------------------------------------------
def test_annotation_links_to_its_genome(catalog):
    annotation = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    assert annotation["derived_from"] == ["Wm82.gnm4.4PTR"]


def test_assembly_conventions_are_inherited_down_the_edge(catalog):
    """chromosome_prefix lives on the genome README only. Without inheritance a
    caller cannot turn `Gm12` into `glyma.Wm82.gnm4.Gm12` from an annotation."""
    annotation = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    assert annotation["chromosome_prefix"] == "Gm"
    assert annotation["bioproject"] == "PRJNA19861"
    assert "chromosome_prefix" in annotation["inherited"]


def test_inheritance_never_overwrites_a_published_value(tmp_path):
    """An annotation that publishes its own prefix keeps it."""
    root = str(tmp_path)
    genome = os.path.join(root, "Glycine", "max", "genomes", "Wm82.gnm4.4PTR")
    annotation = os.path.join(
        root, "Glycine", "max", "annotations", "Wm82.gnm4.ann1.T8TQ"
    )
    _write(os.path.join(genome, "README.Wm82.gnm4.4PTR.yml"), GENOME_README)
    _write(
        os.path.join(annotation, "README.Wm82.gnm4.ann1.T8TQ.yml"),
        ANNOTATION_README + "chromosome_prefix: Chr\n",
    )
    catalog = CatalogBuilder(root).build()
    record = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    assert record["chromosome_prefix"] == "Chr"
    assert "chromosome_prefix" not in record.get("inherited", [])


def test_genomes_do_not_derive_from_themselves(catalog):
    genome = _by_path(catalog)["Glycine/max/genomes/Wm82.gnm4.4PTR"]
    assert "derived_from" not in genome


# --- urls and output ---------------------------------------------------------
def test_base_url_points_at_the_datastore(catalog):
    annotation = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    assert annotation["base_url"] == (
        "https://data.legumeinfo.org/Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"
    )


def test_datastore_url_is_overridable(metadata_dir):
    catalog = CatalogBuilder(
        metadata_dir, datastore_url="https://example.org/ds/"
    ).build()
    record = _by_path(catalog)["Glycine/max/genomes/Wm82.gnm4.4PTR"]
    assert record["base_url"].startswith("https://example.org/ds/Glycine")


def test_write_round_trips(metadata_dir, tmp_path):
    out = str(tmp_path / "out" / "catalog.json")
    CatalogBuilder(metadata_dir).write(out)
    with open(out, encoding="utf-8") as handle:
        assert json.load(handle)["stats"]["collections"] == 3


def test_unreadable_readme_degrades_the_record_not_the_build(metadata_dir):
    """One malformed file must not sink a 1,000-collection build — and must not drop
    the collection either. A real datastore collection ships a README that is not
    valid YAML; its files are still real and still worth cataloguing."""
    broken = os.path.join(metadata_dir, "Glycine", "max", "maps", "Bad.map.X")
    _write(os.path.join(broken, "README.Bad.map.X.yml"), "---\n: : not yaml : :\n")
    catalog = CatalogBuilder(metadata_dir).build()
    assert catalog["stats"]["collections"] == 4
    record = _by_path(catalog)["Glycine/max/maps/Bad.map.X"]
    assert record["index_status"] == "unknown"  # no CHECKSUM either
    assert "publication_doi" not in record  # nothing was invented


def test_a_collection_with_only_a_checksum_is_catalogued(metadata_dir):
    """Five genome_alignments collections publish a CHECKSUM and no README. Keying off
    the README alone dropped them, along with the indexed BAMs they hold."""
    wga = os.path.join(
        metadata_dir, "Arachis", "hypogaea", "genome_alignments", "X.wga.AAAA"
    )
    _write(
        os.path.join(wga, "CHECKSUM.X.wga.AAAA.md5"),
        "\n".join(
            [
                "d41d8cd98f00b204e9800998ecf8427e  ./arahy.x.y.AAAA.bam",
                "d41d8cd98f00b204e9800998ecf8427e  ./arahy.x.y.AAAA.bam.bai",
            ]
        ),
    )
    catalog = CatalogBuilder(metadata_dir).build()
    record = _by_path(catalog)["Arachis/hypogaea/genome_alignments/X.wga.AAAA"]
    assert record["index_status"] == "known"
    assert record["files"] == [
        {"n": "arahy.x.y.AAAA.bam", "i": [".bai"], "src": "checksum"}
    ]


def test_a_bare_README_is_still_parsed(metadata_dir):
    """At least one collection ships `README` with no extension; it is still YAML."""
    wga = os.path.join(
        metadata_dir, "Arachis", "ipaensis", "genome_alignments", "Y.wga.BBBB"
    )
    _write(
        os.path.join(wga, "README"),
        "identifier: Y.wga.BBBB\nsynopsis: Genome alignments\npublication_doi: 10.1/x\n",
    )
    catalog = CatalogBuilder(metadata_dir).build()
    record = _by_path(catalog)["Arachis/ipaensis/genome_alignments/Y.wga.BBBB"]
    assert record["synopsis"] == "Genome alignments"
    assert record["publication_doi"] == "10.1/x"


def test_taxa_descriptions_are_ingested(metadata_dir):
    """Common name and abbreviation live only in about_this_collection, so without
    this a consumer cannot map 'soybean' onto a taxon at all."""
    about = os.path.join(metadata_dir, "Glycine", "max", "about_this_collection")
    _write(
        os.path.join(about, "description_Glycine_max.yml"),
        "---\ntaxid: 3847\ngenus: Glycine\nspecies: max\nabbrev: glyma\n"
        "commonName: soybean\nresources:\n  - name: GlycineMine\n    URL: https://x\n",
    )
    catalog = CatalogBuilder(metadata_dir).build()
    entry = catalog["taxa"]["Glycine/max"]
    assert entry["commonName"] == "soybean"
    assert entry["abbrev"] == "glyma"
    assert entry["resources"][0]["name"] == "GlycineMine"


def test_genus_level_description_lists_its_species(metadata_dir):
    """A genus file uses `species` as a LIST; keying it like a species file would
    conflate the two."""
    about = os.path.join(metadata_dir, "Glycine", "GENUS", "about_this_collection")
    _write(
        os.path.join(about, "description_Glycine.yml"),
        "---\ntaxid: 3847\ngenus: Glycine\ncommonName: soybean\n"
        "species:\n  - max\n  - soja\n",
    )
    catalog = CatalogBuilder(metadata_dir).build()
    assert catalog["taxa"]["Glycine"]["species_in_genus"] == ["max", "soja"]
    assert "species" not in catalog["taxa"]["Glycine"]


# --- file resolution for collections that publish no CHECKSUM -------------------------
def test_predicted_files_follow_the_documented_convention(metadata_dir):
    """qtl/gwas/maps ship no CHECKSUM, so their files are constructed from the
    vocabulary in filetypes.yml plus `{abbrev}.{identifier}.{content}.{ext}`."""
    builder = CatalogBuilder(metadata_dir)
    names = builder.predicted_files("qtl", "Demo.qtl.X_1990", "glyma")
    assert names == [
        "glyma.Demo.qtl.X_1990.obo.tsv.gz",
        "glyma.Demo.qtl.X_1990.qtl.tsv.gz",
        "glyma.Demo.qtl.X_1990.qtlmrk.tsv.gz",
        "glyma.Demo.qtl.X_1990.trait.tsv.gz",
    ]


def test_prefixless_types_omit_the_abbrev(metadata_dir):
    """pangenes and genefamilies are genus/family-scoped: files are
    `Cicer.pan2.CMWZ.clust.tsv.gz`, with no abbrev in front."""
    builder = CatalogBuilder(metadata_dir)
    names = builder.predicted_files("pangenes", "Cicer.pan2.CMWZ", "cicar")
    assert all(n.startswith("Cicer.pan2.CMWZ.") for n in names)


def test_unpredictable_types_are_never_guessed(metadata_dir):
    """Supplementary files carry study-specific names with no convention."""
    assert (
        CatalogBuilder(metadata_dir).predicted_files(
            "supplements", "mixed.esm.X_2016", "glyma"
        )
        == []
    )


def test_unknown_type_predicts_nothing(metadata_dir):
    assert (
        CatalogBuilder(metadata_dir).predicted_files("nosuchtype", "X", "glyma") == []
    )


def test_abbrev_prefixed_type_without_an_abbrev_predicts_nothing(metadata_dir):
    """Better an empty list than a filename missing its prefix."""
    assert (
        CatalogBuilder(metadata_dir).predicted_files("qtl", "Demo.qtl.X_1990", None)
        == []
    )


def test_offline_build_labels_files_as_predicted(catalog):
    """Without --verify nothing is probed, so the label must not claim more."""
    qtl = _by_path(catalog)["Glycine/max/qtl/Demo.qtl.X_1990"]
    assert qtl["index_status"] == "inferred"
    assert qtl["files"]
    assert {f["src"] for f in qtl["files"]} == {"predicted"}


def test_checksum_files_are_labelled_authoritative(catalog):
    annotation = _by_path(catalog)["Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"]
    assert annotation["index_status"] == "known"
    assert {f["src"] for f in annotation["files"]} == {"checksum"}


def test_verify_keeps_only_files_that_exist(metadata_dir, monkeypatch):
    """--verify turns prediction into evidence: absent files are dropped and the
    survivors are relabelled, so a consumer can tell the two apart."""
    real = {"glyma.Demo.qtl.X_1990.qtl.tsv.gz", "glyma.Demo.qtl.X_1990.obo.tsv.gz"}
    builder = CatalogBuilder(metadata_dir, verify=True)
    monkeypatch.setattr(
        builder, "url_exists", lambda url: url.rsplit("/", 1)[-1] in real
    )
    catalog = builder.build()
    qtl = _by_path(catalog)["Glycine/max/qtl/Demo.qtl.X_1990"]
    assert qtl["index_status"] == "verified"
    assert {f["n"] for f in qtl["files"]} == real
    assert {f["src"] for f in qtl["files"]} == {"verified"}
    assert catalog["stats"]["probes"] == 4  # all four candidates were asked about
    assert catalog["stats"]["verified_files"] == 2


def test_verify_probes_percent_encode_non_ascii_paths(metadata_dir, monkeypatch):
    """21 collections carry non-ASCII names. Unencoded, urllib raises rather than
    returning a status — which, swallowed, reports every one of their files as absent.
    This regressed once and cost 17 collections their entire file list."""
    asked = []

    def fake_urlopen(request, timeout=None):
        asked.append(request.full_url)
        raise OSError("stop here")

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    CatalogBuilder(metadata_dir).url_exists(
        "https://x/Nicol\u00e1s_Hungria_2005/f.tsv.gz"
    )
    assert asked and "Nicol%C3%A1s" in asked[0]
    assert "\u00e1" not in asked[0]


def test_types_documented_as_unindexed_are_not_probed_for_indexes(metadata_dir):
    """qtl/gwas/maps never publish .fai/.tbi (verified: 72 probes, 0 hits), so their
    files carry no index-unknown marker."""
    builder = CatalogBuilder(metadata_dir)
    assert builder.type_is_indexed("qtl") is False
    assert builder.type_is_indexed("annotations") is True


MARKERS_README = """---
identifier: SoySNP50K.mrk.ABCD
synopsis: Demo marker set
scientific_name: Glycine max
scientific_name_abbrev: glyma
"""


def _add_markers_collection(metadata_dir):
    """A markers collection: no CHECKSUM, and a type that CAN carry index siblings."""
    markers = os.path.join(
        metadata_dir, "Glycine", "max", "markers", "SoySNP50K.mrk.ABCD"
    )
    _write(os.path.join(markers, "README.SoySNP50K.mrk.ABCD.yml"), MARKERS_README)
    return "Glycine/max/markers/SoySNP50K.mrk.ABCD"


def test_offline_index_status_of_an_indexable_type_is_unknown(metadata_dir):
    """A predicted file of a type that can be indexed was never probed, so whether it is
    streamable is unknown -- not "no index"."""
    path = _add_markers_collection(metadata_dir)
    record = _by_path(CatalogBuilder(metadata_dir).build())[path]
    assert record["files"] == [
        {"n": "glyma.SoySNP50K.mrk.ABCD.gff3.gz", "src": "predicted", "i_unknown": True}
    ]


def test_verify_attaches_the_index_siblings_that_exist(metadata_dir, monkeypatch):
    """Under --verify the unknown is resolved: each surviving file is probed for the
    siblings its extension allows, and only the ones found are recorded."""
    path = _add_markers_collection(metadata_dir)
    real = {"glyma.SoySNP50K.mrk.ABCD.gff3.gz", "glyma.SoySNP50K.mrk.ABCD.gff3.gz.tbi"}
    builder = CatalogBuilder(metadata_dir, verify=True)
    monkeypatch.setattr(
        builder, "url_exists", lambda url: url.rsplit("/", 1)[-1] in real
    )
    catalog = builder.build()
    assert _by_path(catalog)[path]["files"] == [
        {"n": "glyma.SoySNP50K.mrk.ABCD.gff3.gz", "src": "verified", "i": [".tbi"]}
    ]
    # qtl: 4 file probes and no index probes (indexed: false); markers: 1 + .tbi/.csi.
    assert catalog["stats"]["probes"] == 7


@pytest.mark.parametrize(
    "name, expected",
    [
        ("x.genome_main.fna.gz", (".fai",)),
        ("x.protein.faa.gz", (".fai",)),
        ("x.bam", (".bai", ".csi")),
        ("x.cram", (".crai",)),
        ("x.gene_models_main.gff3.gz", (".tbi", ".csi")),
        ("x.vcf.gz", (".tbi", ".csi")),
        ("x.unrecognised.bin", INDEX_SUFFIXES),
    ],
)
def test_index_probes_are_limited_to_what_the_extension_allows(name, expected):
    """Asking a .tsv.gz about a .bai is a guaranteed 404. An unrecognised type falls
    back to asking about every sibling rather than none."""
    assert CatalogBuilder.plausible_indexes(name) == expected


def test_a_missing_vocabulary_degrades_to_no_prediction(metadata_dir, monkeypatch):
    """Losing filetypes.yml must not fail a build; it just stops predicting."""
    builder = CatalogBuilder(metadata_dir)
    monkeypatch.setattr(builder, "filetypes", {})
    assert builder.predicted_files("qtl", "Demo.qtl.X_1990", "glyma") == []
