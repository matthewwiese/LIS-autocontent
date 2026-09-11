"""Tests for the whole-store datastore index.

The index is the model every output is built from, so these exercise it directly:
collection discovery, CHECKSUM parsing, file-name parsing, prediction and probing.
What the catalog document makes of it is covered in test_catalog.py.
"""

import os

import pytest

from lis_autocontent.datastore_files import (
    INDEX_SUFFIXES,
    NODE_README_FIELDS,
    README_FIELDS,
    DatastoreFiles,
    DatastoreIndex,
    split_pairwise_parents,
)

# Collections in the `metadata_dir` fixture tree (conftest.py).
GENOME = "Glycine/max/genomes/Wm82.gnm4.4PTR"
ANNOTATION = "Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"
QTL = "Glycine/max/qtl/Demo.qtl.X_1990"
MD5 = "d41d8cd98f00b204e9800998ecf8427e"


def _checksum(*names):
    return "\n".join(f"{MD5}  ./{name}" for name in names)


# --- discovery -----------------------------------------------------------------------
def test_every_collection_is_indexed_not_only_those_with_a_checksum(metadata_dir):
    """Keying off CHECKSUM alone missed ~45% of the store."""
    index = DatastoreIndex(metadata_dir).build()
    assert {path: c.index_status for path, c in index.collections.items()} == {
        ANNOTATION: "known",
        GENOME: "unknown",  # no CHECKSUM, and genomes have no documented vocabulary
        QTL: "inferred",
    }


def test_orphan_collections_are_those_without_a_checksum(metadata_dir):
    assert DatastoreIndex(metadata_dir).build().orphan_collections() == [GENOME, QTL]


def test_rebuilding_does_not_accumulate(metadata_dir):
    index = DatastoreIndex(metadata_dir)
    first = dict(index.build().stats)
    assert index.build().stats == first
    assert len(index.files) == sum(len(c.files) for c in index.collections.values())


def test_node_fields_are_a_subset_of_catalog_fields():
    """The node list is deliberately narrower, but must never drift outside it."""
    assert set(NODE_README_FIELDS) <= set(README_FIELDS)


# --- CHECKSUM parsing ------------------------------------------------------------------
def test_every_file_is_indexed_including_metadata_and_index_siblings(metadata_dir):
    """The index keeps every CHECKSUM line; `data_files` is the filtered view."""
    index = DatastoreIndex(metadata_dir).build()
    annotation = index.collections[ANNOTATION]
    assert len(annotation.files) == 9
    assert [r.filename for r in annotation.data_files] == [
        "glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz",
        "glyma.Wm82.gnm4.ann1.T8TQ.legume.fam3.VLMQ.gfa.tsv.gz",
        "glyma.Wm82.gnm4.ann1.T8TQ.protein_primary.faa.gz",
    ]
    protein = annotation.data_files[2]
    fai = index.companion(protein, ".fai")
    assert fai is not None and fai.md5 == MD5
    assert protein.indexes == [".fai"]


def test_file_names_are_parsed_into_canonical_type_and_extensions(metadata_dir):
    index = DatastoreIndex(metadata_dir).build()
    (gff,) = index.select(canonical_type="gene_models_main", endswith=".gff3.gz")
    assert (gff.canonical_type, gff.extensions, gff.gensp) == (
        "gene_models_main",
        "gff3.gz",
        "glyma",
    )
    assert index.sibling(gff, "protein_primary").extensions == "faa.gz"


def test_a_checksum_not_named_for_its_directory_is_still_read(metadata_dir, write):
    """legume.fam1.M65K publishes CHECKSUM.mixed.fam1.M65K.md5. Looking only for
    CHECKSUM.<directory>.md5 missed it and demoted an authoritative list to a guess."""
    path = "LEGUMES/Fabaceae/genefamilies/legume.fam1.M65K"
    coll = os.path.join(metadata_dir, path)
    write(os.path.join(coll, "README.legume.fam1.M65K.yml"), "identifier: x\n")
    write(
        os.path.join(coll, "CHECKSUM.mixed.fam1.M65K.md5"),
        _checksum("legume.fam1.M65K.info_annot_ahrd.tsv.gz"),
    )
    collection = DatastoreIndex(metadata_dir).build().collections[path]
    assert collection.index_status == "known"
    assert [r.src for r in collection.data_files] == ["checksum"]


def test_every_checksum_in_a_collection_is_read_once(metadata_dir, write):
    """One diversity collection publishes two CHECKSUMs, and they overlap."""
    path = "Glycine/max/diversity/Wm82.gnm2.div.Demo"
    coll = os.path.join(metadata_dir, path)
    vcf = "glyma.Wm82.gnm2.div.Demo.vcf.gz"
    write(
        os.path.join(coll, "CHECKSUM.Wm82.gnm2.div.Demo.md5"),
        _checksum(vcf, vcf + ".tbi"),
    )
    write(
        os.path.join(coll, "CHECKSUM.glyma.Wm82.gnm2.div.Demo.md5"),
        _checksum(vcf, "CHANGES.Wm82.gnm2.div.Demo.txt"),
    )
    collection = DatastoreIndex(metadata_dir).build().collections[path]
    assert [r.relative_path for r in collection.files] == [
        "CHANGES.Wm82.gnm2.div.Demo.txt",
        vcf,
        vcf + ".tbi",
    ]
    assert collection.data_files[0].indexes == [".tbi"]


def test_hidden_file_names_keep_their_leading_dot(metadata_dir, write):
    """Only the `./` prefix is stripped: `./.nextflow.log` is `.nextflow.log`."""
    path = "Arachis/hypogaea/expression/Tifrunner.gnm2.ann2.expr.Demo"
    coll = os.path.join(metadata_dir, path)
    write(
        os.path.join(coll, "CHECKSUM.Tifrunner.gnm2.ann2.expr.Demo.md5"),
        _checksum(".nextflow.log", "BUSCO/.busco.full_table.tsv.gz"),
    )
    collection = DatastoreIndex(metadata_dir).build().collections[path]
    assert [r.relative_path for r in collection.files] == [
        ".nextflow.log",
        "BUSCO/.busco.full_table.tsv.gz",
    ]


@pytest.mark.parametrize(
    "name, metadata",
    [
        ("README.X.gnm1.ABCD.yml", True),
        ("README", True),
        ("CHANGES.X.gnm1.ABCD.txt.gz", True),
        ("BUSCO/MANIFEST.X.gnm1.ABCD.yml", True),
        ("MANIFEST-000002", False),  # six real data files are named like this
        ("glyma.X.gnm1.ABCD.genome_main.fna.gz", False),
    ],
)
def test_metadata_is_matched_on_the_dotted_prefix(name, metadata):
    record = DatastoreFiles("Glycine", "max", "genomes", "X.gnm1.ABCD", name)
    assert record.is_metadata is metadata


def test_pairwise_parents_are_parsed_from_the_file_name():
    assert split_pairwise_parents(
        "glyma.Wm82.gnm4.x.phavu.G19833.gnm2.PXV3.minimap2.paf.gz"
    ) == ["glyma.Wm82.gnm4", "phavu.G19833.gnm2"]
    assert split_pairwise_parents("glyma.Wm82.gnm4.4PTR.genome_main.fna.gz") == []


# --- BUSCO and counts, offline -----------------------------------------------------
def test_metrics_are_parsed_from_the_committed_summary(metadata_dir):
    """The short_summary JSON carries assembly metrics as well as completeness, so
    no .fai fetch is needed and the build stays offline."""
    index = DatastoreIndex(metadata_dir)
    busco, counts = index.busco_summary(os.path.join(metadata_dir, GENOME))
    assert busco["complete_pct"] == 99.4
    assert busco["total"] == 5366
    assert counts["contig_n50"] == 419290
    assert counts["scaffold_n50"] == 49893278
    assert counts["percent_gaps"] == 2.648
    assert index.build().collections[GENOME].busco == busco


def test_busco_summary_absent_is_not_an_error(metadata_dir):
    index = DatastoreIndex(metadata_dir)
    assert index.busco_summary(os.path.join(metadata_dir, ANNOTATION)) == (None, None)


# --- file prediction for collections that publish no CHECKSUM ------------------------
def test_predicted_files_follow_the_documented_convention(metadata_dir):
    """qtl/gwas/maps ship no CHECKSUM, so their files are constructed from the
    vocabulary in filetypes.yml plus `{abbrev}.{identifier}.{content}.{ext}`."""
    index = DatastoreIndex(metadata_dir)
    names = index.predicted_files("qtl", "Demo.qtl.X_1990", "glyma")
    assert names == [
        "glyma.Demo.qtl.X_1990.obo.tsv.gz",
        "glyma.Demo.qtl.X_1990.qtl.tsv.gz",
        "glyma.Demo.qtl.X_1990.qtlmrk.tsv.gz",
        "glyma.Demo.qtl.X_1990.trait.tsv.gz",
    ]


def test_prefixless_types_omit_the_abbrev(metadata_dir):
    """pangenes and genefamilies are genus/family-scoped: files are
    `Cicer.pan2.CMWZ.clust.tsv.gz`, with no abbrev in front."""
    index = DatastoreIndex(metadata_dir)
    names = index.predicted_files("pangenes", "Cicer.pan2.CMWZ", "cicar")
    assert all(n.startswith("Cicer.pan2.CMWZ.") for n in names)


def test_unpredictable_types_are_never_guessed(metadata_dir):
    """Supplementary files carry study-specific names with no convention."""
    index = DatastoreIndex(metadata_dir)
    assert index.predicted_files("supplements", "mixed.esm.X_2016", "glyma") == []


def test_unknown_type_predicts_nothing(metadata_dir):
    assert (
        DatastoreIndex(metadata_dir).predicted_files("nosuchtype", "X", "glyma") == []
    )


def test_abbrev_prefixed_type_without_an_abbrev_predicts_nothing(metadata_dir):
    """Better an empty list than a filename missing its prefix."""
    index = DatastoreIndex(metadata_dir)
    assert index.predicted_files("qtl", "Demo.qtl.X_1990", None) == []


def test_a_missing_vocabulary_degrades_to_no_prediction(metadata_dir, monkeypatch):
    """Losing filetypes.yml must not fail a build; it just stops predicting."""
    index = DatastoreIndex(metadata_dir)
    monkeypatch.setattr(index, "filetypes", {})
    assert index.predicted_files("qtl", "Demo.qtl.X_1990", "glyma") == []


def test_types_documented_as_unindexed_are_not_probed_for_indexes(metadata_dir):
    """qtl/gwas/maps never publish .fai/.tbi (verified: 72 probes, 0 hits), so their
    files carry no index-unknown marker."""
    index = DatastoreIndex(metadata_dir)
    assert index.type_is_indexed("qtl") is False
    assert index.type_is_indexed("annotations") is True


# --- probing ---------------------------------------------------------------------------
def test_verify_probes_percent_encode_non_ascii_paths(metadata_dir, monkeypatch):
    """21 collections carry non-ASCII names. Unencoded, urllib raises rather than
    returning a status — which, swallowed, reports every one of their files as absent.
    This regressed once and cost 17 collections their entire file list."""
    asked = []

    def fake_urlopen(request, timeout=None):  # pylint: disable=unused-argument
        asked.append(request.full_url)
        raise OSError("stop here")

    monkeypatch.setattr("urllib.request.urlopen", fake_urlopen)
    DatastoreIndex(metadata_dir).url_exists("https://x/Nicolás_Hungria_2005/f.tsv.gz")
    assert asked and "Nicol%C3%A1s" in asked[0]
    assert "á" not in asked[0]


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
    assert DatastoreIndex.plausible_indexes(name) == expected


def test_probed_index_siblings_are_indexed_as_files(
    metadata_dir, markers_collection, monkeypatch
):
    """A sibling found by a HEAD request is a real file: `companion` finds it, and it
    is labelled as verified rather than passed off as CHECKSUM-listed."""
    real = {"glyma.SoySNP50K.mrk.ABCD.gff3.gz", "glyma.SoySNP50K.mrk.ABCD.gff3.gz.tbi"}
    index = DatastoreIndex(metadata_dir, verify=True)
    monkeypatch.setattr(index, "url_exists", lambda url: url.rsplit("/", 1)[-1] in real)
    index.build()
    (gff,) = index.collections[markers_collection].data_files
    tbi = index.companion(gff, ".tbi")
    assert (tbi.src, tbi.md5) == ("verified", None)
    assert index.companion(gff, ".csi") is None
