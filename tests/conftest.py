"""Shared fixtures: a miniature datastore-metadata checkout.

The tree exercises the four metadata layers, the index-status states, lineage linking
and field inheritance. No network is touched, which is the property an offline build is
supposed to guarantee.
"""

import json
import os

import pytest

MD5 = "d41d8cd98f00b204e9800998ecf8427e"

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

MARKERS_README = """---
identifier: SoySNP50K.mrk.ABCD
synopsis: Demo marker set
scientific_name: Glycine max
scientific_name_abbrev: glyma
"""

ANNOTATION_CHECKSUM = "\n".join(
    f"{MD5}  ./{name}"
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

GENOME = "Glycine/max/genomes/Wm82.gnm4.4PTR"
ANNOTATION = "Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"
QTL = "Glycine/max/qtl/Demo.qtl.X_1990"
MARKERS = "Glycine/max/markers/SoySNP50K.mrk.ABCD"


def write_file(path, text):
    """Write ``text`` to ``path``, creating its directories."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


@pytest.fixture(name="write")
def fixture_write():
    """The file writer, for tests that add to the tree."""
    return write_file


@pytest.fixture(name="metadata_dir")
def fixture_metadata_dir(tmp_path):
    """A miniature datastore-metadata checkout."""
    root = str(tmp_path)
    genome = os.path.join(root, GENOME)
    annotation = os.path.join(root, ANNOTATION)
    qtl = os.path.join(root, QTL)

    write_file(os.path.join(genome, "README.Wm82.gnm4.4PTR.yml"), GENOME_README)
    write_file(
        os.path.join(genome, "BUSCO", "glyma.Wm82.gnm4.4PTR.busco.short_summary.json"),
        json.dumps(BUSCO),
    )
    write_file(
        os.path.join(annotation, "README.Wm82.gnm4.ann1.T8TQ.yml"), ANNOTATION_README
    )
    write_file(
        os.path.join(annotation, "CHECKSUM.Wm82.gnm4.ann1.T8TQ.md5"),
        ANNOTATION_CHECKSUM,
    )
    write_file(
        os.path.join(annotation, "MANIFEST.Wm82.gnm4.ann1.T8TQ.yml"),
        ANNOTATION_MANIFEST,
    )
    # No CHECKSUM: the real shape of nearly every qtl and gwas collection.
    write_file(os.path.join(qtl, "README.Demo.qtl.X_1990.yml"), QTL_README)
    return root


@pytest.fixture(name="markers_collection")
def fixture_markers_collection(metadata_dir):
    """Add a markers collection -- no CHECKSUM, and a type that CAN carry index
    siblings -- and return its path."""
    write_file(
        os.path.join(metadata_dir, MARKERS, "README.SoySNP50K.mrk.ABCD.yml"),
        MARKERS_README,
    )
    return MARKERS
