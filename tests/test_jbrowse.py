"""The JBrowse plan: what populate-jbrowse2 asks the jbrowse CLI to configure."""

import json
import os

import pytest
import yaml

from lis_autocontent import jbrowse
from lis_autocontent.catalog import CatalogBuilder
from lis_autocontent.datastore_files import DatastoreIndex

MD5 = "d41d8cd98f00b204e9800998ecf8427e"
ANNOTATION = "Cameor.gnm2.ann1.KT6X"
GFF3 = f"pissa.{ANNOTATION}.gene_models_main.gff3.gz"


GENOME = "Cameor.gnm2.ABCD"
FASTA = f"pissa.{GENOME}.genome_main.fna.gz"
DATA = "https://data.legumeinfo.org/Pisum/sativum"


def _tree(write, root, annotation_indexes):
    """One pea genome and its annotation, the annotation indexed as given."""
    write(
        os.path.join(root, "Pisum/GENUS/about_this_collection/description_Pisum.yml"),
        yaml.safe_dump({"genus": "Pisum", "species": ["sativum"]}),
    )
    for ctype, key, names in (
        ("genomes", GENOME, (FASTA, FASTA + ".fai", FASTA + ".gzi")),
        ("annotations", ANNOTATION, (GFF3, *(GFF3 + s for s in annotation_indexes))),
    ):
        base = os.path.join(root, f"Pisum/sativum/{ctype}/{key}")
        write(
            os.path.join(base, f"README.{key}.yml"), yaml.safe_dump({"identifier": key})
        )
        write(
            os.path.join(base, f"CHECKSUM.{key}.md5"),
            "".join(f"{MD5}  ./{name}\n" for name in names),
        )


def _uri(url):
    return {"uri": url, "locationType": "UriLocation"}


def _config(path):
    """A deployed config: the genome, its annotation behind a .tbi, a portal-hosted
    track and a track in a collection the metadata doesn't have."""
    gff3 = f"{DATA}/annotations/{ANNOTATION}/{GFF3}"
    config = {
        "assemblies": [
            {
                "name": "pissa.Cameor.gnm2",
                "sequence": {
                    "adapter": {
                        "type": "BgzipFastaAdapter",
                        "fastaLocation": _uri(f"{DATA}/genomes/{GENOME}/{FASTA}"),
                        "faiLocation": _uri(f"{DATA}/genomes/{GENOME}/{FASTA}.fai"),
                    }
                },
            }
        ],
        "tracks": [
            {
                "type": "FeatureTrack",
                "trackId": GFF3[: -len(".gz")],
                "assemblyNames": ["pissa.Cameor.gnm2"],
                "adapter": {
                    "type": "Gff3TabixAdapter",
                    "gffGzLocation": _uri(gff3),
                    "index": {"location": _uri(f"{gff3}.tbi"), "indexType": "TBI"},
                },
            },
            {
                "type": "QuantitativeTrack",
                "trackId": "portal_rnaseq",
                "assemblyNames": ["pissa.Cameor.gnm2"],
                "adapter": {"type": "BigWigAdapter", "bigWigLocation": _uri("x.bw")},
            },
            {
                "type": "SyntenyTrack",
                "trackId": "moved.paf",
                "assemblyNames": ["pissa.Cameor.gnm2"],
                "adapter": {
                    "type": "PAFAdapter",
                    "pafLocation": _uri(f"{DATA}/synteny/Gone.syn.ZZZZ/moved.paf.gz"),
                },
            },
        ],
    }
    path.write_text(json.dumps(config), encoding="utf-8")
    return str(path)


def test_placements_come_from_the_deployed_config(tmp_path, write):
    root = str(tmp_path / "datastore-metadata")
    _tree(write, root, (".csi",))
    builder = CatalogBuilder(root, jbrowse_config=_config(tmp_path / "config.json"))
    document = builder.build()
    (instance,) = document["jbrowse_instances"].items()
    assert instance[0] == "all-genera" and instance[1]["status"] == "ok"
    assert instance[1]["url"] == "https://all-genera.lis.ncgr.org/tools/jbrowse2/"
    records = {c["id"]: c for c in document["collections"]}
    assert records[GENOME]["jbrowse"] == [
        {"instance": "all-genera", "assemblies": ["pissa.Cameor.gnm2"]}
    ]
    assert records[ANNOTATION]["jbrowse"] == [
        {
            "instance": "all-genera",
            "assemblies": ["pissa.Cameor.gnm2"],
            "tracks": [
                {
                    "id": GFF3[: -len(".gz")],
                    "type": "FeatureTrack",
                    "file": GFF3,
                    "index": "TBI",
                }
            ],
        }
    ]


def test_an_unreadable_config_leaves_its_instance_unavailable(tmp_path, write):
    """Unknown, not absent: no placements are claimed for an instance not read."""
    root = str(tmp_path / "datastore-metadata")
    _tree(write, root, (".csi",))
    broken = tmp_path / "config.json"
    broken.write_text("{not json", encoding="utf-8")
    document = CatalogBuilder(root, jbrowse_config=str(broken)).build()
    assert document["jbrowse_instances"]["all-genera"]["status"] == "unavailable"
    assert not any("jbrowse" in c for c in document["collections"])


def test_the_report_names_every_drift(tmp_path, write):
    root = str(tmp_path / "datastore-metadata")
    _tree(write, root, (".csi",))
    builder = CatalogBuilder(root, jbrowse_config=_config(tmp_path / "config.json"))
    builder.build()
    index = builder.index
    text = jbrowse.report(
        index,
        builder.deployments,
        jbrowse.plan(index, ["Pisum"]),
    )
    assert "| all-genera | ok | 1 | 1 of 3 | 1 | 1 | 0 | 1 |" in text
    assert f"{DATA}/annotations/{ANNOTATION}/{GFF3}.tbi" in text  # CSI-only collection
    assert f"{DATA}/synteny/Gone.syn.ZZZZ/moved.paf.gz" in text
    assert "deployed index type differs from the plan's: 1" in text
    assert "Assemblies planned, not deployed: 0" in text


@pytest.mark.parametrize(
    "indexes, names_csi",
    [((".csi",), True), ((".tbi",), False), ((".tbi", ".csi"), False)],
)
def test_an_annotation_published_with_csi_only_names_its_index(
    tmp_path, write, indexes, names_csi
):
    """Without --indexFile the CLI configures <file>.tbi, which 404s for CSI-only."""
    root = str(tmp_path / "datastore-metadata")
    _tree(write, root, indexes)
    entries = jbrowse.plan(DatastoreIndex(root).build(), ["Pisum"])
    entry = next(e for e in entries if e.kind == "annotation")
    command = jbrowse.command(entry, str(tmp_path / "jbrowse"))
    assert (f"--indexFile {entry.url}.csi " in command) is names_csi
    assert f"--trackId {GFF3[: -len('.gz')]} " in command
