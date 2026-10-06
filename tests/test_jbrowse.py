"""The JBrowse plan: what populate-jbrowse2 asks the jbrowse CLI to configure."""

import os

import pytest
import yaml

from lis_autocontent import jbrowse
from lis_autocontent.datastore_files import DatastoreIndex

MD5 = "d41d8cd98f00b204e9800998ecf8427e"
ANNOTATION = "Cameor.gnm2.ann1.KT6X"
GFF3 = f"pissa.{ANNOTATION}.gene_models_main.gff3.gz"


@pytest.mark.parametrize(
    "indexes, names_csi",
    [((".csi",), True), ((".tbi",), False), ((".tbi", ".csi"), False)],
)
def test_an_annotation_published_with_csi_only_names_its_index(
    tmp_path, write, indexes, names_csi
):
    """Without --indexFile the CLI configures <file>.tbi, which 404s for CSI-only."""
    root = str(tmp_path / "datastore-metadata")
    write(
        os.path.join(root, "Pisum/GENUS/about_this_collection/description_Pisum.yml"),
        yaml.safe_dump({"genus": "Pisum", "species": ["sativum"]}),
    )
    base = os.path.join(root, f"Pisum/sativum/annotations/{ANNOTATION}")
    write(
        os.path.join(base, f"README.{ANNOTATION}.yml"),
        yaml.safe_dump({"identifier": ANNOTATION}),
    )
    write(
        os.path.join(base, f"CHECKSUM.{ANNOTATION}.md5"),
        "".join(f"{MD5}  ./{name}\n" for name in (GFF3, *(GFF3 + s for s in indexes))),
    )
    (entry,) = jbrowse.plan(DatastoreIndex(root).build(), ["Pisum"])
    command = jbrowse.command(entry, str(tmp_path / "jbrowse"))
    assert (f"--indexFile {entry.url}.csi " in command) is names_csi
    assert f"--trackId {GFF3[: -len('.gz')]} " in command
