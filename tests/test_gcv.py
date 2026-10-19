"""populate-gcv: plan instances from a checkout, validate them, render their overlays.

The rendered tree is compared with tests/data/gcv/. After an intended change, run
`LIS_UPDATE_EXPECTED=1 pytest tests/test_gcv.py` and review the diff.
"""

import gzip
import json
import logging
import os

import pytest
import yaml

from lis_autocontent import gcv, gcv_inputs
from lis_autocontent.catalog import CatalogBuilder
from lis_autocontent.datastore_files import DatastoreIndex

EXPECTED = os.path.join(os.path.dirname(__file__), "data", "gcv")
UPDATING = bool(os.environ.get("LIS_UPDATE_EXPECTED"))
MD5 = "d41d8cd98f00b204e9800998ecf8427e"
URL = "https://data.legumeinfo.org"

SOY = "Glycine/max/annotations/Wm82.gnm4.ann1.T8TQ"
COWPEA_ANN1 = "Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann1.zb5D"
COWPEA_ANN2 = "Vigna/unguiculata/annotations/IT97K-499-35.gnm1.ann2.FD7K"


def _checksum(*names):
    return "".join(f"{MD5}  ./{name}\n" for name in names)


def _assembly(write, root, taxon, genome, annotations, fai=True):
    """A genome and its annotations, each with gene models and a legfed_v1_0 file."""
    species, abbrev = taxon
    base = os.path.join(root, species)
    fasta = f"{abbrev}.{genome}.genome_main.fna.gz"
    write(
        os.path.join(base, f"genomes/{genome}/CHECKSUM.{genome}.md5"),
        _checksum(fasta, *([fasta + ".fai"] if fai else [])),
    )
    for annotation in annotations:
        write(
            os.path.join(base, f"annotations/{annotation}/CHECKSUM.{annotation}.md5"),
            _checksum(
                f"{abbrev}.{annotation}.gene_models_main.gff3.gz",
                f"{abbrev}.{annotation}.legfed_v1_0.M65K.gfa.tsv.gz",
                f"{abbrev}.{annotation}.legume.fam3.VLMQ.gfa.tsv.gz",
            ),
        )


@pytest.fixture(name="checkout")
def fixture_checkout(tmp_path, write):
    """Soybean, and cowpea's two annotations of one assembly."""
    root = str(tmp_path / "datastore-metadata")
    _assembly(
        write, root, ("Glycine/max", "glyma"), "Wm82.gnm4.4PTR", ["Wm82.gnm4.ann1.T8TQ"]
    )
    _assembly(
        write,
        root,
        ("Vigna/unguiculata", "vigun"),
        "IT97K-499-35.gnm1.QnBW",
        ["IT97K-499-35.gnm1.ann1.zb5D", "IT97K-499-35.gnm1.ann2.FD7K"],
    )
    return root


def _config(tmp_path, instances, **top):
    """A gcv.yml with ``instances``, written and loaded."""
    config = {
        "template": "v2.7.0-c0",
        "env": {"DOMAIN": "localhost"},
        "instances": instances,
        **top,
    }
    path = tmp_path / "gcv.yml"
    path.write_text(yaml.safe_dump(config, sort_keys=False))
    gcv.load_config.cache_clear()
    return gcv.load_config(str(path))


def _instance(members, family="legfed_v1_0", profile="services", **env):
    return {"profile": profile, "family": family, "env": env, "members": members}


def _plan(checkout, config, names=None):
    return gcv.plan(DatastoreIndex(checkout, datastore_url=URL).build(), config, names)


def _problems(checkout, config):
    with pytest.raises(gcv.GcvError) as err:
        _plan(checkout, config)
    return str(err.value)


def test_a_member_loads_its_gene_models_families_and_genome_index(checkout, tmp_path):
    config = _config(tmp_path, {"glycine": _instance([SOY])})
    (instance,) = _plan(checkout, config)
    (load,) = instance.loads
    assert (load.genus, load.species, load.strain) == ("Glycine", "max", "Wm82")
    assert load.fai_url == (
        f"{URL}/Glycine/max/genomes/Wm82.gnm4.4PTR/"
        "glyma.Wm82.gnm4.4PTR.genome_main.fna.gz.fai"
    )
    assert (
        load.gene_gff_url
        == f"{URL}/{SOY}/glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3.gz"
    )
    assert (
        load.gfa_url
        == f"{URL}/{SOY}/glyma.Wm82.gnm4.ann1.T8TQ.legfed_v1_0.M65K.gfa.tsv.gz"
    )
    assert load.rename is None


def test_the_instance_family_picks_the_family_file(checkout, tmp_path):
    config = _config(tmp_path, {"glycine_fam3": _instance([SOY], family="legume.fam3")})
    (instance,) = _plan(checkout, config)
    assert instance.loads[0].gfa_url.endswith(".legume.fam3.VLMQ.gfa.tsv.gz")


def test_family_files_in_subdirectories_are_ignored(checkout, tmp_path, write):
    checksum = os.path.join(checkout, SOY, "CHECKSUM.Wm82.gnm4.ann1.T8TQ.md5")
    with open(checksum, encoding="utf-8") as handle:
        text = handle.read()
    write(
        checksum,
        text
        + _checksum("bad_gene/glyma.Wm82.gnm4.ann1.T8TQ.legfed_v1_0.M65K.gfa.tsv.gz"),
    )
    config = _config(tmp_path, {"glycine": _instance([SOY])})
    (instance,) = _plan(checkout, config)
    assert "/bad_gene/" not in instance.loads[0].gfa_url


def test_a_derived_gff3_beside_the_gene_models_is_not_taken_for_them(
    checkout, tmp_path, write
):
    checksum = os.path.join(checkout, SOY, "CHECKSUM.Wm82.gnm4.ann1.T8TQ.md5")
    with open(checksum, encoding="utf-8") as handle:
        text = handle.read()
    write(
        checksum,
        text
        + _checksum("glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gcv_genes.gff3.gz"),
    )
    config = _config(tmp_path, {"glycine": _instance([SOY])})
    (instance,) = _plan(checkout, config)
    assert instance.loads[0].gene_gff_url.endswith(".T8TQ.gene_models_main.gff3.gz")


def test_a_genome_whose_checksum_lacks_its_fai_is_refused(tmp_path, write):
    root = str(tmp_path / "datastore-metadata")
    _assembly(
        write,
        root,
        ("Glycine/max", "glyma"),
        "Wm82.gnm4.4PTR",
        ["Wm82.gnm4.ann1.T8TQ"],
        fai=False,
    )
    config = _config(tmp_path, {"glycine": _instance([SOY])})
    assert "lists no glyma.Wm82.gnm4.4PTR.genome_main.fna.gz.fai" in _problems(
        root, config
    )


def test_two_annotations_of_one_assembly_need_renamed_chromosomes(checkout, tmp_path):
    config = _config(tmp_path, {"lis": _instance([COWPEA_ANN1, COWPEA_ANN2])})
    assert "share assembly IT97K-499-35.gnm1.QnBW" in _problems(checkout, config)


def test_renamed_chromosomes_carry_the_annotation_version(checkout, tmp_path):
    members = [{"path": COWPEA_ANN1, "rename_chromosomes": True}, COWPEA_ANN2]
    config = _config(tmp_path, {"lis": _instance(members)})
    (instance,) = _plan(checkout, config)
    assert [load.rename for load in instance.loads] == [("gnm1", "gnm1.ann1"), None]


def test_every_problem_is_reported_at_once(checkout, tmp_path):
    members = [SOY, SOY, "Glycine/max/annotations/Absent.gnm1.ann1.XXXX"]
    config = _config(tmp_path, {"glycine": _instance(members, family="legfed_v9")})
    problems = _problems(checkout, config)
    assert "listed more than once" in problems
    assert "no such collection" in problems
    assert "0 legfed_v9 family files" in problems


def test_only_the_named_instances_are_planned(checkout, tmp_path):
    config = _config(
        tmp_path,
        {"glycine": _instance([SOY]), "lis": _instance([COWPEA_ANN2])},
    )
    assert [i.name for i in _plan(checkout, config, ["lis"])] == ["lis"]
    with pytest.raises(gcv.GcvError, match="no instance vigna"):
        _plan(checkout, config, ["vigna"])


def test_the_config_refuses_uppercase_names_and_derived_keys(tmp_path):
    with pytest.raises(gcv.GcvError, match="must be lowercase"):
        _config(tmp_path, {"Glycine": _instance([SOY])})
    for key in ("COMPOSE_PROJECT_NAME", "MICROSERVICES_PATH", "GCV_PATH", "RDB_URL"):
        with pytest.raises(gcv.GcvError, match=f"sets {key}"):
            _config(tmp_path, {"a": _instance([SOY])}, env={key: "x"})


def test_the_packaged_config_is_valid():
    assert gcv.load_config()["instances"]
    assert gcv.load_app_config()["servers"]


def test_an_instance_web_app_calls_its_own_microservices(checkout, tmp_path):
    config = _config(tmp_path, {"cicer": _instance([SOY], profile="all")})
    (instance,) = _plan(checkout, config)
    env = gcv.env_file(instance)
    assert "GCV_PATH=/cicer/gcv/\n" in env
    assert "GCV_CONFIG=./data/cicer/config.json\n" in env
    app = json.loads(gcv.app_config(instance, gcv.load_app_config()))
    urls = {key: app["servers"][0][key]["url"] for key in gcv.APP_SERVICES}
    assert urls["chromosome"] == "/cicer/chromosome/grpc-web"
    assert all(url.startswith("/cicer/") for url in urls.values())


def test_a_services_instance_has_no_web_app(checkout, tmp_path):
    config = _config(tmp_path, {"glycine": _instance([SOY])})
    (instance,) = _plan(checkout, config)
    assert "GCV_PATH" not in gcv.env_file(instance)
    written = gcv.render([instance], config, str(tmp_path / "gcv"))
    assert not any(path.endswith("config.json") for path in written)


def test_the_rendered_tree_matches_the_reference(checkout, tmp_path):
    members = [{"path": COWPEA_ANN1, "rename_chromosomes": True}, COWPEA_ANN2]
    config = _config(
        tmp_path,
        {
            "glycine": _instance([SOY]),
            "lis": _instance(members, profile="all"),
        },
    )
    out = tmp_path / "gcv"
    written = gcv.render(_plan(checkout, config), config, str(out))
    for path in written:
        relative = os.path.relpath(path, out)
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
        if relative.endswith(gcv.INPUTS_SCRIPT):
            with open(gcv.INPUTS_PATH, encoding="utf-8") as handle:
                assert text == handle.read()
            continue
        expected = os.path.join(EXPECTED, relative)
        if UPDATING:
            os.makedirs(os.path.dirname(expected), exist_ok=True)
            with open(expected, "w", encoding="utf-8") as handle:
                handle.write(text)
        with open(expected, encoding="utf-8") as handle:
            assert text == handle.read(), relative
    assert os.access(out / "data/lis/docker-entrypoint-initdb.d/10-load.sh", os.X_OK)


def test_each_instance_is_served_at_its_own_path(checkout, tmp_path):
    config = _config(tmp_path, {"lis": _instance([SOY]), "lis_fam3": _instance([SOY])})
    instances = _plan(checkout, config)
    assert "MICROSERVICES_PATH=/lis_fam3\n" in gcv.env_file(instances[1])


def _catalog(checkout, config):
    index = DatastoreIndex(checkout, datastore_url=URL).build()
    return gcv.catalog_section(index, config, logging.getLogger("test"))


def test_the_catalog_places_each_annotation_an_instance_loads(checkout, tmp_path):
    members = [{"path": COWPEA_ANN1, "rename_chromosomes": True}, COWPEA_ANN2]
    config = _config(tmp_path, {"vigna": _instance(members, profile="all")})
    instances, placements = _catalog(checkout, config)
    assert instances == {
        "vigna": {
            "url": "https://localhost/vigna/gcv/",
            "source": "lis",
            "status": "ok",
            "family": "legfed_v1_0",
            "services": "https://localhost/vigna/",
        }
    }
    placed = [{"instance": "vigna", "source": "lis"}]
    assert placements == {COWPEA_ANN1: placed, COWPEA_ANN2: placed}


def test_an_instance_that_cannot_be_planned_is_unavailable_and_places_nothing(
    checkout, tmp_path
):
    absent = "Glycine/max/annotations/Absent.gnm1.ann1.XXXX"
    config = _config(
        tmp_path,
        {
            "glycine": _instance([SOY, absent], profile="all"),
            "vigna": _instance([COWPEA_ANN2], profile="all"),
        },
    )
    instances, placements = _catalog(checkout, config)
    assert instances["glycine"]["status"] == "unavailable"
    assert "no such collection" in instances["glycine"]["detail"]
    assert instances["vigna"]["status"] == "ok"
    assert list(placements) == [COWPEA_ANN2]


def test_an_instance_without_the_web_app_is_not_catalogued(checkout, tmp_path):
    config = _config(tmp_path, {"glycine": _instance([SOY])})
    assert _catalog(checkout, config) == ({}, {})


def test_the_catalog_records_gcv_only_when_given_a_gcv_yml(checkout, tmp_path):
    _config(tmp_path, {"glycine": _instance([SOY], profile="all")})
    document = CatalogBuilder(checkout, gcv_config=str(tmp_path / "gcv.yml")).build()
    assert document["gcv_instances"]["glycine"]["status"] == "ok"
    soy = next(c for c in document["collections"] if c["path"] == SOY)
    assert soy["gcv"] == [{"instance": "glycine", "source": "lis"}]
    document = CatalogBuilder(checkout).build()
    assert "gcv_instances" not in document
    assert all("gcv" not in c for c in document["collections"])


def test_an_unreadable_gcv_yml_leaves_gcv_out_of_the_catalog(checkout, tmp_path):
    bad = tmp_path / "bad.yml"
    bad.write_text("instances: [")
    document = CatalogBuilder(checkout, gcv_config=str(bad)).build()
    assert "gcv_instances" not in document


def test_chromosomes_span_each_fai_sequence(tmp_path):
    fai = tmp_path / "genome.fna.gz.fai"
    fai.write_text(
        "glyma.Wm82.gnm4.Gm01\t56831624\t52\t60\t61\nscaf;1\t1000\t0\t60\t61\n"
    )
    out = tmp_path / "chromosomes.gff3"
    gcv_inputs.main(["chromosomes", str(fai), str(out)])
    assert out.read_text().splitlines() == [
        "##gff-version 3",
        "glyma.Wm82.gnm4.Gm01\t.\tchromosome\t1\t56831624\t.\t.\t.\t"
        "ID=glyma.Wm82.gnm4.Gm01;Name=glyma.Wm82.gnm4.Gm01",
        "scaf;1\t.\tchromosome\t1\t1000\t.\t.\t.\tID=scaf%3B1;Name=scaf%3B1",
    ]


def test_genes_are_renamed_and_read_up_to_any_fasta(tmp_path):
    gff = tmp_path / "genes.gff3.gz"
    rows = [
        "##gff-version 3",
        "vigun.IT97K-499-35.gnm1.Vu01\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1",
        "vigun.IT97K-499-35.gnm1.Vu01\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1",
        "##FASTA",
        "vigun.IT97K-499-35.gnm1.Vu01\tsrc\tgene\t1\t9\t.\t+\t.\tID=g2",
    ]
    with gzip.open(gff, "wt") as handle:
        handle.write("\n".join(rows) + "\n")
    out = tmp_path / "genes.gff3"
    gcv_inputs.main(["genes", str(gff), str(out), "--rename", "gnm1", "gnm1.ann1"])
    assert out.read_text().splitlines() == [
        "##gff-version 3",
        "vigun.IT97K-499-35.gnm1.ann1.Vu01\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1",
    ]


def test_a_sequence_without_the_renamed_token_fails_without_writing(tmp_path):
    fai = tmp_path / "genome.fna.gz.fai"
    fai.write_text("vigun.IT97K-499-35.gnm1.Vu01\t10\nunplaced\t5\n")
    out = tmp_path / "chromosomes.gff3"
    with pytest.raises(SystemExit, match="unplaced has no .gnm1."):
        gcv_inputs.main(["chromosomes", str(fai), str(out), "--rename", "gnm1", "x"])
    assert not out.exists()
