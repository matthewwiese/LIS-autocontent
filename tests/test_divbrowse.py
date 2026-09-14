"""Tests for populate-divbrowse: a Divbrowse docker-compose.yml built from metadata.

The reference is legumeinfo/divbrowse's own docker-compose.yml, vendored byte for byte
in tests/data/divbrowse/ (last changed upstream in eaf28a6). For the three collections
it serves, the generated file must be that file with only the service names and data
directories changed to the collection identifiers. A collection the format can't express
must stop the command with its reason and write nothing.
"""

import os
import re

import pytest
import yaml
from click.testing import CliRunner

from lis_autocontent import lis_cli

UPSTREAM = os.path.join(
    os.path.dirname(__file__), "data", "divbrowse", "upstream-docker-compose.yml"
)
MD5 = "d41d8cd98f00b204e9800998ecf8427e"
SOYBEAN = (
    "Wm82.gnm4.div.Song_Hyten_2015",
    "Wm82.gnm5.div.Song_Hyten_2015",
    "Wm82.gnm6.div.Song_Hyten_2015",
)
# Upstream's hand-picked names, and what they become. Nothing else may differ.
RENAMES = {
    "  divbrowse-gnm4:": "  divbrowse-wm82.gnm4.div.song_hyten_2015:",
    "  divbrowse-gnm5:": "  divbrowse-wm82.gnm5.div.song_hyten_2015:",
    "  divbrowse-gnm6:": "  divbrowse-wm82.gnm6.div.song_hyten_2015:",
    "./data/gnm4:": "./data/Wm82.gnm4.div.Song_Hyten_2015:",
    "./data/gnm5:": "./data/Wm82.gnm5.div.Song_Hyten_2015:",
    "./data/gnm6:": "./data/Wm82.gnm6.div.Song_Hyten_2015:",
}


def _checksum(*names):
    return "\n".join(f"{MD5}  ./{name}" for name in names) + "\n"


def _assembly(write, root, strain_gnm, key, prefix, annotations=("ann1.AAAA",)):
    """A Glycine max genome with its annotations, as datastore-metadata lays them out."""
    genome = f"{strain_gnm}.{key}"
    write(
        os.path.join(root, f"Glycine/max/genomes/{genome}/README.{genome}.yml"),
        yaml.safe_dump(
            {
                "identifier": genome,
                "scientific_name_abbrev": "glyma",
                "chromosome_prefix": prefix,
            }
        ),
    )
    for annotation in annotations:
        ann = f"{strain_gnm}.{annotation}"
        write(
            os.path.join(root, f"Glycine/max/annotations/{ann}/README.{ann}.yml"),
            yaml.safe_dump({"identifier": ann, "scientific_name_abbrev": "glyma"}),
        )
        gff3 = f"glyma.{ann}.gene_models_main.gff3.gz"
        write(
            os.path.join(root, f"Glycine/max/annotations/{ann}/CHECKSUM.{ann}.md5"),
            _checksum(gff3, gff3 + ".tbi", f"glyma.{ann}.protein.faa.gz"),
        )


def _diversity(write, root, identifier, files):
    write(
        os.path.join(
            root, f"Glycine/max/diversity/{identifier}/README.{identifier}.yml"
        ),
        yaml.safe_dump({"identifier": identifier, "scientific_name_abbrev": "glyma"}),
    )
    write(
        os.path.join(
            root, f"Glycine/max/diversity/{identifier}/CHECKSUM.{identifier}.md5"
        ),
        _checksum(*files),
    )


@pytest.fixture(name="clone")
def fixture_clone(tmp_path, write):
    """The metadata behind upstream's three services, mirroring the real store:
    gnm5 names its chromosomes Chr, gnm4 and gnm6 name them Gm."""
    root = str(tmp_path / "datastore-metadata")
    for strain_gnm, key, prefix, annotation in (
        ("Wm82.gnm4", "4PTR", "Gm", "ann1.T8TQ"),
        ("Wm82.gnm5", "NRKG", "Chr", "ann1.J7HW"),
        ("Wm82.gnm6", "S97D", "Gm", "ann1.PKSW"),
    ):
        _assembly(write, root, strain_gnm, key, prefix, annotations=(annotation,))
        identifier = f"{strain_gnm}.div.Song_Hyten_2015"
        vcf = f"glyma.{identifier}.vcf.gz"
        _diversity(write, root, identifier, (vcf, vcf + ".tbi"))
    return root


def _run(clone, tmp_path, collections, *options):
    out = tmp_path / "divbrowse" / "docker-compose.yml"
    args = ["--from_github", clone, "--compose_out", str(out)]
    args += ["--log_file", str(tmp_path / "divbrowse.log"), *options]
    for identifier in collections:
        args += ["--collection", identifier]
    result = CliRunner().invoke(lis_cli.populate_divbrowse, args)
    return result, out


def test_the_upstream_file_is_reproduced_with_collection_names(clone, tmp_path):
    result, out = _run(clone, tmp_path, SOYBEAN)
    assert result.exit_code == 0, result.output
    with open(UPSTREAM, encoding="utf-8") as handle:
        expected = handle.read()
    for old, new in RENAMES.items():
        assert expected.count(old) == 1, old  # each rename lands exactly once
        expected = expected.replace(old, new)
    assert out.read_text(encoding="utf-8") == expected


def test_service_names_are_valid_image_names(clone, tmp_path):
    """Compose names each service's image after it, and `docker compose build` rejects
    uppercase; every diversity identifier in the store has some."""
    _, out = _run(clone, tmp_path, SOYBEAN)
    services = yaml.safe_load(out.read_text(encoding="utf-8"))["services"]
    component = re.compile(r"[a-z0-9]+(?:(?:[._]|__|[-]+)[a-z0-9]+)*")
    assert all(component.fullmatch(name) for name in services)
    assert [s["volumes"] for s in services.values()] == [
        [f"./data/{identifier}:/opt/divbrowse"] for identifier in SOYBEAN
    ]


def test_deployment_options_reach_the_file(clone, tmp_path):
    options = (
        "--base_url",
        "https://example.org/divbrowse/",
        "--port",
        "9000",
        "--datastore_url",
        "https://mirror.example.org",
    )
    result, out = _run(clone, tmp_path, SOYBEAN[1:], *options)
    assert result.exit_code == 0, result.output
    services = list(
        yaml.safe_load(out.read_text(encoding="utf-8"))["services"].values()
    )
    assert [s["ports"] for s in services] == [["9000:8080"], ["9001:8080"]]
    environment = dict(entry.split("=", 1) for entry in services[0]["environment"])
    assert (
        environment["BASE_URL"]
        == "https://example.org/divbrowse/Wm82.gnm5.div.Song_Hyten_2015/"
    )
    assert environment["VCF_URL"].startswith(
        "https://mirror.example.org/Glycine/max/diversity/"
    )
    assert environment["GFF3_URL"].startswith(
        "https://mirror.example.org/Glycine/max/annotations/"
    )


def _several_vcfs(write, root):
    _diversity(write, root, "Wm82.gnm4.div.Two_VCF_2020", ("a.vcf.gz", "b.vcf.gz"))
    return "Wm82.gnm4.div.Two_VCF_2020", "has 2 VCF files: a.vcf.gz, b.vcf.gz"


def _no_vcf(write, root):
    _diversity(write, root, "Wm82.gnm4.div.No_VCF_2020", ("glyma.SNPs.txt.gz",))
    return "Wm82.gnm4.div.No_VCF_2020", "has 0 VCF files"


def _unlinked(write, root):
    _diversity(write, root, "CDCFrontier.div.Unlinked_2020", ("x.vcf.gz",))
    return "CDCFrontier.div.Unlinked_2020", "not linked to a genome assembly"


def _two_annotations(write, root):
    _assembly(
        write, root, "Lee.gnm1", "ABCD", "Gm", annotations=("ann1.AAAA", "ann2.BBBB")
    )
    _diversity(write, root, "Lee.gnm1.div.Two_Ann_2020", ("x.vcf.gz",))
    return (
        "Lee.gnm1.div.Two_Ann_2020",
        "has 2 annotations: Lee.gnm1.ann1.AAAA, Lee.gnm1.ann2.BBBB",
    )


def _chromosome_list(write, root):
    """Weilv-9.gnm1's README gives chromosome names where a prefix belongs."""
    _assembly(write, root, "Zh13.gnm1", "CCCC", "chr,Pt,Mt")
    _diversity(write, root, "Zh13.gnm1.div.Prefix_List_2020", ("x.vcf.gz",))
    return (
        "Zh13.gnm1.div.Prefix_List_2020",
        "chromosome_prefix 'chr,Pt,Mt', which is not a single prefix",
    )


def _unknown(write, root):  # pylint: disable=unused-argument
    return "Wm82.gnm9.div.Nowhere_2099", "no diversity collection has this identifier"


@pytest.mark.parametrize(
    "case",
    [_several_vcfs, _no_vcf, _unlinked, _two_annotations, _chromosome_list, _unknown],
)
def test_collections_the_format_cannot_express_stop_the_build(
    clone, tmp_path, write, case
):
    """Alongside valid collections, one unusable collection still writes nothing:
    a compose file with a guessed VCF, annotation or pattern is worse than none."""
    identifier, reason = case(write, clone)
    result, out = _run(clone, tmp_path, (*SOYBEAN, identifier))
    assert result.exit_code != 0
    assert f"{identifier}: " in result.output and reason in result.output
    assert not out.exists()


def test_every_problem_is_reported_at_once(clone, tmp_path, write):
    several, several_reason = _several_vcfs(write, clone)
    unknown, unknown_reason = _unknown(write, clone)
    result, out = _run(clone, tmp_path, (SOYBEAN[0], several, unknown, SOYBEAN[0]))
    assert result.exit_code != 0
    for expected in (several_reason, unknown_reason, "requested more than once"):
        assert expected in result.output
    assert not out.exists()
