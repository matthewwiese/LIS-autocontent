"""populate-divbrowse against the Traefik layout drafted for legumeinfo/divbrowse.

The reference in tests/data/divbrowse/ is that draft. The output must equal it except
for service names and data directories.
"""

import os
import re

import pytest
import yaml
from click.testing import CliRunner

from lis_autocontent import lis_cli
from lis_autocontent.datastore_files import DatastoreIndex
from lis_autocontent.divbrowse import DivbrowseError, compose_file, load_config

REFERENCE = os.path.join(
    os.path.dirname(__file__), "data", "divbrowse", "traefik-docker-compose.yml"
)
MD5 = "d41d8cd98f00b204e9800998ecf8427e"
SOYBEAN = (
    "Wm82.gnm4.div.Song_Hyten_2015",
    "Wm82.gnm5.div.Song_Hyten_2015",
    "Wm82.gnm6.div.Song_Hyten_2015",
)
# Called against aradu1_araip1.gnm1, the two progenitor genomes concatenated.
PEANUT = (
    "aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020",
    "aradu1_araip1.gnm1.div.Otyama_Wilkey_2019",
    "aradu1_araip1.gnm1.div.Clevenger_Korani_2018",
)
# The reference's hand-picked names, and what they become. Nothing else may differ.
RENAMES = {
    "  divbrowse-gnm4:": "  divbrowse-wm82.gnm4.div.song_hyten_2015:",
    "  divbrowse-gnm5:": "  divbrowse-wm82.gnm5.div.song_hyten_2015:",
    "  divbrowse-gnm6:": "  divbrowse-wm82.gnm6.div.song_hyten_2015:",
    "./data/gnm4:": "./data/Wm82.gnm4.div.Song_Hyten_2015:",
    "./data/gnm5:": "./data/Wm82.gnm5.div.Song_Hyten_2015:",
    "./data/gnm6:": "./data/Wm82.gnm6.div.Song_Hyten_2015:",
    "  divbrowse-arahy-otyama-kulkarni-2020:": (
        "  divbrowse-aradu1_araip1.gnm1.div.otyama_kulkarni_2020:"
    ),
    "  divbrowse-arahy-otyama-wilkey-2019:": (
        "  divbrowse-aradu1_araip1.gnm1.div.otyama_wilkey_2019:"
    ),
    "  divbrowse-arahy-clevenger-korani-2018:": (
        "  divbrowse-aradu1_araip1.gnm1.div.clevenger_korani_2018:"
    ),
    "./data/arahy-otyama-kulkarni-2020:": (
        "./data/aradu1_araip1.gnm1.div.Otyama_Kulkarni_2020:"
    ),
    "./data/arahy-otyama-wilkey-2019:": (
        "./data/aradu1_araip1.gnm1.div.Otyama_Wilkey_2019:"
    ),
    "./data/arahy-clevenger-korani-2018:": (
        "./data/aradu1_araip1.gnm1.div.Clevenger_Korani_2018:"
    ),
}


def _checksum(*names):
    return "\n".join(f"{MD5}  ./{name}" for name in names) + "\n"


def _assembly(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    write,
    root,
    strain_gnm,
    key,
    prefix,
    annotations=("ann1.AAAA",),
    taxon=("Glycine/max", "glyma"),
):
    """A genome with its annotations, as datastore-metadata lays them out."""
    species, abbrev = taxon
    genome = f"{strain_gnm}.{key}"
    write(
        os.path.join(root, f"{species}/genomes/{genome}/README.{genome}.yml"),
        yaml.safe_dump(
            {
                "identifier": genome,
                "scientific_name_abbrev": abbrev,
                "chromosome_prefix": prefix,
            }
        ),
    )
    for annotation in annotations:
        ann = f"{strain_gnm}.{annotation}"
        write(
            os.path.join(root, f"{species}/annotations/{ann}/README.{ann}.yml"),
            yaml.safe_dump({"identifier": ann, "scientific_name_abbrev": abbrev}),
        )
        gff3 = f"{abbrev}.{ann}.gene_models_main.gff3.gz"
        write(
            os.path.join(root, f"{species}/annotations/{ann}/CHECKSUM.{ann}.md5"),
            _checksum(gff3, gff3 + ".tbi", f"{abbrev}.{ann}.protein.faa.gz"),
        )


def _diversity(write, root, identifier, files, taxon=("Glycine/max", "glyma")):
    species, abbrev = taxon
    base = os.path.join(root, f"{species}/diversity/{identifier}")
    write(
        os.path.join(base, f"README.{identifier}.yml"),
        yaml.safe_dump({"identifier": identifier, "scientific_name_abbrev": abbrev}),
    )
    write(os.path.join(base, f"CHECKSUM.{identifier}.md5"), _checksum(*files))


def _peanut_diversity(write, root, identifier, vcfs):
    files = [name for vcf in vcfs for name in (vcf, vcf + ".tbi")]
    _diversity(write, root, identifier, files, taxon=("Arachis/hypogaea", "arahy"))


@pytest.fixture(name="clone")
def fixture_clone(tmp_path, write):
    """The metadata behind the reference's services, mirroring the real store: gnm5
    names its chromosomes Chr, gnm4 and gnm6 Gm; peanut's progenitors Aradu.A and
    Araip.B, and Kulkarni publishes a subset VCF beside its main one."""
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
    _assembly(
        write,
        root,
        "V14167.gnm1",
        "SWBf",
        "Aradu.A",
        annotations=("ann1.cxSM",),
        taxon=("Arachis/duranensis", "aradu"),
    )
    _assembly(
        write,
        root,
        "K30076.gnm1",
        "bXJ8",
        "Araip.B",
        annotations=("ann1.J37m",),
        taxon=("Arachis/ipaensis", "araip"),
    )
    kulkarni, wilkey, clevenger = PEANUT
    _peanut_diversity(
        write, root, kulkarni, (f"{kulkarni}.main.vcf.gz", f"{kulkarni}.sub10k.vcf.gz")
    )
    _peanut_diversity(write, root, wilkey, (f"arahy.{wilkey}.snp_chip.vcf.gz",))
    _peanut_diversity(write, root, clevenger, (f"{clevenger}.snp_chip.vcf.gz",))
    return root


def _run(clone, tmp_path, collections, *options):
    out = tmp_path / "divbrowse" / "docker-compose.yml"
    args = ["--from_github", clone, "--compose_out", str(out)]
    args += ["--log_file", str(tmp_path / "divbrowse.log"), *options]
    for identifier in collections:
        args += ["--collection", identifier]
    result = CliRunner().invoke(lis_cli.populate_divbrowse, args)
    return result, out


def test_the_traefik_layout_is_reproduced_with_collection_names(clone, tmp_path):
    result, out = _run(clone, tmp_path, SOYBEAN + PEANUT)
    assert result.exit_code == 0, result.output
    with open(REFERENCE, encoding="utf-8") as handle:
        expected = handle.read()
    for old, new in RENAMES.items():
        assert expected.count(old) == 1, old  # each rename lands exactly once
        expected = expected.replace(old, new)
    assert out.read_text(encoding="utf-8") == expected


def test_service_names_are_valid_image_names(clone, tmp_path):
    """Compose names each service's image after it, and `docker compose build` rejects
    uppercase; every diversity identifier in the store has some."""
    _, out = _run(clone, tmp_path, SOYBEAN + PEANUT)
    services = yaml.safe_load(out.read_text(encoding="utf-8"))["services"]
    services.pop("proxy")
    component = re.compile(r"[a-z0-9]+(?:(?:[._]|__|[-]+)[a-z0-9]+)*")
    assert all(component.fullmatch(name) for name in services)
    assert [s["volumes"] for s in services.values()] == [
        [f"./data/{identifier}:/opt/divbrowse"] for identifier in SOYBEAN + PEANUT
    ]


def test_deployment_options_reach_the_file(clone, tmp_path):
    options = (
        "--host",
        "Glycine=divbrowse.example.org",
        "--datastore_url",
        "https://mirror.example.org",
    )
    result, out = _run(clone, tmp_path, SOYBEAN[1:], *options)
    assert result.exit_code == 0, result.output
    services = yaml.safe_load(out.read_text(encoding="utf-8"))["services"]
    services.pop("proxy")
    service = services["divbrowse-wm82.gnm5.div.song_hyten_2015"]
    assert service["labels"] == {
        "divbrowse.host": "divbrowse.example.org",
        "divbrowse.path": "Wm82.gnm5.div.Song_Hyten_2015",
    }
    assert "ports" not in service  # the proxy is the only way in
    environment = dict(entry.split("=", 1) for entry in service["environment"])
    assert "BASE_URL" not in environment
    assert environment["VCF_URL"].startswith(
        "https://mirror.example.org/Glycine/max/diversity/"
    )
    assert environment["GFF3_URL"].startswith(
        "https://mirror.example.org/Glycine/max/annotations/"
    )


def test_a_genus_without_a_host_stops_the_build(clone):
    index = DatastoreIndex(clone).build()
    with pytest.raises(DivbrowseError) as err:
        compose_file(index, SOYBEAN[:1], hosts={"Arachis": "divbrowse.peanutbase.org"})
    assert "no Divbrowse host for genus Glycine" in str(err.value)


def test_a_combined_reference_needs_every_genome(tmp_path, write):
    root = str(tmp_path / "datastore-metadata")
    identifier = PEANUT[1]
    _peanut_diversity(write, root, identifier, (f"arahy.{identifier}.vcf.gz",))
    with pytest.raises(DivbrowseError) as err:
        compose_file(DatastoreIndex(root).build(), [identifier])
    assert (
        "its combined reference lacks Arachis/duranensis/genomes/V14167.gnm1.SWBf, "
        "Arachis/ipaensis/genomes/K30076.gnm1.bXJ8" in str(err.value)
    )


@pytest.mark.parametrize(
    "text, reason",
    [
        (None, "cannot read"),
        ("hosts: [divbrowse.soybase.org]\n", "hosts must map each genus"),
        (
            "hosts: {}\ncombined_references: {x.gnm1: {note: [a]}}\n",
            "combined reference x.gnm1 needs a genomes list",
        ),
    ],
)
def test_an_unusable_config_is_reported(tmp_path, text, reason):
    path = tmp_path / "divbrowse.yml"
    if text is not None:
        path.write_text(text, encoding="utf-8")
    with pytest.raises(DivbrowseError) as err:
        load_config(str(path))
    assert reason in str(err.value)


@pytest.mark.parametrize("host", ["divbrowse.example.org", "Glycine=bad host"])
def test_malformed_hosts_write_nothing(clone, tmp_path, host):
    result, out = _run(clone, tmp_path, SOYBEAN, "--host", host)
    assert result.exit_code != 0
    assert not out.exists()


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
