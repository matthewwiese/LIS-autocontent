"""Tests for the offline parts of ProcessCollections.

Taxon selection needs no network. Remote existence checks are tested with requests
stubbed out. Everything else fetches from the remote datastore.
"""

import logging
import os
from types import SimpleNamespace

import pytest
import requests
import yaml

from lis_autocontent.process_collections import ProcessCollections


@pytest.fixture(name="clone")
def fixture_clone(tmp_path, write):
    """A datastore-metadata clone with two genera and two non-genus directories."""
    for genus in ("Glycine", "Cicer"):
        write(
            os.path.join(
                tmp_path,
                genus,
                "GENUS",
                "about_this_collection",
                f"description_{genus}.yml",
            ),
            f"genus: {genus}\n",
        )
    write(os.path.join(tmp_path, "LEGUMES", "Fabaceae", "README.md"), "")
    write(os.path.join(tmp_path, ".github", "workflows", "ci.yml"), "")
    return str(tmp_path)


@pytest.fixture(name="parser")
def fixture_parser(clone):
    parser = ProcessCollections(logging.getLogger("test"))
    parser.from_github = clone
    return parser


def test_no_taxa_list_means_every_genus(parser):
    """Genera are discovered from their GENUS descriptions: LEGUMES and .github are
    top-level directories but not genera."""
    assert parser.load_taxa() == [{"genus": "Cicer"}, {"genus": "Glycine"}]


def test_parse_collections_processes_the_taxa_list_it_is_given(
    parser, clone, tmp_path, monkeypatch
):
    taxa = tmp_path / "taxa.yml"
    taxa.write_text(yaml.dump([{"genus": "Cicer", "description": "Chickpea"}]))
    processed = []
    monkeypatch.setattr(parser, "process_taxon", processed.append)
    parser.parse_collections(clone, str(taxa))
    assert processed == [{"genus": "Cicer", "description": "Chickpea"}]


def test_a_directory_given_as_the_taxa_list_falls_back_to_discovery(parser, clone):
    assert parser.load_taxa(clone) == [{"genus": "Cicer"}, {"genus": "Glycine"}]


def test_discovery_needs_a_clone(parser):
    """Reading from the remote datastore there is no directory to discover from."""
    parser.from_github = None
    with pytest.raises(SystemExit):
        parser.load_taxa()


def _stub_head(status, asked):
    """A requests.head replacement. Like the real thing, its response has no body."""

    def fake_head(url, timeout=None):  # pylint: disable=unused-argument
        asked.append(url)
        return SimpleNamespace(status_code=status, text="")

    return fake_head


@pytest.mark.parametrize(
    "status, proteins", [(200, [".protein", ".protein_primary"]), (404, [])]
)
def test_annotations_yield_the_protein_files_that_exist(
    parser, clone, write, monkeypatch, status, proteins
):
    """The regression: a HEAD response has no body, so returning its text made every
    existence check fail, and no protein or protein_primary file ever became a
    DSCensor node or a BLAST database. A file that is absent must still be skipped."""
    key = "CDCFrontier.gnm3.ann1.NPD7"
    write(
        os.path.join(
            clone, "Cicer", "arietinum", "annotations", key, f"README.{key}.yml"
        ),
        f"identifier: {key}\nsynopsis: Chickpea annotation\n",
    )
    asked = []
    monkeypatch.setattr(requests, "head", _stub_head(status, asked))
    parser.add_collections("annotations", "Cicer", "arietinum")
    base = f"https://data.legumeinfo.org/Cicer/arietinum/annotations/{key}/cicar.{key}"
    assert asked == [f"{base}.protein_primary.faa.gz", f"{base}.protein.faa.gz"]
    lookup = "cicar.CDCFrontier.gnm3.ann1"
    assert sorted(parser.files["annotations"]) == sorted(
        [lookup] + [lookup + suffix for suffix in proteins]
    )
