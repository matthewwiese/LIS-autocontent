"""Tests for which taxa a ProcessCollections run covers.

Only taxon selection is tested here: it is offline. Everything past it fetches from
the remote datastore and is exercised by running the CLI.
"""

import logging
import os

import pytest
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


def test_genera_are_discovered_from_their_genus_descriptions(parser):
    """LEGUMES and .github are top-level directories but not genera."""
    assert parser.discover_genera() == ["Cicer", "Glycine"]


def test_no_taxa_list_means_every_genus(parser):
    assert parser.load_taxa() == [{"genus": "Cicer"}, {"genus": "Glycine"}]


def test_a_taxa_list_restricts_the_run(parser, tmp_path):
    taxa = tmp_path / "taxa.yml"
    taxa.write_text(yaml.dump([{"genus": "Cicer", "description": "Chickpea"}]))
    assert parser.load_taxa(str(taxa)) == [
        {"genus": "Cicer", "description": "Chickpea"}
    ]


def test_a_directory_given_as_the_taxa_list_falls_back_to_discovery(parser, clone):
    assert parser.load_taxa(clone) == [{"genus": "Cicer"}, {"genus": "Glycine"}]


def test_discovery_needs_a_clone(parser):
    """Reading from the remote datastore there is no directory to discover from."""
    parser.from_github = None
    with pytest.raises(SystemExit):
        parser.load_taxa()
