"""Guard the vendored filetypes.yml against drift from datastore-specifications.

The network check skips when GitHub is unreachable; structural checks always run.
"""

import urllib.request

import pytest
import yaml

from lis_autocontent.datastore_files import FILETYPES_PATH as FILETYPES

SPEC_TREE = (
    "https://api.github.com/repos/legumeinfo/datastore-specifications/"
    "git/trees/main?recursive=1"
)


@pytest.fixture(name="vocabulary")
def fixture_vocabulary():
    with open(FILETYPES, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def test_every_entry_is_well_formed(vocabulary):
    """Behavioural tests only exercise a few types; this is the guard for the rest."""
    for ctype, spec in vocabulary.items():
        assert isinstance(spec, dict), ctype
        if spec.get("predictable") is False:
            assert "contents" not in spec, f"{ctype}: unpredictable but lists contents"
            continue
        assert spec.get("prefix") in ("abbrev", "none"), ctype
        assert spec.get("contents"), ctype
        assert str(spec.get("extension", "")).startswith("."), ctype


def test_empirically_verified_flags_are_pinned(vocabulary):
    """`prefix` and `indexed` aren't in the spec, so the spec check can't guard them.

    Prefixless sets are named like `Cicer.pan2.CMWZ.*`; qtl, gwas and maps publish no
    index siblings.
    """
    prefixless = {t for t, s in vocabulary.items() if s.get("prefix") == "none"}
    unindexed = {t for t, s in vocabulary.items() if s.get("indexed") is False}
    assert prefixless == {"pangenes", "genefamilies"}
    assert unindexed == {"qtl", "gwas", "maps", "mstmap"}


def _spec_contents():
    try:
        with urllib.request.urlopen(SPEC_TREE, timeout=25) as response:
            import json

            tree = json.load(response).get("tree", [])
    except Exception as err:  # noqa: BLE001
        pytest.skip(f"datastore-specifications unreachable: {err}")
    out = {}
    for entry in tree:
        parts = entry.get("path", "").split("/")
        if (
            entry.get("type") == "blob"
            and len(parts) == 4
            and parts[0] == "Genus"
            and parts[3].endswith(".md")
            and parts[3] != "README.md"
        ):
            out.setdefault(parts[2], set()).add(parts[3][:-3])
    return out


@pytest.mark.parametrize("ctype", ["qtl", "gwas", "maps"])
def test_vocabulary_matches_the_specification(vocabulary, ctype):
    """The types whose collections have no CHECKSUM, so prediction is load-bearing."""
    documented = _spec_contents().get(ctype)
    if not documented:
        pytest.skip(f"specification lists no contents for {ctype}")
    assert set(vocabulary[ctype]["contents"]) == documented, (
        f"{ctype}: filetypes.yml has {sorted(vocabulary[ctype]['contents'])}, "
        f"datastore-specifications documents {sorted(documented)}"
    )
