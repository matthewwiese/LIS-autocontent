"""Guard the vendored file-content vocabulary against drift.

`src/lis_autocontent/filetypes.yml` is a copy of what datastore-specifications documents in its
template tree (Genus/species/<type>/<content>.md). A copy that silently falls behind is
worse than no copy: predictions would quietly stop covering a content type, and the
files would disappear from the catalog with no signal.

The network check is skipped when GitHub is unreachable, so the suite stays runnable
offline; the structural checks always run.
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
    """`prefix` and `indexed` are not in datastore-specifications, so the spec check
    below cannot catch them drifting. Both were verified against the live store:
    prefixless sets are named `Cicer.pan2.CMWZ.*` / `legume.fam1.M65K.*` (an abbrev in
    front yields filenames that 404), and 72 probes across qtl/gwas/maps found no
    index siblings (marking an indexable type unindexed hides streamable data)."""
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
