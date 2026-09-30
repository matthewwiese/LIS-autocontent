# Working on lis-autocontent

Reads a checkout of `legumeinfo/datastore-metadata` and writes what downstream LIS services
deploy from: Jekyll YAML, JBrowse2 and BLAST commands, DSCensor nodes, a whole-store
catalog, and a Divbrowse compose file.

## Artifacts are contracts

Each subcommand is a pure transform of metadata, consumed by a service you can't see. An
unintended change to an artifact's content or shape is a bug, even when tests pass.

Before committing a change to code that feeds an artifact, build it before and after
against one checkout and diff:

```
git clone https://github.com/legumeinfo/datastore-metadata.git
lis-autocontent populate-dscensor --from_github ./datastore-metadata --nodes_out ./before
# apply the change, reinstall, then:
lis-autocontent populate-dscensor --from_github ./datastore-metadata --nodes_out ./after
diff -r ./before ./after
```

Every difference must be intended and named in the commit message.

## Architecture

### Input

datastore-metadata mirrors the Data Store's tree, holding only metadata:

```
<Genus>/<species>/<collection_type>/<collection>/
    README.<collection>.yml       descriptive metadata (a few ship a bare README)
    CHECKSUM.<collection>.md5     every file in the collection, with md5
    MANIFEST.<collection>.yml     per-file descriptions and applications
    BUSCO/*.short_summary.json    completeness and assembly metrics
<Genus>/{GENUS,<species>}/about_this_collection/description_*.yml   taxon metadata
<Genus>/<species>/gene_functions/<abbrev>.traits.yml                curated gene symbols
```

### Two pipelines

| | `DatastoreIndex` (`datastore_files.py`) | `ProcessCollections` (`process_collections.py`) |
| --- | --- | --- |
| Serves | `populate-catalog`, `populate-divbrowse` | `populate-jekyll`, `-jbrowse2`, `-blast`, `-dscensor` |
| Reads | the checkout only; offline | the checkout plus HEAD probes of the live store; the store alone without a checkout |
| Model | typed dataclasses | nested dicts keyed by collection type |
| Role | the base for new work | legacy; outputs already deployed, change conservatively |

They share one constant: `NODE_README_FIELDS` (DSCensor's README fields) must stay a subset
of `README_FIELDS` (the catalog's). A test enforces it.

### `DatastoreIndex.build()`

1. **Discover.** `find_collections` walks the checkout. A directory exactly four levels deep
   (`Genus/species/type/collection`) is a collection if it holds a README, CHECKSUM or
   MANIFEST.
2. **Read the README.** Unreadable yields `{}` and a warning, never a crash. Non-empty
   `README_FIELDS` go to `collection.metadata`.
3. **Resolve files**, most authoritative first:
   - `checksum_files`: every `CHECKSUM.*.md5` in the directory. `src=checksum`,
     `index_status=known`, md5 carried.
   - else `predicted_records`: names built from `filetypes.yml` conventions.
     `src=predicted`, `index_status=inferred`; with `--verify`, HEAD-confirmed as
     `verified`. A type with no convention stays `unknown`, with no files.
4. **Annotate files.** `indexes` from sibling suffixes present in the list
   (`INDEX_SUFFIXES`); `description`/`applications` from MANIFEST; `busco`/`counts` from the
   BUSCO JSON.
5. **Link lineage.** `link_lineage` derives `derived_from` from identifiers
   (`Wm82.gnm4.ann1.T8TQ` -> `Wm82.gnm4`) and copies `INHERITED_FIELDS` such as
   `chromosome_prefix` down from the parent genome.
6. **Store-wide passes.** `pairwise_relationships` (genome pairs parsed from filenames by
   `PAIRWISE_PATTERN`), `descriptions` (taxa), `curated_symbols`, then `count` for stats.

Query the built index with `select`, `sibling`, `companion` and `orphan_collections`.

### Provenance

Nothing is presented as more certain than it is:

| Field | Values |
| --- | --- |
| `DatastoreFiles.src` | `checksum` · `verified` · `predicted` |
| `DatastoreCollection.index_status` | `known` · `verified` · `inferred` · `unknown` |
| `DatastoreFiles.index_unknown` | true when a predicted file's index siblings were never probed |

`unknown` is not "none": collapsing it reports streamable data as unreadable. Likewise
`checksum_files` returns `None` for "could not look" and `[]` for "looked, empty", and
callers depend on the difference.

### Consumers of the index

- `catalog.py` serializes the index to one JSON document and holds no logic of its own.
  Its shape is versioned: change it only with `SCHEMA_VERSION`.
- `divbrowse.py` builds one compose service per diversity collection and refuses any
  collection the format can't express. Per-genus hosts and combined references are
  data in `divbrowse.yml`, package data like `filetypes.yml`.

New metadata belongs on `DatastoreIndex` (a field or a pass) and reaches consumers through
the catalog. Don't re-walk the checkout elsewhere. Change `ProcessCollections` only to keep
its existing outputs correct.

## Code style

black (88 cols), pylint (`max-line-length` 110, `fail-under` 8) and pre-commit enforce the
basics. Beyond those:

- **Records are `@dataclass`es.** Derived values are `@property`; pure helpers are
  `@staticmethod`s or module functions.
- **Vocabularies are module constants** (`SRC_*`, `STATUS_*`, `README_FIELDS`), never
  repeated string literals. Regexes are compiled module constants; name the groups
  you capture.
- **Offline and deterministic.** Sort collections and files. Apart from the catalog's
  `built_at` stamp, an artifact depends only on the checkout: no set iteration order, and
  no network access except behind an explicit flag (`--verify`).
- **Degrade loudly.** Missing or unreadable input yields an empty value and a
  `logger.warning`: never a crash, never a silent default.
- **Validate completely.** Collect every problem, then fail once listing them all
  (`DivbrowseError`). Validators return `(value, problems)`.
- **Lazy logging:** `logger.info("wrote %s", path)`, not f-strings.
- **Narrow pylint disables,** inline on the line that needs one, never file-wide.
- **`lis_cli.py` is wiring:** click options and calls, no logic.
- **Tests:** declare fixtures as `@pytest.fixture(name="x")` over `def fixture_x()`. Name
  tests as sentences stating the behavior
  (`test_orphan_collections_are_those_without_a_checksum`). Build fixture trees with
  conftest's `write`.

`process_collections.py` predates these rules. Conform the lines you touch; leave the rest.

## Comments

- **Terse.** One line, two at most. A docstring is a one-line summary, plus a short
  paragraph only for a contract a caller must know.
- **Why, not what.** Don't restate the code.
- **The code as it is, in the present tense.** No history: not "fixed", "used to",
  "previously", "now", nor what a bug did or how it was found. That goes in the commit
  message.

Bad:

```python
def head_remote(self, url):
    """Returns True if it exists, otherwise False.

    A HEAD response has no body, so returning response.text here gave "" for every
    file that exists, and callers treated each one as missing.
    """
```

Good:

```python
def head_remote(self, url):
    """True if the URL answers 200. HEAD has no body, so decide on status alone."""
```

Bad:

```python
# Types documented as never carrying index siblings need no probing to know
# none of these files is randomly accessible. For the types that CAN carry
# them, an unprobed file's index status is genuinely unknown -- reporting it
# as "not indexed" would mark streamable data unreadable, which is the same
# mistake `index_status` exists to prevent.
```

Good:

```python
# Unprobed siblings are unknown, not absent.
```

## Tests

- **Unit tests** run over fixture trees from `tests/conftest.py`, each pinning one
  metadata edge case. Fixing an edge case means adding the fixture that reproduces it.
- **Golden tests** compare artifacts with reviewed copies in `tests/data/`.

A failing golden test is the point. Read its diff and decide whether the change is
intended; only then run `LIS_UPDATE_EXPECTED=1 pytest` and commit the reviewed
`tests/data/` diff. Never regenerate to turn a test green. Don't write tests for coverage.

## Setup and commands

```
pip install -e '.[dev]'     # Python >=3.10,<4
pre-commit install
pytest                      # 79 tests, ~1s, no network
black src tests
pylint src/lis_autocontent
```

CI runs pylint and pytest on Python 3.10 and 3.13. `populate-jbrowse2` and
`populate-blast` execute `jbrowse`/`makeblastdb` unless given `--cmds_only`.

## Invariants

- Existence checks decide on status code; a HEAD response has no body.
- `--verify` stays off in automation: a probe failing for any reason reads as "absent" and
  drops the file.
- `CHECKSUM` is the only enumeration listing `.fai`/`.tbi` siblings, so random-access flags
  come from nowhere else.
- Compose service names are lowercase, or Docker rejects the image tag.

## Related repositories

- `legumeinfo/datastore-metadata` is the input (`--from_github`, default
  `./datastore-metadata`).
- The catalog release workflow lives in datastore-metadata: it checks this repo out, runs
  `populate-catalog`, and publishes `catalog.json` as a release asset there.

## Commits

Imperative subject. The body says why, and names every artifact that changed, or states
that none did. History belongs here, not in comments.
