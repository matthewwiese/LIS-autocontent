# Working on lis-autocontent

Scrapes the LIS Data Store's metadata and writes the configs and databases that
downstream services deploy from: a Jekyll site, JBrowse2, BLAST, DSCensor, a whole-store
catalog, and a Divbrowse compose file.

## The one rule that matters

**Every subcommand is a pure transform of datastore metadata, and its output is consumed
by something you cannot see.** Treat each artifact as a published contract. A refactor
that "cleans up" the shape of `dscensor_nodes/*.json`, the Jekyll YAML, or the BLAST
command list breaks a consumer, silently, at deploy time.

So: before committing a change to any code that feeds an artifact, prove the artifact
didn't move. Build it before and after against the same `datastore-metadata` checkout and
diff the trees:

```
git clone https://github.com/legumeinfo/datastore-metadata.git   # the input, ~1 min
lis-autocontent populate-dscensor --from_github ./datastore-metadata --nodes_out ./before
git stash && pip install -e . && lis-autocontent populate-dscensor \
    --from_github ./datastore-metadata --nodes_out ./after
diff -r ./before ./after      # must be empty, or every difference must be one you intended
```

An intended difference is fine — an unexplained one is a bug. Say plainly in the commit
message which artifacts changed and why. "No artifact changed" is itself worth stating.

## Layout

| File | Role |
| --- | --- |
| `datastore_files.py` | `DatastoreIndex`: the whole-store model. Reads every collection in a checkout — README metadata, CHECKSUM/MANIFEST file lists, BUSCO, lineage, pairwise relationships — with no network access. **This is the base abstraction; extend it rather than re-reading the store elsewhere.** |
| `catalog.py` | A thin serializer over `DatastoreIndex`, nothing more. Its JSON shape is a contract: change it only together with `SCHEMA_VERSION`. |
| `process_collections.py` | The older scraper path behind `populate-jekyll`, `populate-jbrowse2`, `populate-blast` and `populate-dscensor`. Derived from the original SammyJava scripts; its quirks are load-bearing because its outputs are already deployed. |
| `divbrowse.py` | Builds the Divbrowse compose file. Refuses, with reasons, any collection the format can't express rather than emitting something half-valid. |
| `lis_cli.py` | Click options and wiring only; no logic. |
| `filetypes.yml` | The filetype vocabulary, shipped as package data (see `[tool.setuptools.package-data]`). Not importable Python — it must stay listed there or installs break. |

`CHECKSUM.<key>.md5` is the authoritative file list for a collection, not the naming
conventions and not the store's HTML index: it is the only enumeration that includes the
`.fai`/`.tbi` index siblings, so random-access flags can come from nowhere else.

## Setup and commands

```
pip install -e '.[dev]'     # Python >=3.10,<4
pre-commit install          # end-of-file-fixer, trailing-whitespace, black
pytest                      # 79 tests, ~1s, no network
black src tests             # 88 cols, non-negotiable (pre-commit enforces it)
pylint src/lis_autocontent  # max-line-length 110, fail-under 8.0
```

CI runs pylint and pytest on 3.10 and 3.13. There is no local `actionlint`; lint workflow
files with `docker run --rm -v "$(pwd)":/repo -w /repo rhysd/actionlint:latest -no-color`.

## Tests

Two kinds, and the distinction matters:

- **Unit tests** over fixture trees built by `tests/conftest.py`. Each one pins a specific
  edge case in real datastore metadata — a collection with no CHECKSUM, an unreadable
  README, a pairwise name carrying a program suffix. When you fix a metadata edge case,
  add the fixture that reproduces it.
- **Expected-output ("golden") tests** comparing generated artifacts against reviewed
  copies in `tests/data/`. These are the regression net for the rule above.

When a golden test fails, the failure is the point. Read the diff first and decide whether
the change was intended. Only then:

```
LIS_UPDATE_EXPECTED=1 pytest    # regenerates tests/data/, then review that diff and commit it
```

Never regenerate to make a red test green. A reviewed `tests/data/` diff in the commit is
the evidence that an artifact change was deliberate.

Avoid tests written for coverage's own sake. A test earns its place by pinning a metadata
edge case or an artifact's shape.

## Gotchas found the hard way

- **HEAD requests have no body.** `head_remote` must return a boolean from the status code;
  returning `response.text` makes every probe falsy, which once silently dropped 263
  protein collections from DSCensor and BLAST for two years. Any "file absent" path
  deserves this suspicion.
- **`--verify` is off by default and should stay off in automation.** It costs ~1,650 HEAD
  requests, and a probe that fails for *any* reason (a blip, a rate limit) reads as "file
  absent" and drops the file — so a verified build is less reproducible than an offline one.
- **`index_status` is four-state** (known/verified/inferred/unknown). ~45% of collections
  are `unknown`; collapsing that to "no indexes" reports streamable data as unreadable.
- Docker Compose service names must be lowercase, or the image tag is rejected.

## Related repositories

- **`legumeinfo/datastore-metadata`** — the input. Every subcommand reads a checkout of it
  via `--from_github` (default `./datastore-metadata`).
- The **catalog release workflow lives in datastore-metadata**, not here: it checks this
  repo out, runs `populate-catalog`, and publishes `catalog.json` as a release asset there.
  The catalog is built *from* that repo's contents, so building it there avoids a
  cross-repository token and a second run to chase.

## Conventions

Commit subjects are imperative and lowercase-ish after the first word ("Add populate-catalog
action", "Publish the catalog as a release instead of a workflow artifact"). Explain *why*
in the body, and name any artifact that changed. Comments in this codebase explain
reasoning, not mechanics — match that; don't narrate what the next line plainly does.
