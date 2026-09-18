# lis-autocontent
Populates various configs and databases for deployment from the datastore-meta data github.

## Reqiurements

1. JBrowse2 (https://jbrowse.org/jb2/docs/combined/)
2. NCBI-BLAST+ (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
3. Python3.10+

## Quick Install With pip

1. Create a virtual environment for python3. `python3 -m venv lis_autocontent_env` (optional)
2. Source environment. `. ./lis_autocontent_env/bin/activate` (optional)
3. pip install from the repo. `cd LIS-autocontent;pip install .`

## Developer Install

1. Clone repository. `git clone https://github.com/legumeinfo/LIS-autocontent.git`
2. Create a virtual environment for python3 `python3 -m venv lis_autocontent_env`  (optional)
3. Source virtual environment `. ./lis_autocontent_env/bin/activate`  (optional)
4. Install Black and pre-commit for git hooks. `pip install black pre-commit`
5. Initialize pre-commit (if you haven't already). `pip install pre-commit;pre-commit install`
6. Install the package in editable mode with its dev dependencies. `pip install -e '.[dev]'`
7. Run the tests. `pytest`

   Some tests compare generated artifacts (the catalog, DSCensor nodes, BLAST and
   JBrowse2 commands, Jekyll YAML) with reviewed copies in `tests/data/`. If one fails
   because the change to that artifact is intended, regenerate the copies with
   `LIS_UPDATE_EXPECTED=1 pytest` and review the diff before committing it.

## Run

```
(lis_autocontent_env) $ lis-autocontent --help
Usage: lis-autocontent [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  populate-blast     CLI entry for populate-blast
  populate-catalog   CLI entry for populate-catalog
  populate-divbrowse CLI entry for populate-divbrowse
  populate-dscensor  CLI entry for populate-dscensor
  populate-jbrowse2  CLI entry for populate-jbrowse2
  populate-jekyll    CLI entry for populate-jekyll
```

## Building the catalog

`populate-catalog` describes the **whole** datastore in one JSON document, built
from a `datastore-metadata` checkout with **no network access and no taxa list**:

```
(lis_autocontent_env) $ lis-autocontent populate-catalog --from_github ./datastore-metadata \
                            --catalog_out ./autocontent/catalog.json
```

It joins the four metadata layers the datastore specification guarantees --
`README` (provenance, DOI, taxonomy, assembly conventions), `CHECKSUM` (the file
list), `MANIFEST` (per-file descriptions and application tags) and the committed
BUSCO summaries (completeness *and* assembly counts) -- then links annotations to
the genomes they derive from and inherits assembly conventions along that edge.

The reading is done by `DatastoreIndex`
([`src/lis_autocontent/datastore_files.py`](src/lis_autocontent/datastore_files.py)),
an in-memory index of every collection and file in the checkout that Python callers
can query directly; `populate-catalog` serializes it.

Against the full store this is roughly 1,000 collections and 5,000 files in about
a second. Feed the result to DSCensor with `dscensor --catalog ./autocontent/catalog.json`.

### Collections that publish no CHECKSUM

Roughly 45% of collections (nearly every `qtl`, `gwas` and `maps` one) publish no
CHECKSUM, and CHECKSUM is the only enumeration that lists a collection's files. For
those, `populate-catalog` builds the file list from the datastore's own documented
convention -- `{abbrev}.{identifier}.{content}.{ext}`, with the content vocabulary
vendored in [`src/lis_autocontent/filetypes.yml`](src/lis_autocontent/filetypes.yml) from
[datastore-specifications](https://github.com/legumeinfo/datastore-specifications).

Add `--verify` to confirm each constructed filename with a HEAD request before it
enters the catalog:

```
lis-autocontent populate-catalog --from_github ./datastore-metadata --verify
```

The difference is recorded, never assumed. Every file carries a `src`, and every
collection an `index_status`:

| `index_status` | file `src` | meaning |
|---|---|---|
| `known` | `checksum` | listed in CHECKSUM; authoritative |
| `verified` | `verified` | built from the convention, HEAD-confirmed to exist |
| `inferred` | `predicted` | built from the convention, not checked |
| `unknown` | — | no CHECKSUM and no documented convention |

Measured on the full store: 1,651 predicted filenames, 1,283 confirmed by `--verify`
in ~85s. Consumers (DSCensor, the legumista MCP server, JBrowse/BLAST config) read
these labels and report them; **none of them re-derives filenames**, so the convention
lives in one place.

Two properties worth preserving if you change it:

* **The build never touches the network unless you pass `--verify`.** That is what lets
  it re-run on every metadata push. Assembly metrics come from the committed BUSCO summaries, which
  already carry `Scaffold N50`, `Number of contigs` and `Total length`; no `.fai`
  fetch is required.
* **`index_status` is tri-state.** `"known"` means the collection published a
  CHECKSUM and its file list is authoritative. `"unknown"` means it did not, so an
  empty file list says nothing about whether indexes exist. Roughly 45% of
  collections (nearly every `qtl` and `gwas` one) are `"unknown"`. Collapsing that
  into "no indexes" would report streamable data as unreadable.

### Building it on GitHub Actions

`.github/workflows/populate-catalog.yml` builds `catalog.json` for the whole store from
the current `datastore-metadata` `main`, publishes it as a GitHub release asset, and then
notifies an external endpoint. Trigger it from the repository's Actions tab, or on the
command line:

```
gh workflow run populate-catalog.yml --repo legumeinfo/LIS-autocontent
```

To trigger it from another repository's workflow, dispatch it with a token that has
`actions:write` on this repository (a fine-grained token needs the Actions permission, a
classic one the `repo` scope). The calling repository's own `GITHUB_TOKEN` only grants
access to the repository it runs in, so it cannot dispatch this workflow:

```yaml
- name: Rebuild the LIS catalog
  env:
    GH_TOKEN: ${{ secrets.AUTOCONTENT_DISPATCH_TOKEN }}
  run: gh workflow run populate-catalog.yml --repo legumeinfo/LIS-autocontent
```

The workflow has to exist on the default branch to be dispatchable at all; it then runs
against whichever ref the dispatch names (`--ref`, defaulting to the default branch).

The workflow builds offline, without `--verify`, so the same `datastore-metadata` commit
always gives the same catalog: a probe that fails for any reason reads as "file absent"
and would drop that file from what gets published.

Each run publishes its own release, tagged with the build time and the
`datastore-metadata` commit that determined the content (for example
`catalog-20260918-120435-04d9d86`), carrying `catalog.json` as its only asset. Every past
catalog stays fetchable at its own tag, and the newest one is always at a stable URL:

```
https://github.com/legumeinfo/LIS-autocontent/releases/latest/download/catalog.json
```

Unlike a workflow artifact, a release asset downloads without authentication, so any
consumer can fetch it with plain `curl` from a public repository. The run needs
`contents: write` to create the release. Seconds are part of the tag so that a re-run
against the same `datastore-metadata` commit gets its own release: `gh release create`
fails on an existing tag rather than overwriting a published catalog.

Once the release is published, the workflow POSTs JSON to the URL in the `WEBHOOK_URL`
repository variable (default `https://example.com`, which rejects POSTs, so runs fail at
that step until the variable is set):

```json
{
  "event": "catalog-built",
  "repository": "legumeinfo/LIS-autocontent",
  "commit": "<LIS-autocontent commit>",
  "datastore_metadata_commit": "<datastore-metadata commit>",
  "run_url": "https://github.com/legumeinfo/LIS-autocontent/actions/runs/<run id>",
  "release": {
    "tag": "catalog-20260918-120435-04d9d86",
    "url": "https://github.com/legumeinfo/LIS-autocontent/releases/tag/<tag>",
    "asset": {
      "name": "catalog.json",
      "download_url": "https://github.com/legumeinfo/LIS-autocontent/releases/download/<tag>/catalog.json",
      "latest_download_url": "https://github.com/legumeinfo/LIS-autocontent/releases/latest/download/catalog.json",
      "size": 2411008,
      "sha256": "<digest of the published file>"
    }
  }
}
```

The receiver fetches `release.asset.download_url` directly, with no token. If the
`WEBHOOK_SECRET` secret is set, the body is signed the way GitHub signs its webhooks:
`X-Hub-Signature-256: sha256=<HMAC-SHA256 of the raw body, keyed by the secret>`. A failed
build fails the run before anything is published or sent. A webhook the endpoint rejects
also fails the run, but the release is already published by then and stays published;
only timeouts, refused connections and 408/429/5xx responses are retried.

## Building a Divbrowse compose file

`populate-divbrowse` writes a [Divbrowse](https://github.com/legumeinfo/divbrowse)
`docker-compose.yml` with one service per diversity collection, offline from a
`datastore-metadata` checkout:

```
lis-autocontent populate-divbrowse --from_github ./datastore-metadata \
    --collection Wm82.gnm4.div.Song_Hyten_2015 \
    --collection Wm82.gnm5.div.Song_Hyten_2015 \
    --collection Wm82.gnm6.div.Song_Hyten_2015 \
    --compose_out ../divbrowse/docker-compose.yml
```

Each service gets the collection's VCF, the `gene_models_main` GFF3 of its assembly's
annotation, a `CHROM_PATTERN` built from the genome's `chromosome_prefix`, and a
`BASE_URL` under `--base_url` (default `https://divbrowse.soybase.org`); host ports start
at `--port` (default 8080). The services build from the Divbrowse `Dockerfile`, so write
the file to the root of a divbrowse checkout. Service names are the collection
identifier in lowercase, because Docker image names must be lowercase; data directories
keep the identifier as-is.

A collection this format can't express stops the command with the reason, and nothing
is written: several VCFs or none, no linked genome assembly, more than one annotation,
or a `chromosome_prefix` that isn't a single prefix.

### Building it on GitHub Actions

`.github/workflows/populate-divbrowse.yml` builds the compose file for the collections
listed in [`.github/divbrowse-collections.txt`](.github/divbrowse-collections.txt) (one
identifier per line, in service order) against the current `datastore-metadata` `main`,
stores it as a workflow artifact, and then notifies an external endpoint. Trigger it from
the repository's Actions tab, or on the command line:

```
gh workflow run populate-divbrowse.yml --repo legumeinfo/LIS-autocontent
```

Triggering it from another repository's workflow works the same way as for the catalog
above, and needs the same `actions:write` token.

Once the artifact is uploaded, the workflow POSTs JSON to the URL in the
`WEBHOOK_URL` repository variable (default `https://example.com`, which rejects
POSTs, so runs fail at that step until the variable is set):

```json
{
  "event": "divbrowse-compose-built",
  "repository": "legumeinfo/LIS-autocontent",
  "commit": "<LIS-autocontent commit>",
  "datastore_metadata_commit": "<datastore-metadata commit>",
  "run_url": "https://github.com/legumeinfo/LIS-autocontent/actions/runs/<run id>",
  "collections": ["Wm82.gnm4.div.Song_Hyten_2015", "..."],
  "artifact": {
    "id": 1234,
    "name": "docker-compose.yml",
    "url": "https://github.com/legumeinfo/LIS-autocontent/actions/runs/<run id>/artifacts/1234",
    "download_url": "https://api.github.com/repos/legumeinfo/LIS-autocontent/actions/artifacts/1234/zip",
    "sha256": "<digest of the uploaded file>"
  }
}
```

GitHub artifacts can't be downloaded anonymously, so the receiver fetches
`artifact.download_url` with its own GitHub token that can read this repository;
`artifact.url` is the same artifact for a signed-in browser. If the
`WEBHOOK_SECRET` secret is set, the body is signed the way GitHub signs its
webhooks: `X-Hub-Signature-256: sha256=<HMAC-SHA256 of the raw body, keyed by the secret>`.

A collection that can't be built fails the run before anything is uploaded or sent. A
webhook the endpoint rejects also fails the run; only timeouts, refused connections and
408/429/5xx responses are retried.
