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
6. Build package locally. `python setup.py build` (optional)

## Run

```
(lis_autocontent_env) $ lis-autocontent --help
Usage: lis-autocontent [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  populate-blast     CLI entry for populate-blast
  populate-catalog   CLI entry for populate-catalog
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
