"""Build a single catalog artifact describing the whole LIS Data Store.

This is the "pass 1" half of the census: everything derivable from the
datastore-metadata checkout alone, with **no network access**. It reads the four
metadata layers the datastore specification guarantees --

    README.<collection>.yml     provenance, DOI, taxonomy, assembly conventions
    CHECKSUM.<collection>.md5   the authoritative file list
    MANIFEST.<collection>.yml   a description (and application tags) per file
    BUSCO/*.short_summary.json  completeness scores and assembly counts

-- and emits one JSON document that downstream consumers (DSCensor, the
legumista MCP server, any future client) can hold in memory.

Why a separate module from ``process_collections``: that class scrapes the store
to drive JBrowse/BLAST/Jekyll and needs the network to do it. The catalog needs
none, so welding the two together would make a sub-second job inherit a
multi-minute one. Keeping them apart lets the catalog rebuild on every metadata
push while enrichment runs on its own cadence.

Two details worth knowing before changing this file:

* ``CHECKSUM`` is the *only* enumeration anywhere that lists the ``.fai``/``.tbi``
  index siblings -- the store's HTML index and its JSON API both filter those
  extensions out. That is why file-level random-access flags come from here and
  can come from nowhere else.
* Roughly 45% of collections publish no CHECKSUM at all (nearly every ``qtl`` and
  ``gwas`` one). Those get ``index_status: "unknown"`` and an empty file list --
  never ``"none"``. Collapsing "we looked and there is no index" into "we could
  not look" is a correctness bug, not a cosmetic one.
"""

import json
import logging
import os
import subprocess
import re
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import yaml

SCHEMA_VERSION = 1

# Vendored copy of the datastore file-content vocabulary; see filetypes.yml.
FILETYPES_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "filetypes.yml"
)

# How a file entry came to be known. Recorded per file so consumers (DSCensor, the
# legumista MCP server, JBrowse/BLAST config) can report provenance without
# re-deriving it -- they stay dumb readers of what this build decided.
SRC_CHECKSUM = "checksum"  # listed in CHECKSUM.md5; authoritative
SRC_VERIFIED = "verified"  # constructed by convention, HEAD-confirmed to exist
SRC_PREDICTED = "predicted"  # constructed by convention, not checked

PROBE_TIMEOUT = int(os.environ.get("LIS_CATALOG_PROBE_TIMEOUT", "20"))
PROBE_WORKERS = int(os.environ.get("LIS_CATALOG_PROBE_WORKERS", "8"))

# Index siblings, mapped to the reader that consumes the file they belong to.
INDEX_SUFFIXES = (".fai", ".tbi", ".csi", ".bai", ".crai")
# Companions that are not themselves data and never determine access on their own.
COMPANION_SUFFIXES = (".gzi", ".md5")
METADATA_PREFIXES = ("README.", "MANIFEST.", "CHANGES.", "CHECKSUM.", "LINKOUTS.")

# README fields copied verbatim onto a collection. Everything here is published
# by the datastore spec; the point of the list is that nothing gets silently
# dropped the way `taxid` was before.
README_FIELDS = (
    "synopsis",
    "description",
    "scientific_name",
    "taxid",
    "scientific_name_abbrev",
    "genotype",
    "publication_doi",
    "publication_title",
    "citation",
    "license",
    "public_access_level",
    "keywords",
    "chromosome_prefix",
    "supercontig_prefix",
    "bioproject",
    "genbank_accession",
    "sraproject",
    "genetic_map",
    "expression_unit",
    "dataset_doi",
    "related_to",
    "source",
    "dataset_release_date",
)

# The subset of README_FIELDS worth carrying on a per-file DSCensor node. Nodes are
# per-file and numerous, so the long prose (description, citation, provenance) stays in
# the catalog and out of the nodes. Shared with process_collections so the two paths
# cannot drift apart.
NODE_README_FIELDS = (
    "taxid",
    "scientific_name_abbrev",
    "publication_doi",
    "publication_title",
    "license",
    "chromosome_prefix",
    "supercontig_prefix",
    "bioproject",
    "genetic_map",
    "expression_unit",
    "genotype",
)

# Fields an annotation may inherit from the genome it derives from. These live on
# the genome README only, so a flat node list cannot express them -- following the
# derived_from edge is the whole reason the catalog is a graph and not a table.
INHERITED_FIELDS = ("chromosome_prefix", "supercontig_prefix", "bioproject")


class CatalogBuilder:
    """Assemble the catalog from a datastore-metadata checkout."""

    def __init__(self, metadata_dir, logger=None, datastore_url=None, verify=False):
        self.metadata_dir = os.path.abspath(metadata_dir)
        self.logger = logger or logging.getLogger("catalog")
        self.datastore_url = (datastore_url or "https://data.legumeinfo.org").rstrip(
            "/"
        )
        self.collections = []
        self._busco_available = 0
        # When True, every predicted filename is confirmed with a HEAD request before
        # it enters the catalog. Off by default so a build stays offline and fast;
        # CI turns it on. The difference is visible in each file's `src`, never silent.
        self.verify = verify
        self.filetypes = self._load_filetypes()
        self.stats = {
            "collections": 0,
            "files": 0,
            "indexed_files": 0,
            "with_doi": 0,
            "index_unknown": 0,
            "with_busco": 0,
            "curated_symbols": 0,
            "described_taxa": 0,
            "pairwise_files": 0,
            "paired_genomes": 0,
            "predicted_files": 0,
            "verified_files": 0,
            "probes": 0,
        }

    # ------------------------------------------------------------------ helpers
    def _read_text(self, path):
        """Read a metadata file, or return "" when it is absent or unreadable.

        Metadata files are individually optional in practice, so a miss must
        degrade the record rather than abort the build.
        """
        try:
            with open(path, encoding="utf-8", errors="replace") as handle:
                return handle.read()
        except OSError:
            return ""

    def _read_yaml(self, path):
        """Parse a datastore YAML file into its first document, or {}."""
        text = self._read_text(path)
        if not text.strip():
            return {}
        try:
            docs = [
                doc
                for doc in yaml.safe_load_all(text)
                if isinstance(doc, dict)  # skip stray list/None documents
            ]
        except yaml.YAMLError as err:
            self.logger.warning("could not parse %s: %s", path, err)
            return {}
        return docs[0] if docs else {}

    def source_commit(self):
        """The datastore-metadata commit this catalog was built from.

        Recorded so a consumer can tell how stale its copy is. A checkout sitting
        hundreds of commits behind produces a catalog that under-reports what is
        streamable, and the only defence is making the provenance visible.
        """
        try:
            out = subprocess.run(
                ["git", "-C", self.metadata_dir, "rev-parse", "HEAD"],
                capture_output=True,
                check=False,
                timeout=30,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        if out.returncode != 0:
            return None
        return out.stdout.decode("utf-8", "replace").strip() or None

    def _load_filetypes(self):
        """The vendored content vocabulary, or {} if it cannot be read.

        A missing vocabulary degrades to the previous behaviour (no prediction), it
        does not fail the build.
        """
        try:
            with open(FILETYPES_PATH, encoding="utf-8") as handle:
                loaded = yaml.safe_load(handle)
        except (OSError, yaml.YAMLError) as err:
            self.logger.warning("could not read %s: %s", FILETYPES_PATH, err)
            return {}
        return loaded if isinstance(loaded, dict) else {}

    def predicted_files(self, ctype, identifier, abbrev):
        """Filenames a collection of this type is expected to publish.

        Returns [] when the type has no documented vocabulary or is explicitly
        unpredictable. Prediction is a default, not a guarantee -- callers must label
        the result (see SRC_*).
        """
        spec = self.filetypes.get(ctype)
        if not isinstance(spec, dict) or spec.get("predictable") is False:
            return []
        contents = spec.get("contents") or []
        extension = spec.get("extension") or ""
        if spec.get("prefix") == "none":
            stem = identifier
        elif abbrev:
            stem = "{0}.{1}".format(abbrev, identifier)
        else:
            return []  # abbrev-prefixed type with no abbrev: nothing safe to build
        return ["{0}.{1}{2}".format(stem, content, extension) for content in contents]

    def type_is_indexed(self, ctype):
        """Whether this collection type ever publishes index siblings."""
        spec = self.filetypes.get(ctype)
        return not (isinstance(spec, dict) and spec.get("indexed") is False)

    def url_exists(self, url):
        """HEAD one datastore URL.

        The path is percent-encoded first: 21 collections carry non-ASCII names
        (Nicolás, Valdés-López, Cortés, ...), and urllib raises UnicodeEncodeError on
        those rather than returning a status -- which, caught below, would silently
        report every one of their files as absent.
        """
        request = urllib.request.Request(
            urllib.parse.quote(url, safe=":/?#[]@!$&'()*+,;="),
            method="HEAD",
            headers={"User-Agent": "lis-autocontent/1.0"},
        )
        try:
            with urllib.request.urlopen(request, timeout=PROBE_TIMEOUT) as response:
                return 200 <= getattr(response, "status", 200) < 300
        except Exception:  # noqa: BLE001 - a 404 (or any failure) means "not there"
            return False

    def attach_indexes(self, base_url, files):
        """Probe each file for the index siblings its extension allows.

        Only reached under --verify, and only for types that can carry indexes, so the
        cost is a few dozen requests. Without it those files would have to be reported
        as index-status-unknown.
        """
        jobs = []
        for entry in files:
            for suffix in self.plausible_indexes(entry["n"]):
                jobs.append((entry, suffix))
        if not jobs:
            return
        urls = ["{0}/{1}{2}".format(base_url.rstrip("/"), e["n"], s) for e, s in jobs]
        with ThreadPoolExecutor(max_workers=PROBE_WORKERS) as pool:
            found = list(pool.map(self.url_exists, urls))
        self.stats["probes"] += len(urls)
        for (entry, suffix), ok in zip(jobs, found):
            if ok:
                entry.setdefault("i", []).append(suffix)

    @staticmethod
    def plausible_indexes(name):
        """Index suffixes worth probing for this filename.

        A file's own type decides which siblings are even possible: .fai pairs with
        FASTA, .bai/.csi with BAM, .crai with CRAM, .tbi/.csi with a bgzipped tabbed
        file. Asking a .tsv.gz about a .bai is a guaranteed 404.
        """
        low = name.lower()
        if low.endswith(
            (
                ".fa",
                ".fa.gz",
                ".fasta",
                ".fasta.gz",
                ".fna",
                ".fna.gz",
                ".faa",
                ".faa.gz",
            )
        ):
            return (".fai",)
        if low.endswith(".bam"):
            return (".bai", ".csi")
        if low.endswith(".cram"):
            return (".crai",)
        if low.endswith(
            (
                ".vcf.gz",
                ".bcf",
                ".gff.gz",
                ".gff3.gz",
                ".gtf.gz",
                ".bed.gz",
                ".sam.gz",
                ".tsv.gz",
                ".txt.gz",
            )
        ):
            return (".tbi", ".csi")
        return INDEX_SUFFIXES

    def verify_files(self, base_url, names):
        """Keep the names that actually exist. Probes run concurrently."""
        if not names:
            return []
        urls = ["{0}/{1}".format(base_url.rstrip("/"), n) for n in names]
        with ThreadPoolExecutor(max_workers=PROBE_WORKERS) as pool:
            found = list(pool.map(self.url_exists, urls))
        self.stats["probes"] += len(urls)
        return [n for n, ok in zip(names, found) if ok]

    # ------------------------------------------------------------------ parsing
    def checksum_files(self, coll_dir, identifier):
        """File names from CHECKSUM, or None when the collection publishes none.

        ``None`` and ``[]`` mean different things here and callers rely on it:
        None is "we could not look", [] is "we looked and it is empty".
        """
        path = os.path.join(coll_dir, "CHECKSUM.{0}.md5".format(identifier))
        if not os.path.isfile(path):
            return None
        names = []
        for line in self._read_text(path).splitlines():
            parts = line.split(None, 1)
            if len(parts) == 2:
                names.append(parts[1].strip().lstrip("./"))
        return names

    def manifest_descriptions(self, coll_dir, identifier):
        """Map file name -> {description, applications} from MANIFEST."""
        path = os.path.join(coll_dir, "MANIFEST.{0}.yml".format(identifier))
        text = self._read_text(path)
        if not text.strip():
            return {}
        try:
            loaded = list(yaml.safe_load_all(text))
        except yaml.YAMLError as err:
            self.logger.warning("could not parse %s: %s", path, err)
            return {}
        out = {}
        for doc in loaded:
            entries = doc if isinstance(doc, list) else [doc]
            for entry in entries:
                if not isinstance(entry, dict) or not entry.get("name"):
                    continue
                record = {}
                description = str(entry.get("description") or "").strip()
                # "MISSING" is the spec's placeholder for an undescribed file.
                if description and description != "MISSING":
                    record["description"] = description
                applications = entry.get("applications")
                if applications:
                    record["applications"] = [str(a) for a in applications]
                if record:
                    out[str(entry["name"])] = record
        return out

    def busco_summary(self, coll_dir):
        """BUSCO scores and assembly counts, both from the committed summary.

        The short_summary JSON carries ``Scaffold N50``/``Number of contigs``/
        ``Total length`` alongside the completeness scores, so assembly metrics
        need no ``.fai`` fetch and pass 1 stays fully offline.
        """
        busco_dir = os.path.join(coll_dir, "BUSCO")
        if not os.path.isdir(busco_dir):
            return None, None
        summary = None
        for name in sorted(os.listdir(busco_dir)):
            if name.endswith(".json") and "short_summary" in name:
                summary = os.path.join(busco_dir, name)
                break
        if not summary:
            return None, None
        try:
            with open(summary, encoding="utf-8") as handle:
                results = json.load(handle).get("results", {})
        except (OSError, ValueError) as err:
            self.logger.warning("could not parse %s: %s", summary, err)
            return None, None
        if not results:
            return None, None

        markers = float(results.get("n_markers") or 0)
        busco = {
            "lineage": results.get("domain"),
            "total": int(markers),
            "complete_pct": results.get("Complete"),
            "single_copy_pct": results.get("Single copy"),
            "duplicate_pct": results.get("Multi copy"),
            "fragmented_pct": results.get("Fragmented"),
            "missing_pct": results.get("Missing"),
        }
        counts = {}
        for key, field in (
            ("scaffolds", "Number of scaffolds"),
            ("contigs", "Number of contigs"),
            ("total_length", "Total length"),
            ("scaffold_n50", "Scaffold N50"),
            ("contig_n50", "Contigs N50"),
        ):
            raw = results.get(field)
            if raw in (None, ""):
                continue
            try:
                counts[key] = int(raw)
            except (TypeError, ValueError):
                continue
        gaps = results.get("Percent gaps")
        if isinstance(gaps, str) and gaps.endswith("%"):
            try:
                counts["percent_gaps"] = float(gaps.rstrip("%"))
            except ValueError:
                pass
        return busco, (counts or None)

    def classify_files(self, names):
        """Split a CHECKSUM listing into data files, each tagged with its indexes."""
        present = set(names)
        files = []
        for name in sorted(names):
            base = os.path.basename(name)
            if name.endswith(INDEX_SUFFIXES + COMPANION_SUFFIXES):
                continue
            if base.startswith(METADATA_PREFIXES):
                continue
            indexes = [s for s in INDEX_SUFFIXES if name + s in present]
            record = {"n": name}
            if indexes:
                record["i"] = indexes
            files.append(record)
        return files

    # ------------------------------------------------------------------ building
    def collection_record(self, coll_dir, genus, species, ctype, identifier, readme):
        """Assemble one collection record from its metadata files."""
        names = self.checksum_files(coll_dir, identifier)
        manifest = self.manifest_descriptions(coll_dir, identifier)
        busco, counts = self.busco_summary(coll_dir)

        if names is not None:
            files = self.classify_files(names)
            for entry in files:
                entry["src"] = SRC_CHECKSUM
            index_status = "known"
        else:
            # No CHECKSUM: build the file list from the documented convention instead.
            base_url = "{0}/{1}/{2}/{3}/{4}".format(
                self.datastore_url, genus, species, ctype, identifier
            )
            predicted = self.predicted_files(
                ctype, identifier, readme.get("scientific_name_abbrev")
            )
            if not predicted:
                files, index_status = [], "unknown"
            elif self.verify:
                predicted = self.verify_files(base_url, predicted)
                files = [{"n": n, "src": SRC_VERIFIED} for n in predicted]
                index_status = "verified"
                self.stats["verified_files"] += len(files)
            else:
                files = [{"n": n, "src": SRC_PREDICTED} for n in predicted]
                index_status = "inferred"
                self.stats["predicted_files"] += len(files)
            # Types documented as never carrying index siblings need no probing to know
            # none of these files is randomly accessible. For the types that CAN carry
            # them, an unprobed file's index status is genuinely unknown -- reporting it
            # as "not indexed" would mark streamable data unreadable, which is the same
            # mistake `index_status` exists to prevent. Verified: a markers .gff3.gz
            # constructed this way does have a .tbi.
            if files and self.type_is_indexed(ctype):
                if self.verify:
                    self.attach_indexes(base_url, files)
                else:
                    for entry in files:
                        entry["i_unknown"] = True

        for record in files:
            extra = manifest.get(os.path.basename(record["n"]))
            if extra:
                record.update(extra)

        record = {
            "path": "{0}/{1}/{2}/{3}".format(genus, species, ctype, identifier),
            "id": identifier,
            "type": ctype,
            "genus": genus,
            "species": species,
            "base_url": "{0}/{1}/{2}/{3}/{4}".format(
                self.datastore_url, genus, species, ctype, identifier
            ),
            "index_status": index_status,
            "files": files,
        }
        for field in README_FIELDS:
            value = readme.get(field)
            if value not in (None, "", [], {}):
                record[field] = value
        # NOTE: assembly metrics are parsed but deliberately NOT emitted yet.
        # `busco_summary()` above still ingests both the completeness scores and the
        # assembly counts from the committed BUSCO summary -- that code is live and
        # tested, and needs no network. Only the inclusion in the catalog is deferred,
        # so re-enabling it is uncommenting these four lines. Anything downstream that
        # consumes `busco`/`counts` is commented out with a matching NOTE.
        # if busco:
        #     record["busco"] = busco
        # if counts:
        #     record["counts"] = counts
        if busco:
            self._busco_available += 1  # build diagnostic; independent of emission
        del busco, counts  # parsed above; see NOTE -- remove when metrics are enabled
        return record

    def pairwise_relationships(self):
        """Genome pairs, derived from the filenames of synteny and alignment files.

        A pairwise file is named ``<A>.x.<B>[.<epoch>].<KEY>.<ext>`` and is stored ONCE,
        under whichever genome is the reference. So a genome's relationships live in two
        places: its own collection (as A) and other species' collections (as B). A
        consumer that reads only a genome's own collection silently misses half its
        synteny -- and, for a genome with no collection of its own (Medicago has none),
        all of it.

        Deriving the graph here means every consumer gets it for free and none of them
        has to re-implement the filename parse. Self-comparisons carry a whole-genome-
        duplication epoch instead of a partner (``...x.glyma.Wm82.gnm2.old_duplication``);
        the epoch is kept rather than collapsed, because which duplication a block came
        from is the informative part.
        """
        pattern = re.compile(
            r"^(?P<a>[a-z]{4,6}\.[A-Za-z0-9_-]+\.gnm\d+)"
            r"\.x\."
            r"(?P<b>[a-z]{4,6}\.[A-Za-z0-9_-]+\.gnm\d+)"
            r"(?P<epoch>\.[A-Za-z0-9_]+)?"
            r"\.(?P<key>[A-Za-z0-9]{4})"
            r"\.(?P<ext>gff3\.gz|paf\.gz|bam)$"
        )
        pairs = []
        for record in self.collections:
            if record["type"] not in ("synteny", "genome_alignments"):
                continue
            for entry in record.get("files", []):
                match = pattern.match(entry["n"])
                if not match:
                    continue
                epoch = (match.group("epoch") or "").lstrip(".")
                pair = {
                    "a": match.group("a"),
                    "b": match.group("b"),
                    "kind": "synteny" if record["type"] == "synteny" else "alignment",
                    "format": match.group("ext"),
                    "collection": record["path"],
                    "file": entry["n"],
                    "url": "{0}/{1}".format(record["base_url"].rstrip("/"), entry["n"]),
                }
                if epoch:
                    pair["epoch"] = epoch  # whole-genome-duplication label
                if match.group("a") == match.group("b"):
                    pair["self"] = True
                if entry.get("i"):
                    pair["i"] = entry["i"]
                pairs.append(pair)
        pairs.sort(key=lambda p: (p["a"], p["b"], p["file"]))
        return pairs

    def descriptions(self):
        """Species- and genus-level metadata from ``about_this_collection``.

        ``description_<Genus>[_<species>].yml`` carries the common name, NCBI taxid,
        the abbreviation the datastore uses everywhere, a prose description, and links
        to related resources (mines, browsers). None of it is on a collection README,
        so without this a consumer cannot answer "what is soybean called here" or map
        a common name onto a taxon at all.
        """
        out = {}
        for dirpath, _dirnames, filenames in os.walk(self.metadata_dir):
            if os.sep + ".git" in dirpath or not dirpath.endswith(
                "about_this_collection"
            ):
                continue
            for name in sorted(filenames):
                if not (name.startswith("description_") and name.endswith(".yml")):
                    continue
                doc = self._read_yaml(os.path.join(dirpath, name))
                if not doc:
                    continue
                genus = str(doc.get("genus") or "").strip()
                if not genus:
                    continue
                # A species-level file names one species; a genus-level file lists the
                # genus's species instead. Key them apart rather than conflating.
                raw_species = doc.get("species")
                species = (
                    str(raw_species).strip()
                    if isinstance(raw_species, (str, int))
                    else ""
                )
                key = f"{genus}/{species}" if species else genus
                record = {}
                for field in ("taxid", "abbrev", "commonName", "description", "genus"):
                    value = doc.get(field)
                    if value not in (None, "", [], {}):
                        record[field] = value
                if species:
                    record["species"] = species
                elif isinstance(raw_species, list) and raw_species:
                    record["species_in_genus"] = [str(x) for x in raw_species]
                resources = doc.get("resources")
                if isinstance(resources, list) and resources:
                    record["resources"] = [
                        {
                            k: r.get(k)
                            for k in ("name", "URL", "description")
                            if r.get(k)
                        }
                        for r in resources
                        if isinstance(r, dict)
                    ]
                out[key] = record
        return out

    def curated_symbols(self):
        """Curated gene symbol -> gene id, from ``gene_functions/<abbrev>.traits.yml``.

        These directories carry no README, so the collection walk skips them, but the
        traits files ARE tracked in datastore-metadata. Folding them in lets a consumer
        resolve a symbol like ``GmNARK`` without a network round-trip, and the file
        also carries the gene's own publication DOI.
        """
        symbols = {}
        for dirpath, _dirnames, filenames in os.walk(self.metadata_dir):
            if os.sep + ".git" in dirpath or not dirpath.endswith("gene_functions"):
                continue
            for name in sorted(filenames):
                if not name.endswith(".traits.yml"):
                    continue
                abbrev = name[: -len(".traits.yml")]
                text = self._read_text(os.path.join(dirpath, name))
                try:
                    docs = list(yaml.safe_load_all(text))
                except yaml.YAMLError as err:
                    self.logger.warning("could not parse %s: %s", name, err)
                    continue
                index = symbols.setdefault(abbrev, {})
                for doc in docs:
                    if not isinstance(doc, dict):
                        continue
                    gene = (doc.get("gene_model_full_id") or "").strip()
                    if not gene:
                        continue
                    dois = [
                        ref.get("doi")
                        for ref in (doc.get("references") or [])
                        if isinstance(ref, dict) and ref.get("doi")
                    ]
                    entry = {"gene": gene}
                    if dois:
                        entry["doi"] = dois[0]
                    synopsis = (doc.get("phenotype_synopsis") or "").strip()
                    if synopsis:
                        entry["synopsis"] = synopsis
                    for symbol in doc.get("gene_symbols") or []:
                        index[str(symbol).strip().lower()] = entry
        return {k: v for k, v in symbols.items() if v}

    def walk(self):
        """Find every collection directory and build its record.

        A collection qualifies on *any* metadata file, not just a README. Some
        collections publish only a CHECKSUM -- five genome_alignments directories do,
        and they hold indexed BAMs -- so keying off the README alone silently dropped
        them along with their randomly accessible files.
        """
        root = self.metadata_dir
        for dirpath, _dirnames, filenames in os.walk(root):
            if os.sep + ".git" in dirpath:
                continue
            relative = os.path.relpath(dirpath, root).split(os.sep)
            if len(relative) != 4:  # Genus/species/type/collection
                continue
            genus, species, ctype, identifier = relative
            readme_file = self._readme_filename(filenames)
            has_metadata = readme_file or any(
                f.startswith(("CHECKSUM.", "MANIFEST.")) for f in filenames
            )
            if not has_metadata:
                continue
            readme = {}
            if readme_file:
                readme = self._read_yaml(os.path.join(dirpath, readme_file))
                if not readme:
                    self.logger.warning("unreadable README in %s", dirpath)
            self.collections.append(
                self.collection_record(
                    dirpath, genus, species, ctype, identifier, readme
                )
            )

    @staticmethod
    def _readme_filename(filenames):
        """The collection's README, tolerating the extension-less form.

        At least one collection ships a bare ``README`` rather than
        ``README.<collection>.yml``; the contents are YAML either way.
        """
        for name in sorted(filenames):
            if name.startswith("README.") and name.endswith(".yml"):
                return name
        return "README" if "README" in filenames else None

    def link_lineage(self):
        """Add ``derived_from`` edges and inherit assembly conventions along them.

        An annotation identifier embeds the assembly it was called against
        (``Wm82.gnm4.ann1.T8TQ`` -> ``Wm82.gnm4``), so the edge is derivable
        without any extra metadata. ``chromosome_prefix`` and friends live on the
        genome README only, which is what makes the inheritance necessary.
        """
        genomes = {}
        for record in self.collections:
            if record["type"] != "genomes":
                continue
            stem = ".".join(record["id"].split(".")[:2])  # strain.gnmN
            genomes.setdefault((record["genus"], record["species"], stem), record)

        for record in self.collections:
            parts = record["id"].split(".")
            if len(parts) < 3 or not parts[1].startswith("gnm"):
                continue
            stem = ".".join(parts[:2])
            parent = genomes.get((record["genus"], record["species"], stem))
            if parent is None or parent is record:
                continue
            record["derived_from"] = [parent["id"]]
            for field in INHERITED_FIELDS:
                if field not in record and field in parent:
                    record[field] = parent[field]
                    record.setdefault("inherited", []).append(field)

    def build(self):
        """Build the whole catalog. Returns the document as a dict."""
        self.collections = []
        self._busco_available = 0
        self.walk()
        self.link_lineage()
        self.collections.sort(key=lambda record: record["path"])

        # with_busco is a build diagnostic counted at parse time, so it stays true
        # while the metrics themselves are withheld from the record (see NOTE above).
        self.stats["with_busco"] = self._busco_available
        self.stats["collections"] = len(self.collections)
        for record in self.collections:
            self.stats["files"] += len(record["files"])
            self.stats["indexed_files"] += sum(1 for f in record["files"] if "i" in f)
            if record.get("publication_doi"):
                self.stats["with_doi"] += 1
            if record["index_status"] == "unknown":
                self.stats["index_unknown"] += 1

        symbols = self.curated_symbols()
        self.stats["curated_symbols"] = sum(len(v) for v in symbols.values())
        descriptions = self.descriptions()
        self.stats["described_taxa"] = len(descriptions)
        pairs = self.pairwise_relationships()
        self.stats["pairwise_files"] = len(pairs)
        self.stats["paired_genomes"] = len({g for p in pairs for g in (p["a"], p["b"])})
        return {
            "schema": SCHEMA_VERSION,
            "built_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "source_commit": self.source_commit(),
            "datastore_url": self.datastore_url,
            "stats": dict(self.stats),
            "taxa": descriptions,
            "pairwise": pairs,
            "gene_symbols": symbols,
            "collections": self.collections,
        }

    def write(self, out_path, indent=None):
        """Build and write the catalog. Returns the document."""
        catalog = self.build()
        parent = os.path.dirname(os.path.abspath(out_path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        separators = None if indent else (",", ":")
        with open(out_path, "w", encoding="utf-8") as handle:
            json.dump(catalog, handle, indent=indent, separators=separators)
        self.logger.info(
            "wrote %s (%s collections, %s files, %s indexed, %s without CHECKSUM)",
            out_path,
            self.stats["collections"],
            self.stats["files"],
            self.stats["indexed_files"],
            self.stats["index_unknown"],
        )
        return catalog
