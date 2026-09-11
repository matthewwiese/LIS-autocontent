"""Index every collection and file in the datastore from a datastore-metadata checkout.

Each collection directory carries a CHECKSUM.<key>.md5 listing every file it
contains, including files in subdirectories such as BUSCO/. Treating that
manifest as the authoritative file list removes the need to construct URLs from
hardcoded naming conventions, and surfaces files the conventions never covered.
It is also the *only* enumeration that lists the ``.fai``/``.tbi`` index siblings --
the store's HTML index and its JSON API both filter those out -- so file-level
random-access flags can come from nowhere else.

The index covers the whole store, with no network access. It reads the four
metadata layers the datastore specification guarantees --

    README.<key>.yml            provenance, DOI, taxonomy, assembly conventions
    CHECKSUM.<key>.md5          the authoritative file list
    MANIFEST.<key>.yml          a description (and application tags) per file
    BUSCO/*.short_summary.json  completeness scores and assembly counts

-- plus the genus/species descriptions and the curated gene symbols.

Roughly 45% of collections (nearly every qtl, gwas and maps one) publish no
CHECKSUM. Their file lists are predicted from the documented naming convention in
filetypes.yml, optionally HEAD-confirmed with ``verify=True``, and every file says
how it came to be known (SRC_*). Every collection likewise records how its file
list was obtained in ``index_status``: collapsing "we could not look" into "we
looked and there is nothing" would report streamable data as unreadable.

Datastore file names follow gensp.<collection_key>.<canonical_type>.<extensions>,
where collection_key is the collection directory name. Pairwise collections
(synteny, genome_alignments) instead use
gensp1.strain1.gnm.x.gensp2.strain2.gnm.<key>.<extensions>.
"""

import json
import logging
import os
import pathlib
import re
import subprocess
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

import yaml

DEFAULT_DATASTORE_URL = "https://data.legumeinfo.org"

# Vendored copy of the datastore file-content vocabulary; see filetypes.yml.
FILETYPES_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "filetypes.yml"
)

# How a file came to be known. Recorded per file so consumers (DSCensor, the
# legumista MCP server, JBrowse/BLAST config) can report provenance without
# re-deriving it -- they stay dumb readers of what this build decided.
SRC_CHECKSUM = "checksum"  # listed in CHECKSUM.md5; authoritative
SRC_VERIFIED = "verified"  # constructed by convention, HEAD-confirmed to exist
SRC_PREDICTED = "predicted"  # constructed by convention, not checked

# How a collection's file list was obtained.
STATUS_KNOWN = "known"  # from CHECKSUM
STATUS_VERIFIED = "verified"  # predicted, then HEAD-confirmed
STATUS_INFERRED = "inferred"  # predicted, not checked
STATUS_UNKNOWN = "unknown"  # no CHECKSUM and no documented convention

PROBE_TIMEOUT = int(os.environ.get("LIS_CATALOG_PROBE_TIMEOUT", "20"))
PROBE_WORKERS = int(os.environ.get("LIS_CATALOG_PROBE_WORKERS", "8"))

# Index siblings, mapped to the reader that consumes the file they belong to.
INDEX_SUFFIXES = (".fai", ".tbi", ".csi", ".bai", ".crai")
# Companions that are not themselves data and never determine access on their own.
COMPANION_SUFFIXES = (".gzi", ".md5")
# Matched with the dot: real data files are named e.g. "MANIFEST-000002".
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
# derived_from edge is the whole reason the index is a graph and not a table.
INHERITED_FIELDS = ("chromosome_prefix", "supercontig_prefix", "bioproject")

# <A>.x.<B>[.<epoch>].<KEY>.<ext>; see DatastoreIndex.pairwise_relationships.
PAIRWISE_PATTERN = re.compile(
    r"^(?P<a>[a-z]{4,6}\.[A-Za-z0-9_-]+\.gnm\d+)"
    r"\.x\."
    r"(?P<b>[a-z]{4,6}\.[A-Za-z0-9_-]+\.gnm\d+)"
    r"(?P<epoch>\.[A-Za-z0-9_]+)?"
    r"\.(?P<key>[A-Za-z0-9]{4})"
    r"\.(?P<ext>gff3\.gz|paf\.gz|bam)$"
)

EMPTY_VALUES = (None, "", [], {})


@dataclass
class DatastoreFiles:  # pylint: disable=too-many-instance-attributes
    """One file in a collection, with its position in the datastore.

    ``src`` says how the file came to be known; only files listed in a CHECKSUM
    carry an md5. ``indexes`` holds the index suffixes present alongside a data
    file, and ``index_unknown`` marks a predicted file whose siblings were never
    probed.
    """

    genus: str
    species: str
    collection_type: str
    collection_key: str
    relative_path: str
    md5: str = None
    src: str = SRC_CHECKSUM
    canonical_type: str = None
    extensions: str = ""
    gensp: str = None
    parents: list = field(default_factory=list)
    indexes: list = field(default_factory=list)
    index_unknown: bool = False
    description: str = None
    applications: list = field(default_factory=list)

    @property
    def filename(self):
        """Base name of the file, without any subdirectory."""
        return self.relative_path.rsplit("/", 1)[-1]

    @property
    def infraspecies(self):
        """Strain portion of the collection key."""
        return self.collection_key.split(".")[0]

    @property
    def is_metadata(self):
        """True for README, CHECKSUM, MANIFEST, CHANGES and LINKOUTS files."""
        return self.filename == "README" or self.filename.startswith(METADATA_PREFIXES)

    @property
    def is_companion(self):
        """True for index siblings (.fai, .tbi, ...) and companions (.gzi, .md5)."""
        return self.relative_path.endswith(INDEX_SUFFIXES + COMPANION_SUFFIXES)

    @property
    def is_data(self):
        """True for files that are neither metadata nor companions of another file."""
        return not (self.is_metadata or self.is_companion)

    @property
    def collection_path(self):
        """Path of the containing collection, relative to the datastore root."""
        return (
            f"{self.genus}/{self.species}/{self.collection_type}/{self.collection_key}"
        )

    def url(self, datastore_url):
        """Absolute URL for this file in the remote datastore."""
        return f"{datastore_url}/{self.collection_path}/{self.relative_path}"


@dataclass
class DatastoreCollection:  # pylint: disable=too-many-instance-attributes
    """One collection directory: its README metadata, its files and their provenance.

    ``metadata`` holds the non-empty README_FIELDS, plus any INHERITED_FIELDS taken
    from the genome it derives from (named in ``inherited``). ``files`` is every
    file, metadata and index siblings included, sorted by path.
    """

    genus: str
    species: str
    collection_type: str
    collection_key: str
    metadata: dict = field(default_factory=dict)
    index_status: str = STATUS_UNKNOWN
    files: list = field(default_factory=list)
    busco: dict = None
    counts: dict = None
    derived_from: list = field(default_factory=list)
    inherited: list = field(default_factory=list)

    @property
    def path(self):
        """Path of the collection, relative to the datastore root."""
        return (
            f"{self.genus}/{self.species}/{self.collection_type}/{self.collection_key}"
        )

    @property
    def data_files(self):
        """Files that are neither metadata nor index siblings/companions."""
        return [record for record in self.files if record.is_data]

    def url(self, datastore_url):
        """Absolute URL for this collection in the remote datastore."""
        return f"{datastore_url}/{self.path}"


def split_on_key(basename, collection_key):
    """Split a file name into canonical type and extensions.

    Tries the full collection key first, then just its trailing token, which is
    how pairwise alignment files embed the key. Returns (None, "") when neither
    appears, which marks an ancillary file such as a log or usage policy.
    """
    for marker in (collection_key, collection_key.split(".")[-1]):
        token = f".{marker}."
        if token in basename:
            rest = basename.split(token, 1)[1].split(".")
            return rest[0], ".".join(rest[1:])
    return None, ""


def split_pairwise_parents(basename):
    """Return the two assembly prefixes of a pairwise file, or an empty list.

    Names of the form gensp1.strain1.gnmN.x.gensp2.strain2.gnmM.key.ext encode
    both sides of an alignment around a literal '.x.' separator.
    """
    parts = basename.split(".")
    if "x" not in parts:
        return []
    index = parts.index("x")
    left = parts[:index]
    right = parts[index + 1 : index + 4]
    if len(left) < 3 or len(right) < 3:
        return []
    return [".".join(left[:3]), ".".join(right[:3])]


class DatastoreIndex:  # pylint: disable=too-many-instance-attributes,too-many-public-methods
    """Every collection and file in the datastore, from a datastore-metadata checkout."""

    STAT_KEYS = (
        "collections",
        "files",
        "indexed_files",
        "with_doi",
        "index_unknown",
        "with_busco",
        "curated_symbols",
        "described_taxa",
        "pairwise_files",
        "paired_genomes",
        "predicted_files",
        "verified_files",
        "probes",
    )

    def __init__(self, root, logger=None, datastore_url=None, verify=False):
        self.root = pathlib.Path(os.path.abspath(root))
        self.logger = logger or logging.getLogger(__name__)
        self.datastore_url = (datastore_url or DEFAULT_DATASTORE_URL).rstrip("/")
        # When True, every predicted filename is confirmed with a HEAD request before
        # it enters the index. Off by default so a build stays offline and fast. The
        # difference is visible in each file's `src`, never silent.
        self.verify = verify
        self.filetypes = self._load_filetypes()
        self.files = []  # every file record, across all collections
        self.collections = {}  # collection path -> DatastoreCollection, sorted by path
        self.taxa = {}
        self.gene_symbols = {}
        self.pairwise = []
        self.stats = dict.fromkeys(self.STAT_KEYS, 0)

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
        """The datastore-metadata commit this index was built from.

        Recorded so a consumer can tell how stale its copy is. A checkout sitting
        hundreds of commits behind under-reports what is streamable, and the only
        defence is making the provenance visible.
        """
        try:
            out = subprocess.run(
                ["git", "-C", str(self.root), "rev-parse", "HEAD"],
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

        A missing vocabulary degrades to no prediction; it does not fail the build.
        """
        try:
            with open(FILETYPES_PATH, encoding="utf-8") as handle:
                loaded = yaml.safe_load(handle)
        except (OSError, yaml.YAMLError) as err:
            self.logger.warning("could not read %s: %s", FILETYPES_PATH, err)
            return {}
        return loaded if isinstance(loaded, dict) else {}

    # ---------------------------------------------------------------- prediction
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
            stem = f"{abbrev}.{identifier}"
        else:
            return []  # abbrev-prefixed type with no abbrev: nothing safe to build
        return [f"{stem}.{content}{extension}" for content in contents]

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
        except Exception:  # pylint: disable=broad-exception-caught
            return False  # a 404 (or any failure) means "not there"

    def verify_files(self, base_url, names):
        """Keep the names that actually exist. Probes run concurrently."""
        if not names:
            return []
        urls = [f"{base_url.rstrip('/')}/{name}" for name in names]
        with ThreadPoolExecutor(max_workers=PROBE_WORKERS) as pool:
            found = list(pool.map(self.url_exists, urls))
        self.stats["probes"] += len(urls)
        return [name for name, ok in zip(names, found) if ok]

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

    def probe_indexes(self, collection, records):
        """Records for the index siblings of ``records`` that a HEAD request finds.

        Only reached under verify, and only for types that can carry indexes, so the
        cost is a few dozen requests. Without it those files would have to be reported
        as index-status-unknown.
        """
        names = [
            f"{record.relative_path}{suffix}"
            for record in records
            for suffix in self.plausible_indexes(record.relative_path)
        ]
        found = self.verify_files(collection.url(self.datastore_url), names)
        return [self.file_record(collection, name, src=SRC_VERIFIED) for name in found]

    def predicted_records(self, collection):
        """File records for a collection without a CHECKSUM, from the naming convention.

        Also sets the collection's index_status to say how far the prediction was
        checked.
        """
        ctype = collection.collection_type
        names = self.predicted_files(
            ctype,
            collection.collection_key,
            collection.metadata.get("scientific_name_abbrev"),
        )
        if not names:
            collection.index_status = STATUS_UNKNOWN
            return []
        if self.verify:
            names = self.verify_files(collection.url(self.datastore_url), names)
            src, collection.index_status = SRC_VERIFIED, STATUS_VERIFIED
            self.stats["verified_files"] += len(names)
        else:
            src, collection.index_status = SRC_PREDICTED, STATUS_INFERRED
            self.stats["predicted_files"] += len(names)
        records = [self.file_record(collection, name, src=src) for name in names]
        # Types documented as never carrying index siblings need no probing to know
        # none of these files is randomly accessible. For the types that CAN carry
        # them, an unprobed file's index status is genuinely unknown -- reporting it
        # as "not indexed" would mark streamable data unreadable, which is the same
        # mistake `index_status` exists to prevent. Verified: a markers .gff3.gz
        # constructed this way does have a .tbi.
        if records and self.type_is_indexed(ctype):
            if self.verify:
                records.extend(self.probe_indexes(collection, records))
            else:
                for record in records:
                    record.index_unknown = True
        return records

    # ------------------------------------------------------------------- parsing
    @staticmethod
    def file_record(collection, relative_path, md5=None, src=SRC_CHECKSUM):
        """A DatastoreFiles record in ``collection``, with its name parsed."""
        record = DatastoreFiles(
            genus=collection.genus,
            species=collection.species,
            collection_type=collection.collection_type,
            collection_key=collection.collection_key,
            relative_path=relative_path,
            md5=md5,
            src=src,
        )
        if not record.is_metadata:
            basename = record.filename
            record.canonical_type, record.extensions = split_on_key(
                basename, collection.collection_key
            )
            record.gensp = basename.split(".")[0]
            record.parents = split_pairwise_parents(basename)
        return record

    def parse_checksum(self, text, collection):
        """Turn one manifest's contents into DatastoreFiles records.

        Only a leading ``./`` is stripped from each name. Stripping every leading dot
        and slash would turn ``./.nextflow.log`` into ``nextflow.log``, a file that
        does not exist.
        """
        records = []
        for line in text.splitlines():
            fields = line.split(None, 1)  # the name is everything after the md5
            if len(fields) < 2:
                continue
            md5, relative_path = fields[0], fields[1].strip()
            if relative_path.startswith("./"):
                relative_path = relative_path[2:]
            records.append(self.file_record(collection, relative_path, md5=md5))
        return records

    def checksum_files(self, coll_dir, collection):
        """Records for every file the collection's CHECKSUM lists, or None without one.

        ``None`` and ``[]`` mean different things here and callers rely on it:
        None is "we could not look", [] is "we looked and it is empty".

        Every ``CHECKSUM.*.md5`` in the directory is read, not just the one named for
        the directory: ``legume.fam1.M65K`` publishes ``CHECKSUM.mixed.fam1.M65K.md5``,
        and one diversity collection publishes two. A file listed twice is kept once.
        """
        paths = sorted(pathlib.Path(coll_dir).glob("CHECKSUM.*.md5"))
        if not paths:
            return None
        records = {}
        for path in paths:
            for record in self.parse_checksum(self._read_text(path), collection):
                records.setdefault(record.relative_path, record)
        return list(records.values())

    def manifest_descriptions(self, coll_dir, identifier):
        """Map file name -> {description, applications} from MANIFEST."""
        path = os.path.join(coll_dir, f"MANIFEST.{identifier}.yml")
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
        need no ``.fai`` fetch and the build stays fully offline.
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
        for key, summary_field in (
            ("scaffolds", "Number of scaffolds"),
            ("contigs", "Number of contigs"),
            ("total_length", "Total length"),
            ("scaffold_n50", "Scaffold N50"),
            ("contig_n50", "Contigs N50"),
        ):
            raw = results.get(summary_field)
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

    # ------------------------------------------------------------------ building
    def find_collections(self):
        """Yield (genus, species, collection_type, collection_key, dirpath, filenames).

        A collection qualifies on *any* metadata file. Keying off CHECKSUM alone
        misses ~45% of the store (nearly every qtl, gwas and maps collection); keying
        off README alone drops the five genome_alignments collections that publish
        only a CHECKSUM, along with the indexed BAMs they hold.
        """
        root = str(self.root)
        for dirpath, _dirnames, filenames in os.walk(root):
            if os.sep + ".git" in dirpath:
                continue
            relative = os.path.relpath(dirpath, root).split(os.sep)
            if len(relative) != 4:  # Genus/species/type/collection
                continue
            has_metadata = self._readme_filename(filenames) or any(
                name.startswith(("CHECKSUM.", "MANIFEST.")) for name in filenames
            )
            if has_metadata:
                yield (*relative, dirpath, filenames)

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

    def read_readme(self, dirpath, filenames):
        """The collection's parsed README, or {} when it is missing or unreadable."""
        readme_file = self._readme_filename(filenames)
        if not readme_file:
            return {}
        readme = self._read_yaml(os.path.join(dirpath, readme_file))
        if not readme:
            self.logger.warning("unreadable README in %s", dirpath)
        return readme

    def build_collection(self, coll_dir, collection, readme):
        """Fill in one collection's metadata, files and metrics. Returns it."""
        for name in README_FIELDS:
            value = readme.get(name)
            if value not in EMPTY_VALUES:
                collection.metadata[name] = value

        records = self.checksum_files(coll_dir, collection)
        if records is not None:
            collection.index_status = STATUS_KNOWN
        else:
            records = self.predicted_records(collection)

        present = {record.relative_path for record in records}
        manifest = self.manifest_descriptions(coll_dir, collection.collection_key)
        for record in records:
            if record.is_data:
                record.indexes = [
                    suffix
                    for suffix in INDEX_SUFFIXES
                    if record.relative_path + suffix in present
                ]
            extra = manifest.get(record.filename, {})
            record.description = extra.get("description")
            record.applications = extra.get("applications", [])
        collection.files = sorted(records, key=lambda record: record.relative_path)
        collection.busco, collection.counts = self.busco_summary(coll_dir)
        return collection

    @staticmethod
    def link_lineage(collections):
        """Add ``derived_from`` edges and inherit assembly conventions along them.

        An annotation identifier embeds the assembly it was called against
        (``Wm82.gnm4.ann1.T8TQ`` -> ``Wm82.gnm4``), so the edge is derivable
        without any extra metadata. ``chromosome_prefix`` and friends live on the
        genome README only, which is what makes the inheritance necessary.
        """
        genomes = {}
        for collection in collections:
            if collection.collection_type != "genomes":
                continue
            stem = ".".join(collection.collection_key.split(".")[:2])  # strain.gnmN
            genomes.setdefault((collection.genus, collection.species, stem), collection)

        for collection in collections:
            parts = collection.collection_key.split(".")
            if len(parts) < 3 or not parts[1].startswith("gnm"):
                continue
            stem = ".".join(parts[:2])
            parent = genomes.get((collection.genus, collection.species, stem))
            if parent is None or parent is collection:
                continue
            collection.derived_from = [parent.collection_key]
            for name in INHERITED_FIELDS:
                if name not in collection.metadata and name in parent.metadata:
                    collection.metadata[name] = parent.metadata[name]
                    collection.inherited.append(name)

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
        pairs = []
        for collection in self.collections.values():
            if collection.collection_type not in ("synteny", "genome_alignments"):
                continue
            for record in collection.data_files:
                match = PAIRWISE_PATTERN.match(record.relative_path)
                if not match:
                    continue
                epoch = (match.group("epoch") or "").lstrip(".")
                pair = {
                    "a": match.group("a"),
                    "b": match.group("b"),
                    "kind": (
                        "synteny"
                        if collection.collection_type == "synteny"
                        else "alignment"
                    ),
                    "format": match.group("ext"),
                    "collection": collection.path,
                    "file": record.relative_path,
                    "url": record.url(self.datastore_url),
                }
                if epoch:
                    pair["epoch"] = epoch  # whole-genome-duplication label
                if match.group("a") == match.group("b"):
                    pair["self"] = True
                if record.indexes:
                    pair["i"] = list(record.indexes)
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
        for dirpath, _dirnames, filenames in os.walk(str(self.root)):
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
                for doc_field in (
                    "taxid",
                    "abbrev",
                    "commonName",
                    "description",
                    "genus",
                ):
                    value = doc.get(doc_field)
                    if value not in EMPTY_VALUES:
                        record[doc_field] = value
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
        for dirpath, _dirnames, filenames in os.walk(str(self.root)):
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

    def count(self):
        """Fill in the build diagnostics that summarise the finished index."""
        stats = self.stats
        stats["collections"] = len(self.collections)
        for collection in self.collections.values():
            data = collection.data_files
            stats["files"] += len(data)
            stats["indexed_files"] += sum(1 for record in data if record.indexes)
            stats["with_doi"] += bool(collection.metadata.get("publication_doi"))
            stats["index_unknown"] += collection.index_status == STATUS_UNKNOWN
            # Counted from what was parsed, so it stays truthful while the catalog
            # withholds the metrics themselves (see the NOTE in catalog.py).
            stats["with_busco"] += collection.busco is not None
        stats["curated_symbols"] = sum(len(v) for v in self.gene_symbols.values())
        stats["described_taxa"] = len(self.taxa)
        stats["pairwise_files"] = len(self.pairwise)
        stats["paired_genomes"] = len(
            {genome for pair in self.pairwise for genome in (pair["a"], pair["b"])}
        )

    def build(self):
        """Read the whole checkout and populate the index. Returns self."""
        self.stats = dict.fromkeys(self.STAT_KEYS, 0)
        found = []
        for genus, species, ctype, key, dirpath, filenames in self.find_collections():
            collection = DatastoreCollection(genus, species, ctype, key)
            readme = self.read_readme(dirpath, filenames)
            found.append(self.build_collection(dirpath, collection, readme))
        self.link_lineage(found)
        found.sort(key=lambda collection: collection.path)
        self.collections = {collection.path: collection for collection in found}
        self.files = [record for collection in found for record in collection.files]
        self.pairwise = self.pairwise_relationships()
        self.taxa = self.descriptions()
        self.gene_symbols = self.curated_symbols()
        self.count()
        self.logger.info(
            "Indexed %s files across %s collections",
            len(self.files),
            len(self.collections),
        )
        return self

    # ------------------------------------------------------------------- queries
    def select(self, collection_type=None, canonical_type=None, endswith=None):
        """Filter indexed files by collection type, canonical type, or file suffix.

        Predicted files are included; check ``src`` when only CHECKSUM-listed or
        HEAD-confirmed files will do.
        """
        results = []
        for record in self.files:
            if collection_type and record.collection_type != collection_type:
                continue
            if canonical_type and record.canonical_type != canonical_type:
                continue
            if endswith and not record.relative_path.endswith(endswith):
                continue
            results.append(record)
        return results

    def _collection_files(self, record):
        collection = self.collections.get(record.collection_path)
        return collection.files if collection else []

    def sibling(self, record, canonical_type):
        """Find a file of another canonical type in the same collection."""
        for other in self._collection_files(record):
            if other.canonical_type == canonical_type:
                return other
        return None

    def companion(self, record, suffix):
        """Find an index or companion file, such as a .fai or .bai, for a record."""
        target = f"{record.relative_path}{suffix}"
        for other in self._collection_files(record):
            if other.relative_path == target:
                return other
        return None

    def orphan_collections(self):
        """Collections that publish no CHECKSUM, so their file lists are not authoritative."""
        return [
            path
            for path, collection in self.collections.items()
            if collection.index_status != STATUS_KNOWN
        ]
