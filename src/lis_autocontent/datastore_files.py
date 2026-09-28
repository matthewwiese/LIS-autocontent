"""Index every collection and file in the datastore from a datastore-metadata checkout.

Offline. CHECKSUM.<key>.md5 is the authoritative file list and the only listing of
.fai/.tbi siblings. Collections without one get files predicted from filetypes.yml,
labelled by ``src`` and ``index_status``.
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

# How a file came to be known, recorded so consumers needn't re-derive it.
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

# Siblings that make a data file randomly accessible.
INDEX_SUFFIXES = (".fai", ".tbi", ".csi", ".bai", ".crai")
# Companions that are not themselves data and never determine access on their own.
COMPANION_SUFFIXES = (".gzi", ".md5")
# Matched with the dot: real data files are named e.g. "MANIFEST-000002".
METADATA_PREFIXES = ("README.", "MANIFEST.", "CHANGES.", "CHECKSUM.", "LINKOUTS.")

# README fields copied verbatim onto a collection.
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

# README_FIELDS carried on each DSCensor node; long prose stays in the catalog.
# Shared with process_collections so the two paths can't drift.
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

# Genome README fields an annotation inherits along derived_from.
INHERITED_FIELDS = ("chromosome_prefix", "supercontig_prefix", "bioproject")

# <A>.x.<B>[.<epoch>].<KEY>[.<program>].<ext>; see DatastoreIndex.pairwise_relationships.
PAIRWISE_PATTERN = re.compile(
    r"^(?P<a>[a-z]{4,6}\.[A-Za-z0-9_-]+\.gnm\d+)"
    r"\.x\."
    r"(?P<b>[a-z]{4,6}\.[A-Za-z0-9_-]+\.gnm\d+)"
    r"(?P<epoch>\.[A-Za-z0-9_]+)?"
    r"\.(?P<key>[A-Za-z0-9]{4})"
    r"(?:\.[A-Za-z0-9_-]+)?"  # program, e.g. Cicer's "...PXV3.minimap2.paf.gz"
    r"\.(?P<ext>gff3\.gz|paf\.gz|bam)$"
)

EMPTY_VALUES = (None, "", [], {})


@dataclass
class DatastoreFiles:  # pylint: disable=too-many-instance-attributes
    """One file in a collection.

    ``src`` says how it came to be known; only CHECKSUM files carry an md5.
    ``index_unknown`` marks a predicted file whose siblings were never probed.
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
    """One collection directory: its README metadata, files and provenance.

    ``files`` includes metadata and index siblings, sorted by path. ``metadata``
    includes INHERITED_FIELDS from its genome, named in ``inherited``.
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
    """(canonical_type, extensions) around the collection key, or (None, "") if absent.

    Also matches the key's last token, which is how pairwise files embed it.
    """
    for marker in (collection_key, collection_key.split(".")[-1]):
        token = f".{marker}."
        if token in basename:
            rest = basename.split(token, 1)[1].split(".")
            return rest[0], ".".join(rest[1:])
    return None, ""


def split_pairwise_parents(basename):
    """The two assembly prefixes of a pairwise ``A.x.B`` file name, or []."""
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
        # HEAD-confirm predicted files; off by default so a build stays offline.
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
        """A metadata file's text, or "" when it is absent or unreadable."""
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
        """The datastore-metadata commit this index was built from, or None."""
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
        """The vendored content vocabulary, or {} if it can't be read."""
        try:
            with open(FILETYPES_PATH, encoding="utf-8") as handle:
                loaded = yaml.safe_load(handle)
        except (OSError, yaml.YAMLError) as err:
            self.logger.warning("could not read %s: %s", FILETYPES_PATH, err)
            return {}
        return loaded if isinstance(loaded, dict) else {}

    # ---------------------------------------------------------------- prediction
    def predicted_files(self, ctype, identifier, abbrev):
        """Filenames a collection of this type is expected to publish, or [].

        A prediction, not a guarantee: callers must label the result (SRC_*).
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
        """HEAD one datastore URL, percent-encoded since urllib rejects non-ASCII."""
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
        """Index suffixes possible for this file's type; others are certain 404s."""
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
        """Records for the index siblings of ``records`` that a HEAD request finds."""
        names = [
            f"{record.relative_path}{suffix}"
            for record in records
            for suffix in self.plausible_indexes(record.relative_path)
        ]
        found = self.verify_files(collection.url(self.datastore_url), names)
        return [self.file_record(collection, name, src=SRC_VERIFIED) for name in found]

    def predicted_records(self, collection):
        """Convention-predicted file records; sets the collection's index_status."""
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
        # Unprobed siblings are unknown, not absent.
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
        """DatastoreFiles records from one CHECKSUM's text.

        Strips only a leading ``./``, so ``./.nextflow.log`` stays ``.nextflow.log``.
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
        """Records for every file listed by the directory's CHECKSUM.*.md5 files.

        None means there is no CHECKSUM and [] that it lists nothing; callers depend on
        the difference. A file listed twice is kept once.
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
        """(busco, counts) from the committed BUSCO short_summary, or (None, None)."""
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
        """Yield (genus, species, type, key, dirpath, filenames) for each collection.

        A collection is any Genus/species/type/collection directory holding a README,
        CHECKSUM or MANIFEST; no one of them is present in every collection.
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
        """The collection's README file name, accepting a bare ``README``, or None."""
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
        """Add ``derived_from`` edges and copy INHERITED_FIELDS down them.

        The parent is in the identifier: ``Wm82.gnm4.ann1.T8TQ`` -> ``Wm82.gnm4``.
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
        """Genome pairs parsed from synteny and alignment file names.

        Use this over a genome's own collection: each file is stored once, under its
        reference genome. Self-comparisons carry a duplication ``epoch``.
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
        """Genus- and species-level taxon metadata from ``about_this_collection``."""
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
        """Curated gene symbol -> gene, from ``gene_functions/<abbrev>.traits.yml``.

        Walked separately, as these directories hold no collection metadata.
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
            # Counts parsed summaries, whether or not the catalog emits them.
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
        """Filter indexed files by collection type, canonical type or suffix.

        Includes predicted files; filter on ``src`` for listed or confirmed ones only.
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
