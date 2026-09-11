
"""Index every datastore file by reading CHECKSUM manifests rather than guessing URLs.

Each collection directory carries a CHECKSUM.<key>.md5 listing every file it
contains, including files in subdirectories such as BUSCO/. Treating that
manifest as the authoritative file list removes the need to construct URLs from
hardcoded naming conventions, and surfaces files the conventions never covered.

Datastore file names follow gensp.<collection_key>.<canonical_type>.<extensions>,
where collection_key is the collection directory name. Pairwise collections
(synteny, genome_alignments) instead use
gensp1.strain1.gnm.x.gensp2.strain2.gnm.<key>.<extensions>.
"""

import pathlib
import pathlib
from dataclasses import dataclass, field


METADATA_PREFIXES = ("README", "CHECKSUM", "MANIFEST", "CHANGES")


@dataclass
class DatastoreFiles:
    """One entry from a CHECKSUM manifest, with its position in the datastore."""

    genus: str
    species: str
    collection_type: str
    collection_key: str
    relative_path: str
    md5: str
    canonical_type: str = None
    extensions: str = ""
    gensp: str = None
    parents: list = field(default_factory=list)

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
        """True for README, CHECKSUM, MANIFEST and CHANGES files."""
        return self.filename.startswith(METADATA_PREFIXES)

    @property
    def collection_path(self):
        """Path of the containing collection, relative to the datastore root."""
        return (
            f"{self.genus}/{self.species}/{self.collection_type}/{self.collection_key}"
        )

    def url(self, datastore_url):
        """Absolute URL for this file in the remote datastore."""
        return f"{datastore_url}/{self.collection_path}/{self.relative_path}"


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


class DatastoreIndex:
    """Every file in the datastore, indexed from CHECKSUM manifests."""

    def __init__(self, root, logger=None):
        self.root = pathlib.Path(root)
        self.logger = logger
        self.files = []
        self.collections = {}

    def log(self, level, message):
        """Emit through the caller's logger when one was supplied."""
        if self.logger:
            getattr(self.logger, level)(message)

    def find_checksums(self):
        """Yield (genus, species, collection_type, collection_key, path) for every manifest."""
        for path in sorted(self.root.glob("*/*/*/*/CHECKSUM.*.md5")):
            genus, species, collection_type, collection_key = path.parts[-5:-1]
            yield genus, species, collection_type, collection_key, path

    def parse_checksum(self, text, genus, species, collection_type, collection_key):
        """Turn one manifest's contents into DatastoreFiles records."""
        records = []
        for line in text.splitlines():
            fields = line.split()
            if len(fields) < 2:
                continue
            md5, raw = fields[0], fields[1]
            relative_path = raw[2:] if raw.startswith("./") else raw
            record = DatastoreFiles(
                genus=genus,
                species=species,
                collection_type=collection_type,
                collection_key=collection_key,
                relative_path=relative_path,
                md5=md5,
            )
            if not record.is_metadata:
                basename = record.filename
                record.canonical_type, record.extensions = split_on_key(
                    basename, collection_key
                )
                record.gensp = basename.split(".")[0]
                record.parents = split_pairwise_parents(basename)
            records.append(record)
        return records

    def build(self):
        """Read every manifest under the root and populate the index."""
        for genus, species, ctype, key, path in self.find_checksums():
            text = path.read_text(encoding="utf-8", errors="replace")
            records = self.parse_checksum(text, genus, species, ctype, key)
            self.files.extend(records)
            self.collections[f"{genus}/{species}/{ctype}/{key}"] = records
        self.log(
            "info",
            f"Indexed {len(self.files)} files across {len(self.collections)} collections",
        )
        return self

    def select(self, collection_type=None, canonical_type=None, endswith=None):
        """Filter indexed files by collection type, canonical type, or file suffix."""
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

    def sibling(self, record, canonical_type):
        """Find a file of another canonical type in the same collection."""
        for other in self.collections.get(record.collection_path, []):
            if other.canonical_type == canonical_type:
                return other
        return None

    def companion(self, record, suffix):
        """Find an index or companion file, such as a .fai or .bai, for a record."""
        target = f"{record.relative_path}{suffix}"
        for other in self.collections.get(record.collection_path, []):
            if other.relative_path == target:
                return other
        return None

    def orphan_collections(self):
        """Collection directories that hold a README but no CHECKSUM manifest."""
        indexed = set(self.collections)
        orphans = []
        for readme in self.root.glob("*/*/*/*/README.*.yml"):
            key = "/".join(readme.parts[-5:-1])
            if key not in indexed:
                orphans.append(key)
        return sorted(orphans)
