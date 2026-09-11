"""Serialize the datastore index as a single catalog document.

``DatastoreIndex`` (datastore_files.py) does the reading: every collection in a
datastore-metadata checkout, its README metadata, its files and how each came to be
known, with no network access. This module flattens that into the one JSON document
that downstream consumers (DSCensor, the legumista MCP server, any future client) can
hold in memory. Its shape is a published contract: change it only together with
SCHEMA_VERSION.

Why a separate module from ``process_collections``: that class scrapes the store
to drive JBrowse/BLAST/Jekyll and needs the network to do it. The catalog needs
none, so welding the two together would make a sub-second job inherit a
multi-minute one. Keeping them apart lets the catalog rebuild on every metadata
push while enrichment runs on its own cadence.

File entries use short keys because the document ships whole: ``n`` is the path
within the collection, ``i`` the index suffixes present beside it, ``src`` how the
file came to be known, and ``i_unknown`` marks a predicted file whose index siblings
were never probed. Collections that publish no CHECKSUM get an ``index_status``
other than ``"known"`` -- never an empty list passed off as authoritative.
"""

import json
import logging
import os
from datetime import datetime, timezone

from .datastore_files import DatastoreIndex

SCHEMA_VERSION = 1


class CatalogBuilder:
    """Build the catalog document from a datastore-metadata checkout."""

    def __init__(self, metadata_dir, logger=None, datastore_url=None, verify=False):
        self.logger = logger or logging.getLogger("catalog")
        self.index = DatastoreIndex(
            metadata_dir,
            logger=self.logger,
            datastore_url=datastore_url,
            verify=verify,
        )

    @property
    def stats(self):
        """Build diagnostics from the most recent build."""
        return self.index.stats

    @staticmethod
    def file_entry(record):
        """The catalog entry for one data file."""
        entry = {"n": record.relative_path, "src": record.src}
        if record.indexes:
            entry["i"] = list(record.indexes)
        if record.index_unknown:
            entry["i_unknown"] = True
        if record.description:
            entry["description"] = record.description
        if record.applications:
            entry["applications"] = list(record.applications)
        return entry

    def collection_record(self, collection):
        """The catalog record for one collection."""
        record = {
            "path": collection.path,
            "id": collection.collection_key,
            "type": collection.collection_type,
            "genus": collection.genus,
            "species": collection.species,
            "base_url": collection.url(self.index.datastore_url),
            "index_status": collection.index_status,
            "files": [self.file_entry(f) for f in collection.data_files],
        }
        record.update(collection.metadata)
        if collection.derived_from:
            record["derived_from"] = list(collection.derived_from)
        if collection.inherited:
            record["inherited"] = list(collection.inherited)
        # NOTE: assembly metrics are parsed but deliberately NOT emitted yet.
        # The index still ingests both the completeness scores and the assembly counts
        # from the committed BUSCO summary (collection.busco / collection.counts) --
        # that code is live and tested, and needs no network. Only the inclusion in
        # the catalog is deferred, so re-enabling it is uncommenting these four lines.
        # Anything downstream that consumes `busco`/`counts` is commented out with a
        # matching NOTE.
        # if collection.busco:
        #     record["busco"] = collection.busco
        # if collection.counts:
        #     record["counts"] = collection.counts
        return record

    def build(self):
        """Build the whole catalog. Returns the document as a dict."""
        index = self.index.build()
        return {
            "schema": SCHEMA_VERSION,
            "built_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "source_commit": index.source_commit(),
            "datastore_url": index.datastore_url,
            "stats": dict(index.stats),
            "taxa": index.taxa,
            "pairwise": index.pairwise,
            "gene_symbols": index.gene_symbols,
            "collections": [
                self.collection_record(c) for c in index.collections.values()
            ],
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
