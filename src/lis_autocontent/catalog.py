"""Serialize a DatastoreIndex as one catalog document.

The shape is a published contract. Consumers refuse a schema they don't know, so
an optional key may be added under the same SCHEMA_VERSION; anything else changes
it. File entries
use short keys: ``n`` path, ``i`` index suffixes, ``src`` provenance, ``i_unknown``
siblings never probed.
"""

import json
import logging
import os
from datetime import datetime, timezone

from . import gcv, jbrowse
from .datastore_files import DatastoreIndex

SCHEMA_VERSION = 1


class CatalogBuilder:
    """Build the catalog document from a datastore-metadata checkout."""

    def __init__(  # pylint: disable=too-many-arguments,too-many-positional-arguments
        self,
        metadata_dir,
        logger=None,
        datastore_url=None,
        verify=False,
        jbrowse_config=None,
        gcv_config=None,
    ):
        self.logger = logger or logging.getLogger("catalog")
        # all-genera's deployed config.json; without it JBrowse is left out entirely.
        self.jbrowse_config = jbrowse_config
        # A gcv.yml; without it GCV is left out entirely.
        self.gcv_config = gcv_config
        self.deployments = []
        self.placements = {}
        self.gcv_placements = {}
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
        if collection.path in self.placements:
            record["jbrowse"] = self.placements[collection.path]
        if collection.path in self.gcv_placements:
            record["gcv"] = self.gcv_placements[collection.path]
        # NOTE: BUSCO metrics are parsed but deliberately withheld; uncomment to emit.
        # Downstream consumers of `busco`/`counts` carry a matching NOTE.
        # if collection.busco:
        #     record["busco"] = collection.busco
        # if collection.counts:
        #     record["counts"] = collection.counts
        return record

    def build(self):
        """Build the whole catalog. Returns the document as a dict."""
        index = self.index.build()
        instances = None
        if self.jbrowse_config:
            self.deployments = [
                jbrowse.read_deployment(
                    index, jbrowse.INSTANCE, jbrowse.INSTANCE_URL, self.jbrowse_config
                )
            ]
            instances, self.placements = jbrowse.catalog_section(self.deployments)
        gcv_instances = self.gcv_section(index)
        document = {
            "schema": SCHEMA_VERSION,
            "built_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "source_commit": index.source_commit(),
            "datastore_url": index.datastore_url,
            "stats": dict(index.stats),
            "taxa": index.taxa,
            "pairwise": index.pairwise,
            "gene_symbols": index.gene_symbols,
        }
        if instances is not None:
            document["jbrowse_instances"] = instances
        if gcv_instances is not None:
            document["gcv_instances"] = gcv_instances
        document["collections"] = [
            self.collection_record(c) for c in index.collections.values()
        ]
        return document

    def gcv_section(self, index):
        """gcv_instances, setting the placements; None, leaving GCV out, without one."""
        if not self.gcv_config:
            return None
        try:
            config = gcv.load_config(self.gcv_config)
            instances, self.gcv_placements = gcv.catalog_section(
                index, config, self.logger
            )
        except gcv.GcvError as err:
            self.logger.warning("leaving GCV out of the catalog: %s", err)
            return None
        return instances

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
