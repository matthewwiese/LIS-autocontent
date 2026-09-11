"""CLI for interacting with the ProcessCollections class."""

#!/usr/bin/env python3

import sys
import logging
import click
from .catalog import CatalogBuilder
from .process_collections import ProcessCollections


def setup_logging(log_file, log_level, process):
    """initializes a logger object with a common format"""
    log_level = getattr(
        logging, log_level.upper(), logging.INFO
    )  # set provided or set INFO
    msg_format = "%(asctime)s|%(name)s|[%(levelname)s]: %(message)s"
    logging.basicConfig(format=msg_format, datefmt="%m-%d %H:%M", level=log_level)
    log_handler = logging.FileHandler(log_file, mode="w")
    formatter = logging.Formatter(msg_format)
    log_handler.setFormatter(formatter)
    logger = logging.getLogger(
        f"{process}"
    )  # sets what will be printed for the log process
    logger.addHandler(log_handler)
    return logger


@click.command()
@click.option(
    "--taxa_list",
    default="../_data/taxon_list.yml",
    help="""Taxa.yml file. (Default: ../_data/taxon_list.yml)""",
)
@click.option(
    "--collections_out", default="../_data/taxa/", help="""Output for collections."""
)
@click.option(
    "--from_github",
    default="./datastore-metadata",
    help="""Path to datastore-metadata github directory. (Default: ./datastore-metadata).""",
)
@click.option(
    "--log_file",
    default="./populate-jekyll.log",
    help="""Log file to output messages. (default: ./populate-jekyll.log)""",
)
@click.option(
    "--log_level",
    default="INFO",
    help="""Log Level to output messages. (default: INFO)""",
)
def populate_jekyll(taxa_list, collections_out, from_github, log_file, log_level):
    """CLI entry for populate-jekyll"""
    logger = setup_logging(log_file, log_level, "populate-jekyll")
    logger.info("Processing Collections...")
    parser = ProcessCollections(logger, out_dir=collections_out)  # initialize class
    logger.info("Outputting Collections...")
    parser.parse_collections(taxa_list, from_github)  # parse_collections


@click.command()
@click.option(
    "--nodes_out",
    default="./autocontent",
    help="""Output for dscensor nodes.""",
)
@click.option(
    "--from_github",
    default="./datastore-metadata",
    help="""Path to datastore-metadata github directory. (Default: ./datastore-metadata).""",
)
@click.option(
    "--log_file",
    default="./populate-dscensor.log",
    help="""Log file to output messages. (default: ./populate-dscensor.log)""",
)
@click.option(
    "--log_level",
    default="INFO",
    help="""Log Level to output messages. (default: INFO)""",
)
def populate_dscensor(nodes_out, from_github, log_file, log_level):
    """CLI entry for populate-dscensor"""
    logger = setup_logging(log_file, log_level, "populate-dscensor")
    parser = ProcessCollections(logger, out_dir=nodes_out)
    logger.info("Processing Collections...")
    parser.parse_collections(from_github)
    logger.info("Creating DSCensor Nodes...")
    parser.populate_dscensor(nodes_out)


@click.command()
@click.option("--jbrowse_url", help="""URL hosting JBrowse2""")
@click.option(
    "--taxa_list",
    default="../_data/taxon_list.yml",
    help="""Taxa.yml file. (Default: ../_data/taxon_list.yml)""",
)
@click.option(
    "--datastore_url",
    default="https://data.legumeinfo.org",
    help="""URL hosting datastore formatted files.""",
)
@click.option(
    "--jbrowse_out",
    default="./autocontent",
    help="""Output directory for Jbrowse2. (Default: ./autocontent)""",
)
@click.option(
    "--from_github",
    default="./datastore-metadata",
    help="""Path to datastore-metadata github directory. (Default: ./datastore-metadata).""",
)
@click.option(
    "--cmds_only",
    is_flag=True,
    help="""Output commands only. Do not run Jbrowse2 just output the commands that would be run.""",
)
@click.option(
    "--log_file",
    default="./populate-jbrowse2.log",
    help="""Log file to output messages. (default: ./populate-jbrowse2.log)""",
)
@click.option(
    "--log_level",
    default="INFO",
    help="""Log Level to output messages. (default: INFO)""",
)
def populate_jbrowse2(
    jbrowse_url,
    datastore_url,
    taxa_list,
    jbrowse_out,
    from_github,
    cmds_only,
    log_file,
    log_level,
):
    """CLI entry for populate-jbrowse2"""
    logger = setup_logging(log_file, log_level, "populate-jbrowse2")
    if not jbrowse_url:
        logger.error("--jbrowse_url required for populate-jbrowse2")
        sys.exit(1)
    parser = ProcessCollections(
        logger,
        jbrowse_url=jbrowse_url,
        datastore_url=datastore_url,
        out_dir=jbrowse_out,
    )  # initialize class
    logger.info("Processing Collections...")
    parser.parse_collections(taxa_list, from_github)  # parse_collections
    logger.info("Creating JBrowse2 Config...")
    parser.populate_jbrowse2(jbrowse_out, cmds_only)  # populate JBrowse2


@click.command()
@click.option(
    "--taxa_list",
    default="../_data/taxon_list.yml",
    help="""Taxa.yml file. (Default: ../_data/taxon_list.yml)""",
)
@click.option(
    "--blast_out",
    default="./autocontent",
    help="""Output directory for BLAST DBs. (Default: ./autocontent)""",
)
@click.option(
    "--from_github",
    default="./datastore-metadata",
    help="""Path to datastore-metadata github directory. (Default: ./datastore-metadata).""",
)
@click.option(
    "--cmds_only",
    is_flag=True,
    help="""Output commands only. Do not run makeblastdb just output the commands that would be run.""",
)
@click.option(
    "--log_file",
    default="./populate-blast.log",
    help="""Log file to output messages. (default: ./populate-blast.log)""",
)
@click.option(
    "--log_level",
    default="INFO",
    help="""Log Level to output messages. (default: INFO)""",
)
def populate_blast(taxa_list, blast_out, from_github, cmds_only, log_file, log_level):
    """CLI entry for populate-blast"""
    logger = setup_logging(log_file, log_level, "populate-blast")
    parser = ProcessCollections(logger, out_dir=blast_out)  # initialize class
    logger.info(f"Processing Collections from {taxa_list}")
    parser.parse_collections(taxa_list, from_github)  # parse_collections
    logger.info("Creating BLAST DBs...")
    parser.populate_blast(blast_out, cmds_only)  # populate BLAST


@click.command()
@click.option(
    "--from_github",
    default="./datastore-metadata",
    help="""Path to datastore-metadata github directory. (Default: ./datastore-metadata).""",
)
@click.option(
    "--catalog_out",
    default="./autocontent/catalog.json",
    help="""Output path for the catalog document. (Default: ./autocontent/catalog.json)""",
)
@click.option(
    "--datastore_url",
    default="https://data.legumeinfo.org",
    help="""URL the catalog should build file links against.""",
)
@click.option(
    "--indent",
    default=0,
    help="""JSON indent. 0 (default) writes the compact form intended for shipping.""",
)
@click.option(
    "--verify",
    is_flag=True,
    help="""Confirm every convention-derived filename with a HEAD request before it
    enters the catalog. Off by default so the build stays offline and takes ~2s; on, it
    costs ~1,650 requests and ~90s but marks those files 'verified' rather than
    'predicted'. Collections that publish a CHECKSUM are unaffected either way.""",
)
@click.option(
    "--log_file",
    default="./populate-catalog.log",
    help="""Log file to output messages. (default: ./populate-catalog.log)""",
)
@click.option(
    "--log_level",
    default="INFO",
    help="""Log Level to output messages. (default: INFO)""",
)
def populate_catalog(
    from_github, catalog_out, datastore_url, indent, verify, log_file, log_level
):
    """CLI entry for populate-catalog

    Builds the whole-store catalog from a datastore-metadata checkout. Unlike the
    other subcommands this needs no network and no taxa list: it describes every
    collection the checkout knows about, in one document.
    """
    logger = setup_logging(log_file, log_level, "populate-catalog")
    builder = CatalogBuilder(
        from_github, logger=logger, datastore_url=datastore_url, verify=verify
    )
    logger.info(
        "Building catalog from %s (%s)...",
        from_github,
        (
            "verifying predicted files"
            if verify
            else "offline; predicted files unverified"
        ),
    )
    builder.write(catalog_out, indent=indent or None)
