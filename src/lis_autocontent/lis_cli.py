"""CLI for interacting with the ProcessCollections class."""

#!/usr/bin/env python3

import os
import subprocess
import sys
import logging
import click
from .catalog import CatalogBuilder
from .datastore_files import DatastoreIndex
from .divbrowse import DivbrowseError, compose_file, load_config
from . import jbrowse
from .jbrowse import JBrowseError, jekyll_url
from .process_collections import ProcessCollections


def parse_pairs(form):
    """A click callback turning repeated KEY=VALUE options into a dict."""

    def parse(_ctx, _param, values):
        pairs = {}
        for value in values:
            key, sep, item = value.partition("=")
            if not (sep and key and item):
                raise click.BadParameter(f"{value!r} is not {form}")
            pairs[key] = item
        return pairs

    return parse


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
    "--jbrowse_url",
    default=None,
    help="""Base URL of the JBrowse 2 instance the resource links open.
    (Default: jekyll_instance in jbrowse.yml)""",
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
def populate_jekyll(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    taxa_list, collections_out, jbrowse_url, from_github, log_file, log_level
):
    """CLI entry for populate-jekyll"""
    logger = setup_logging(log_file, log_level, "populate-jekyll")
    try:
        jbrowse_url = jbrowse_url or jekyll_url()
    except JBrowseError as err:
        raise click.ClickException(str(err)) from err
    logger.info("Processing Collections...")
    parser = ProcessCollections(
        logger, jbrowse_url=jbrowse_url, out_dir=collections_out
    )
    logger.info("Outputting Collections...")
    parser.parse_collections(from_github, taxa_list)


@click.command()
@click.option(
    "--taxa_list",
    default=None,
    help="""Taxa.yml file restricting the genera processed. (Default: every genus in --from_github)""",
)
@click.option(
    "--nodes_out",
    default="./dscensor_nodes",
    help="""Output directory for DSCensor nodes. (Default: ./dscensor_nodes)""",
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
def populate_dscensor(taxa_list, nodes_out, from_github, log_file, log_level):
    """CLI entry for populate-dscensor"""
    logger = setup_logging(log_file, log_level, "populate-dscensor")
    parser = ProcessCollections(logger, out_dir=nodes_out)  # initialize class
    logger.info("Processing Collections...")
    parser.parse_collections(from_github, taxa_list)
    logger.info("Creating DSCensor Nodes...")
    parser.populate_dscensor(nodes_out)  # populate DSCensor nodes


@click.command()
@click.option(
    "--jbrowse_url",
    help="""Unused: populate-jekyll writes the JBrowse resource links. Accepted so
    existing invocations keep working.""",
)
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
    default="./jbrowse_out",
    help="""Output directory for Jbrowse2. (Default: ./jbrowse_out)""",
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
    """CLI entry for populate-jbrowse2

    Prints or runs the jbrowse commands that build LIS's JBrowse 2 config, planned
    offline from a datastore-metadata checkout.
    """
    logger = setup_logging(log_file, log_level, "populate-jbrowse2")
    if jbrowse_url:
        logger.warning("--jbrowse_url is unused; populate-jekyll writes those links")
    if not os.path.isdir(from_github):
        raise click.ClickException(
            f"{from_github} is not a datastore-metadata checkout"
        )
    if taxa_list and not os.path.exists(taxa_list):
        raise click.ClickException(f"taxon list {taxa_list} does not exist")
    index = DatastoreIndex(
        from_github, logger=logger, datastore_url=datastore_url
    ).build()
    entries = jbrowse.plan(index, jbrowse.genera(index, taxa_list))
    os.makedirs(jbrowse_out, exist_ok=True)
    for entry in entries:
        cmd = jbrowse.command(entry, jbrowse_out)
        if cmds_only:
            click.echo(cmd)
        else:
            subprocess.check_call(cmd, shell=True, executable="/bin/bash")
    logger.info("%s JBrowse entries for %s", len(entries), jbrowse_out)


@click.command()
@click.option(
    "--taxa_list",
    default="../_data/taxon_list.yml",
    help="""Taxa.yml file. (Default: ../_data/taxon_list.yml)""",
)
@click.option(
    "--blast_out",
    default="./blast_out",
    help="""Output directory for BLAST DBs. (Default: ./blast_out)""",
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
    parser.parse_collections(from_github, taxa_list)
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
    "--jbrowse_config",
    "jbrowse_configs",
    multiple=True,
    metavar="INSTANCE=PATH",
    callback=parse_pairs("INSTANCE=PATH"),
    help="""A JBrowse instance's deployed config.json, by its id in jbrowse.yml.
    Repeatable. With any given, the catalog records where each collection appears in
    LIS's JBrowse instances; an instance without one is marked unavailable.""",
)
@click.option(
    "--jbrowse_report",
    default=None,
    help="""Where to write a Markdown summary of the JBrowse instances: what each
    serves, the tracks that drift from the Data Store, and the plan vs the instance
    populate-jbrowse2 builds. Needs --jbrowse_config.""",
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
def populate_catalog(  # pylint: disable=too-many-arguments,too-many-positional-arguments
    from_github,
    catalog_out,
    datastore_url,
    indent,
    verify,
    jbrowse_configs,
    jbrowse_report,
    log_file,
    log_level,
):
    """CLI entry for populate-catalog

    Builds the whole-store catalog from a datastore-metadata checkout. Unlike the
    other subcommands this needs no network and no taxa list: it describes every
    collection the checkout knows about, in one document.
    """
    logger = setup_logging(log_file, log_level, "populate-catalog")
    builder = CatalogBuilder(
        from_github,
        logger=logger,
        datastore_url=datastore_url,
        verify=verify,
        jbrowse_configs=jbrowse_configs,
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
    try:
        builder.write(catalog_out, indent=indent or None)
    except JBrowseError as err:
        raise click.ClickException(str(err)) from err
    if jbrowse_report and builder.deployments:
        index = builder.index
        text = jbrowse.report(
            index,
            builder.deployments,
            jbrowse.plan(index, jbrowse.genera(index)),
            jbrowse.load_config()["built_instance"],
        )
        with open(jbrowse_report, "w", encoding="utf-8") as handle:
            handle.write(text)
        logger.info("wrote %s", jbrowse_report)


@click.command()
@click.option(
    "--collection",
    "collections",
    multiple=True,
    required=True,
    help="""Diversity collection to serve, by identifier (e.g. Wm82.gnm4.div.Song_Hyten_2015).
    Repeat for several; services are written in the order given.""",
)
@click.option(
    "--from_github",
    default="./datastore-metadata",
    help="""Path to datastore-metadata github directory. (Default: ./datastore-metadata).""",
)
@click.option(
    "--compose_out",
    default="./docker-compose.yml",
    help="""Where to write the compose file. Its services build from the Dockerfile beside
    it and its proxy reads traefik/ beside it, so write it to the root of a
    legumeinfo/divbrowse checkout. (Default: ./docker-compose.yml)""",
)
@click.option(
    "--datastore_url",
    default="https://data.legumeinfo.org",
    help="""URL the VCF and GFF3 links are built against.""",
)
@click.option(
    "--host",
    "hosts",
    multiple=True,
    metavar="GENUS=HOSTNAME",
    callback=parse_pairs("GENUS=HOSTNAME"),
    help="""Public hostname for a genus's services, served at http://<host>/<collection>/.
    Repeatable; adds to or overrides the defaults in divbrowse.yml.""",
)
@click.option(
    "--log_file",
    default="./populate-divbrowse.log",
    help="""Log file to output messages. (default: ./populate-divbrowse.log)""",
)
@click.option(
    "--log_level",
    default="INFO",
    help="""Log Level to output messages. (default: INFO)""",
)
def populate_divbrowse(
    collections,
    from_github,
    compose_out,
    datastore_url,
    hosts,
    log_file,
    log_level,
):
    """CLI entry for populate-divbrowse

    Writes a Divbrowse docker-compose.yml with one service per diversity collection,
    behind a Traefik proxy, built offline from a datastore-metadata checkout. A collection the format can't
    express stops the command without writing anything.
    """
    logger = setup_logging(log_file, log_level, "populate-divbrowse")
    index = DatastoreIndex(
        from_github, logger=logger, datastore_url=datastore_url
    ).build()
    try:
        defaults = load_config()["hosts"]
        text = compose_file(index, collections, hosts={**defaults, **hosts})
    except DivbrowseError as err:
        logger.error("not writing %s:\n%s", compose_out, err)
        raise click.ClickException(
            f"cannot build a Divbrowse service for:\n{err}"
        ) from err
    os.makedirs(os.path.dirname(os.path.abspath(compose_out)), exist_ok=True)
    with open(compose_out, "w", encoding="utf-8") as handle:
        handle.write(text)
    logger.info("wrote %s with %s services", compose_out, len(collections))
