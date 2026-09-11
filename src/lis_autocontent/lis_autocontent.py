"""Entry point executable for the CLI. Gets added to PATH on python setup.py install"""

import click
from . import lis_cli


@click.group()
def cli():
    """Click Group Object to Attach Commands."""
    pass


def cli_entry_point():
    """CLI Eentry Point for package build and install."""
    cli.add_command(lis_cli.populate_jekyll)
    cli.add_command(lis_cli.populate_jbrowse2)
    cli.add_command(lis_cli.populate_blast)
    cli.add_command(lis_cli.populate_dscensor)
    cli.add_command(lis_cli.populate_catalog)
    cli()  # invoke cli


if __name__ == "__main__":
    cli_entry_point()
