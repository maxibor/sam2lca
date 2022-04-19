#!/usr/bin/env python3

import click
from sam2lca import __version__
from sam2lca.main import sam2lca, update_database, list_available_db
from sam2lca.acc2tax import base_map_config
from pathlib import Path


@click.group()
@click.version_option(__version__)
@click.pass_context
@click.option(
    "-d",
    "--dbdir",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
    default=f"{str(Path.home())}/.sam2lca",
    show_default=True,
    help="Directory to store taxonomy databases",
)
def cli(ctx, dbdir):
    """\b
    sam2lca: Lowest Common Ancestor on SAM/BAM/CRAM alignment files
    Author: Maxime Borry, Alexander Huebner
    Contact: <maxime_borry[at]eva.mpg.de>
    Homepage & Documentation: github.com/maxibor/sam2lca
    """
    ctx.ensure_object(dict)
    ctx.obj["dbdir"] = dbdir
    pass


@cli.command()
@click.pass_context
@click.argument("sam", type=click.Path(exists=True))
@click.option(
    "-t",
    "--taxonomy",
    type=str,
    default="ncbi",
    show_default=True,
    help="Taxonomy database to use",
)
@click.option(
    "-a",
    "--acc2tax",
    type=str,
    default="nucl",
    show_default=True,
    help="acc2tax database to use",
)
@click.option(
    "-i",
    "--identity",
    type=float,
    default=0.8,
    show_default=True,
    help="Minimum identity",
)
@click.option(
    "-l",
    "--length",
    type=int,
    default=30,
    show_default=True,
    help="Minimum alignment length",
)
@click.option(
    "-c",
    "--conserved",
    is_flag=True,
    help="Ignore reads mapping in ultraconserved regions",
)
@click.option(
    "-p",
    "--process",
    type=int,
    default=2,
    show_default=True,
    help="Number of process for parallelization",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
    default=None,
    help="sam2lca output file",
)
@click.option(
    "-b",
    "--bam_out",
    is_flag=True,
    help="Write BAM output file with XT tag for TAXID",
)
def analyze(ctx, no_args_is_help=True, **kwargs):
    """\b
    Run the sam2lca analysis

    SAM: path to SAM/BAM/CRAM alignment file
    """
    sam2lca(**kwargs, **ctx.obj)


@cli.command()
@click.pass_context
@click.option(
    "--taxonomy",
    type=str,
    default="ncbi",
    show_default=True,
    help="Name of taxonomy database to create",
)
@click.option(
    "--taxo_names",
    type=click.Path(readable=True, dir_okay=False, file_okay=True),
    default=None,
    help="names.dmp file for Taxonomy database (optional). Only needed for custom taxonomy database (non ncbi)",
)
@click.option(
    "--taxo_nodes",
    type=click.Path(readable=True, dir_okay=False, file_okay=True),
    default=None,
    help="nodes.dmp file for Taxonomy database (optional). Only needed for custom taxonomy database (non ncbi)",
)
@click.option(
    "--taxo_merged",
    type=click.Path(readable=True, dir_okay=False, file_okay=True),
    default=None,
    help="merged.dmp file for Taxonomy database (optional). Only needed for custom taxonomy database (non ncbi)",
)
@click.option(
    "--acc2tax",
    type=click.Choice(
        list(base_map_config["mapfiles"].keys()) + ["json"],
        case_sensitive=False,
    ),
    default=None,
    show_default=True,
    help="Type of acc2tax mapping database to build. Choose json to build custom acc2tax database from JSON file",
)
@click.option(
    "--acc2tax_json",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
    default=None,
    show_default=True,
    help="(Optional) JSON file for specifying extra acc2tax mappings",
)
def update_db(ctx, no_args_is_help=True, **kwargs):
    """Download/prepare acc2tax and taxonomy databases"""
    update_database(**ctx.obj, **kwargs)


@cli.command()
@click.pass_context
def list_db(ctx, no_args_is_help=True):
    """List available taxonomy and acc2tax databases"""
    list_available_db(**ctx.obj, verbose=True)


if __name__ == "__main__":
    cli()
