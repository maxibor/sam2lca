#!/usr/bin/env python3

import click
from sam2lca import __version__
from sam2lca.main import sam2lca, update_database
from pathlib import Path


@click.group()
@click.version_option(__version__)
@click.pass_context
@click.option(
    "-m",
    "--mappings",
    type=click.Choice(["nucl", "prot", "test"], case_sensitive=False),
    default="nucl",
    show_default=True,
    help="Mapping type of accession to TAXID",
)
@click.option(
    "-d",
    "--dbdir",
    type=click.Path(writable=True, dir_okay=True, file_okay=False),
    default=f"{str(Path.home())}/.sam2lca",
    show_default=True,
    help="Directory to store taxonomy databases",
)
def cli(ctx, mappings, dbdir):
    """\b
    sam2lca: Last Common Ancestor on SAM/BAM/CRAM alignment files
    Author: Maxime Borry
    Contact: <borry[at]shh.mpg.de>
    Homepage & Documentation: github.com/maxibor/sam2lca
    """
    ctx.ensure_object(dict)

    ctx.obj["mappings"] = mappings
    ctx.obj["dbdir"] = dbdir
    pass


@cli.command()
@click.pass_context
@click.argument("sam", type=click.Path(exists=True))
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
    "-t", "--tree", type=click.Path(exists=True), help="Optional Newick Taxonomy Tree"
)
@click.option(
    "-o",
    "--output",
    type=click.Path(writable=True, dir_okay=False, file_okay=True),
    default=None,
    help="sam2lca output file",
)
def analyze(ctx, no_args_is_help=True, **kwargs):
    """\b
    Run the sam2lca analysis

    SAM: path to SAM/BAM/CRAM alignment file
    """
    sam2lca(**kwargs, **ctx.obj)


@cli.command()
@click.pass_context
@click.option("-n", "--ncbi", is_flag=True, help="Update NCBI taxonomy tree")
def update_db(ctx, no_args_is_help=True, **kwargs):
    """Download/prepare mappings and taxonomy databases"""
    update_database(**ctx.obj, **kwargs)


if __name__ == "__main__":
    cli()
