#!/usr/bin/env python3

import string
import click
from sam2lca import __version__
from sam2lca.main import sam2lca, update_database, list_available_db
from sam2lca.utils import taxo_ranks
from pathlib import Path


from click import option, Option, UsageError


class MutuallyExclusiveOption(Option):
    # Credits goes to Stan Chang for this code snippet
    # https://gist.github.com/stanchan/bce1c2d030c76fe9223b5ff6ad0f03db

    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(kwargs.pop("mutually_exclusive", []))
        help = kwargs.get("help", "")
        if self.mutually_exclusive:
            ex_str = ", ".join(self.mutually_exclusive)
            kwargs["help"] = help + (
                " NOTE: This argument is mutually exclusive with "
                " arguments: [" + ex_str + "]."
            )
        super(MutuallyExclusiveOption, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.mutually_exclusive.intersection(opts) and self.name in opts:
            raise UsageError(
                "Illegal usage: `{}` is mutually exclusive with "
                "arguments `{}`.".format(self.name, ", ".join(self.mutually_exclusive))
            )

        return super(MutuallyExclusiveOption, self).handle_parse_result(ctx, opts, args)


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
@option(
    "-i",
    "--identity",
    cls=MutuallyExclusiveOption,
    type=click.FloatRange(0, 1),
    default=0.8,
    show_default=True,
    help="Minimum identity threshold",
    mutually_exclusive=["distance"],
)
@option(
    "-d",
    "--distance",
    cls=MutuallyExclusiveOption,
    type=int,
    default=None,
    help="Edit distance threshold",
    mutually_exclusive=["identity"],
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
    help="sam2lca output file. Default: [basename].sam2lca.*",
)
@click.option(
    "-b",
    "--bam_out",
    is_flag=True,
    help="Write BAM output file with XT tag for TAXID",
)
@click.option(
    "-r",
    "--bam_split_rank",
    type=click.Choice(taxo_ranks, case_sensitive=False),
    help="Write BAM output file split by TAXID at rank -r. To use in combination with -b/--bam_out",
)
@click.option(
    "-n",
    "--bam_split_read",
    type=int,
    default=50,
    show_default=True,
    help="Minimum numbers of reads to write BAM split by TAXID. To use in combination with -b/--bam_out",
)
def analyze(ctx, no_args_is_help=True, **kwargs):
    """\b
    Run the sam2lca analysis

    SAM: path to SAM/BAM/CRAM alignment file
    """
    sam2lca(**kwargs, **ctx.obj)


@cli.command(help="Download/prepare acc2tax and taxonomy databases")
@click.pass_context
@click.option(
    "-t",
    "--taxonomy",
    type=str,
    default="ncbi",
    show_default=True,
    help="Name of taxonomy database to create (ncbi | gtdb)",
)
@click.option(
    "--taxo_names",
    type=click.Path(readable=True, dir_okay=False, file_okay=True),
    default=None,
    help="names.dmp file for Taxonomy database (optional). Only needed for custom taxonomy database (non ncbi or gtdb)",
)
@click.option(
    "--taxo_nodes",
    type=click.Path(readable=True, dir_okay=False, file_okay=True),
    default=None,
    help="nodes.dmp file for Taxonomy database (optional). Only needed for custom taxonomy database (non ncbi or gtdb)",
)
@click.option(
    "--taxo_merged",
    type=click.Path(readable=True, dir_okay=False, file_okay=True),
    default=None,
    help="merged.dmp file for Taxonomy database (optional). Only needed for custom taxonomy database (non ncbi or gtdb)",
)
@click.option(
    "-a",
    "--acc2tax",
    type=str,
    default=None,
    show_default=True,
    help="Type of acc2tax mapping database to build."
)
@click.option(
    "--acc2tax_json",
    type=click.Path(),
    default="https://raw.githubusercontent.com/maxibor/sam2lca/master/data/acc2tax.json",
    show_default=True,
    help="JSON file for specifying extra acc2tax mappings",
)
def update_db(ctx, no_args_is_help=True, **kwargs):
    update_database(**ctx.obj, **kwargs)


@cli.command(help="List available taxonomy and acc2tax databases")
@click.pass_context
def list_db(ctx, no_args_is_help=True):
    list_available_db(**ctx.obj, verbose=True)


if __name__ == "__main__":
    cli()
