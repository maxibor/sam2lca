#!/usr/bin/env python3

import click
from sam2lca import __version__
from sam2lca.main import sam2lca


@click.command()
@click.version_option(__version__)
@click.argument('sam', type=click.Path(exists=True))
@click.option('-m',
              '--mappings',
              type=click.Choice(['nucl_gb',
                                 'nucl_wgs',
                                 'pdb',
                                 'prot'],
                                case_sensitive=False),
              default='nucl_gb',
              show_default=True,
              help='Mapping type of accession to TAXID')
@click.option('-t',
              '--tree',
              type=click.Path(exists=True),
              help='Optional Newick Taxonomy Tree')
@click.option('-p',
              '--process',
              type=int,
              default=2,
              show_default=True,
              help='Number of process for parallelization')
@click.option('-u',
              '--update',
              is_flag=True,
              help='Update local copy of NCBI taxonomy database')
def cli(no_args_is_help=True, **kwargs):
    """\b
    sam2lca: Last Common Ancestor on SAM/BAM/CRAM alignment files
    Author: Maxime Borry
    Contact: <borry[at]shh.mpg.de>
    Homepage & Documentation: github.com/maxibor/sam2lca

    SAM: path to SAM/BAM/CRAM alignment file
    """
    sam2lca(**kwargs)


if __name__ == "__main__":
    cli()
