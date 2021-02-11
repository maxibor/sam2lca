from setuptools import setup, find_packages
import codecs
import os.path


def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()


def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name='sam2lca',
    version=get_version("sam2lca/__init__.py"),
    description='Last Commont Ancestor on SAM/BAM/CRAM alignment files',
    long_description=open("README.md").read(),
    url='https://github.com/maxibor/sam2lca',
    long_description_content_type="text/markdown",
    license='GNU-GPLv3',
    python_requires=">=3.6",
    install_requires=[
        'click',
        'pysam',
        'tqdm',
        'xopen',
        'python-rocksdb',
        'ete3',
        'pandas'
    ],
    packages=find_packages(include=['sam2lca']),
    entry_points={
        'console_scripts': [
            'sam2lca = sam2lca.cli:cli'
        ]
    },
    include_package_data=True,
    package_data={"accession": ["data/test.accession2taxid.gz"],
                  "md5":["data/test.accession2taxid.gz.md5"]
    }
)
