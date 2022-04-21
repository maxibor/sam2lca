# Contributing

We welcome any contributions !

To further develop sam2lca, or add documentation, please read further:

## Clone the sam2lca repository, and checkout the dev branch

```bash
git clone git@github.com:maxibor/sam2lca.git
git checkout dev
```

## Install and activate the development environment

```bash
conda env create -f environment.yml
conda activate sam2lca
```

## Install sam2lca with pip in editable mode

```bash
pip install -e .
```

### Build the documentation

```bash
cd docs
make html
```

The docs are built in the `docs/build/html` directory
