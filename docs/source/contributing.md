# Contributing

We welcome any contributions !

To further develop sam2lca, or add documentation, please read further:

## 1. Clone the sam2lca repository, and checkout the dev branch

```bash
git clone git@github.com:maxibor/sam2lca.git
git checkout dev
```

## 2. Install and activate the sam2lca conda environment

```bash
conda env create -f environment.yml
conda activate sam2lca
```

> For all the following steps, we strongly suggest you to work in the `sam2lca` conda environment

## 3. Install sam2lca with pip in editable mode

```bash
pip install -e .
```

## Run the unit and integration tests

```bash
pytest -s -vv --script-launch-mode=subprocess
```

## Build the documentation

```bash
cd docs
make html
```

The docs are built in the `docs/build/html` directory

### Claim your sticker

Thanks for contributing to sam2lca !
If you want to spread the word about sam2lca, please get in touch with me to claim your sticker (maxime_borry[at]eva.mpg.de) !
