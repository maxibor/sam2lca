# Changelog

All notable changes to [sam2lca](https://github.com/maxibor/sam2lca) will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.0.0

### Added

- Added sam2lca tutorial
- Add Custom acc2tax with json
- Total number of reads is computed early on to provide progress bar
- unit and integration testing with PyTest
- Total Descendant read counts for each taxon
- GTDB taxonomy and acc2tax added
- 18s acc2tax added
- XN and XR flag in bam output for, respectively, Taxon name and rank
- Add edit distance threshold filtering

### Changed

- Code refactoring for speedup, making use of multithreading on shared dictionaries
- Improve logging and replace prints statements with logging.info
- ete3 has been replaced by taxopy
- Unclassified TAXID is now `12908` by default
- RocksDB params changed
- TAXID of LCA is only attributed to alignment segments passing threshold
- Identity threshold is now a range selection in CLI.

## Dependencies

- ete3 removed
- ordered-set removed
- pytest-console-scripts 1.3.1
- scipy
- Jinja2 pinned version ot 3.1 (see [RTD issue](https://github.com/readthedocs/readthedocs.org/issues/9038) )

### Removed
