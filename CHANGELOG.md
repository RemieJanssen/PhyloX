# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Types of changes:
- Added for new features.
- Changed for changes in existing functionality.
- Deprecated for soon-to-be removed features.
- Removed for now removed features.
- Fixed for any bug fixes.
- Security in case of vulnerabilities.

## [0.0.3] - (Unreleased)

### Added

- Network properties methods: b2_balance, count_cherries_and_reticulated_cherries and blob_properties. (#11)
- Tree-child-sequence functions. (#12)
- Usage of NetworkX referenced in README. (#14)
- Cherry picking hybridization workflow works again, including tests. (#15)
- Documentation using sphinx. (#18)

### Changed

- Cleaned up generators: beta-splitting trees, add_edges, caterpillar, balanced tree, lgt, zods, random tree-child sequence. (#14)

## [0.0.2] - (2023-06-19)

### Added

- Release and version tag explanation. (#3)
- Automated tests. (#4)

### Changed

- Run black and isort. (#4)

## [0.0.1] - (2023-06-18)

### Changed

- cleaned up DiNetwork class.
- cleaned up basic rearrangement methods.
- cleaned up isomorphism testing.
- cleand up checking for classes, including some cherry-picking methods.

### Added

- Migrated messy old code to this repo.
- Automorphism count: count_automorphisms.
- tests for DiNetwork, network classes, rearrangement moves, isomorphism checking, mcmc network generation.
- configuration of pytest with --cov.
