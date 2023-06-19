# PhyloX

PhyloX is a python package with tools for constructing, manipulating, and analyzing phylogenetic networks.

## Install

Install as pypi package phylox:
```
pip install phylox
```

## Usage

You can load the package methods with `import phylox` in python.

## Development

### Linting

precommit is yet to be configured, for now, simply run black and isort.

### Release

set new version number in master branch
 - CHANGELOG.md

release current version
```
git checkout main
git tag [version number]
git checkout release
git merge master
python -m build
python -m twine upload --repository pypi dist/*
```
set new version number in master branch
 - pyproject.toml
