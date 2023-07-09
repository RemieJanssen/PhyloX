# PhyloX

PhyloX is a python package with tools for constructing, manipulating, and analyzing phylogenetic networks.

## Use of NetworkX
The implementation of PhyloX is based on NetworkX (NetworkX is distributed with the [3-clause BSD license](https://networkx.org/documentation/stable/index.html#license)):

> Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring network structure, dynamics, and function using NetworkX”, in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008

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
 - pyproject.toml

release current version
```
git checkout main
git tag [version number]
git checkout release
git merge main
python -m build
python -m twine upload --repository pypi dist/*
```