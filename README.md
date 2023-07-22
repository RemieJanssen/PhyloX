# PhyloX

PhyloX is a python package with tools for constructing, manipulating, and analyzing phylogenetic networks.

For more information, see the [documentation](https://phylox.readthedocs.io).

## Citing PhyloX

For now, simply refer to the github page to cite PhyloX. Registering the software for a DOI is still on the to do list.

### Use of NetworkX
The implementation of PhyloX is based on NetworkX (NetworkX is distributed with the [3-clause BSD license](https://networkx.org/documentation/stable/index.html#license)):

> Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring network structure, dynamics, and function using NetworkX”, in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008

### Citing specific functions
When citing PhyloX, you are most likely also using specific methods, which can be traced back to their original papers. Please take care to cite the original papers as well. A reference to the original paper should be found in the documentation of the method, or of the module containing the method.

## Install

Install as pypi package phylox:
```
pip install phylox
```

## Usage

You can load the package methods with `import phylox` in python.

## Documentation

Documentation is in the docs folder, and is created uses sphinx.
to build the documentation, go to the docs folder and do:
```
make html
```
the docs will be in `docs/build/html`.

## Development

### Linting

precommit is yet to be configured, for now, simply run black and isort.

### Release

set new version number in master branch
 - CHANGELOG.md
 - pyproject.toml

release current version
```
git checkout release
git merge main
git tag [version number]
git push
```