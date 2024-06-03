[![DOI](https://zenodo.org/badge/650049974.svg)](https://zenodo.org/doi/10.5281/zenodo.10122403)
[![JOSS](https://joss.theoj.org/papers/d63bf1aae710fd400a2aba922b359cd7/status.svg)](https://joss.theoj.org/papers/d63bf1aae710fd400a2aba922b359cd7)
[![Documentation Status](https://readthedocs.org/projects/phylox/badge/?version=latest)](https://phylox.readthedocs.io/en/latest/?badge=latest)


# PhyloX

PhyloX is a python package with tools for constructing, manipulating, and analyzing phylogenetic networks.

Installation as pypi package phylox is simple via pip or conda:
```
pip install phylox
```
You can load the package methods with `import phylox` in python.

For more information, see the [documentation](https://phylox.readthedocs.io).

## Citing PhyloX

For now, simply refer to the github page to cite PhyloX. Registering the software for a DOI is still on the to do list.

### Use of NetworkX
The implementation of PhyloX is based on NetworkX (NetworkX is distributed with the [3-clause BSD license](https://networkx.org/documentation/stable/index.html#license)):

> Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring network structure, dynamics, and function using NetworkX”, in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008

### Citing specific functions
When citing PhyloX, you are most likely also using specific methods, which can be traced back to their original papers. Please take care to cite the original papers as well. A reference to the original paper should be found in the documentation of the method, or of the module containing the method.

## Development

For a development version, simply pull the project and in the home of the project do:
```
pip install -e .
```
This installs the phylox package from the source. When you change things in the source, the package gets updated as well.

### Release

set new version number in master branch
 - CHANGELOG.md
 - pyproject.toml

release current version
```
git checkout release
git merge main
git tag [version number]
git push --atomic origin release [version number]
```

### Linting

precommit is yet to be configured, for now, simply run black and isort.

### Documentation

Documentation is in the docs folder, and is created uses sphinx.

#### Requirements
You may need to install the requirements from `docs/requirements.txt` first. Make sure to use python<=3.11.*, for example:
```
conda create -n phylox-sphinx
conda activate phylox-sphinx
python install python==3.11.*
pip install -r docs/requirements.txt
```

#### Creating documentation
to build the documentation, go to the docs folder and run:
```
make html
```
the docs will be in `docs/build/html`.

If you re-run the build, you can first remove the old autosummary files. If you do not, it will not update them.



