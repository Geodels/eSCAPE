# gSCAPE

_Global Landscape Evolution Model_

## Dependencies

### FillIT

+ [numpy](http://numpy.org)
+ [voropy](https://github.com/nschloe/voropy)
+ fortran compiler, preferably [gfortran](https://gcc.gnu.org/wiki/GFortran)

### Notebooks

To run the Examples (jupyter notebooks) it is recommended to use the `gscape-docker` image which has all the libraries already installed.

[https://hub.docker.com/u/geodels/dashboard/](https://hub.docker.com/u/geodels/dashboard/)

## Installation

```bash
python setup.py install
```

## Usage

Two classes are included as part of the **FillIT** package:

+ `depressionFillingZhou`
+ `depressionFillingBarnes`

These classes share similar methods for both structured and unstructured mesh and can be easily interchanged.

For regular grids, the following `classes` are available:
+ `depressionFillingZhou`
+ `depressionFillingBarnes`

To call one of these classes, you will typically do as follows:

``` python
import fillit as pitfill

# Class initialisation
pitClass = pitfill.depressionFillingBarnes(dem)

# Performing pit filling
fillZ = pitClass.performPitFillingSruct()

```
