[![DOI](https://zenodo.org/badge/146602973.svg)](https://zenodo.org/badge/latestdoi/146602973)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00964/status.svg)](https://doi.org/10.21105/joss.00964)

<div align="center">
    <img width=1000 src="https://github.com/Geodels/eSCAPE/blob/master/images/bg.jpg" alt="eSCAPE" title="eSCAPE Model"</img>
</div>

## Overview

**eSCAPE** is a parallel TIN-based landscape evolution model, built to simulate topography dynamic at various space and time scales. The model accounts for hillslope processes (soil creep using **linear** diffusion), fluvial incision (**stream power law**), spatially and temporally varying tectonics (horizontal & vertical displacements) and climatic forces (temporal and spatial precipitation changes and/or sea-level fluctuations).

## Getting started

For installation information and documentation visit our github [**eSCAPE-model** website](https://escape-model.github.io) or the [**wiki page**](https://github.com/Geodels/eSCAPE/wiki) which provides a quick guide on the installation dependencies and **Docker** settings.

API documentation is available from the [eSCAPE-API](https://escape-api.github.io/index.html) website.

A set of examples are available in the [eSCAPE-demo](https://github.com/Geodels/eSCAPE-demo) repository.

The easiest way to get started is with the Docker container
[https://hub.docker.com/u/geodels/](https://hub.docker.com/u/geodels/) using [Kitematic](https://docs.docker.com/kitematic/userguide/). Once **Kitematic** is installed on your computer, open it and look for **Geodels escape-docker** via the *search* menu.

If you want to install it yourself, you can follow the steps provided in the [wiki](https://github.com/Geodels/eSCAPE/wiki/Installation-on-HPC) page.

<div align="center">
<img src="https://github.com/Geodels/eSCAPE-hugo/blob/master/static/images/earthdepo1.gif" width="600"/>
</div>

## The specs...

The model is based on the following approaches:
* an adaptation of the implicit, parallelizable method for calculating drainage area for both single (D8) and multiple flow direction (Dinf) from [**Richardson & Perron (2014)**](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326),
* the methods developed in [**pyBadlands**](https://github.com/badlands-model/pyBadlands_serial) ([**Salles et al. (2018)**](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195557)) for marine sedimentation.

### Community driven

To join the eSCAPE User Group on Slack, send an email request to: <a href="MAILTO:tristan.salles@sydney.edu.au?subject=eSCAPE User Group&body=Please send me an invite to join the eSCAPE User Group">Join eSCAPE User Group</a>

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/lgpl-3.0.en.html>.

<div align="center">
    <img width=1000 src="https://github.com/Geodels/eSCAPE/blob/master/images/escapezoom.jpg" alt="eSCAPE" title="eSCAPE Model"</img>
</div>


## Dependencies

1- Scientific computing libraries:
+ [numpy](http://www.numpy.org)
+ [scipy](https://www.scipy.org)
+ [mpi4py](https://mpi4py.readthedocs.io/en/stable/)
+ [petsc4py](https://petsc4py.readthedocs.io/en/stable/)
+ fortran compiler, preferably [gfortran](https://gcc.gnu.org/wiki/GFortran)
+ [fillit](https://github.com/Geodels/fillit)

2- Reading/Writing/Parsing libraries:
+ [ruamel.yaml](https://yaml.readthedocs.io/en/latest/)
+ [pandas](https://pandas.pydata.org)
+ [meshio](https://github.com/nschloe/meshio)
+ [h5py](https://www.h5py.org)
+ [meshplex](https://github.com/nschloe/meshplex)

### Docker

To use/test **eSCAPE** quickly, it is recommended to use the `Geodels escape-docker` image that is shipped with all the required libraries.

[https://hub.docker.com/u/geodels/](https://hub.docker.com/u/geodels/)

## Installation

Once the libraries mentioned above have been installed, **eSCAPE** will need to be cloned and compiled using the following:

```bash
git clone https://github.com/Geodels/eSCAPE.git
python setup.py install
```

An example on how to install it on a _HPC server_ is provided in the [**wiki**](https://github.com/Geodels/eSCAPE/wiki/Installation-on-HPC) page.

### Testing the installation

A `test` model is available in the [eSCAPE-demo](https://github.com/Geodels/eSCAPE-demo) repository and can be used to test your **eSCAPE** installation. It provides two comparisons related to:

+ the _expected values_ that you should obtained if your installation is successful
+ runtime for both _serial_ and _parallel_ simulations

This test should take less than 1 minute. An example on how to visualise **eSCAPE** outputs in [**Paraview**](https://www.paraview.org/download/) is also provided.


### Performance

The code parallelisation relies on `petsc` and scalability tests are still on-going. Our goal is to be able to run models with more than 10 millions nodes and over several hundreds of CPUs.

## Tutorials

To get some additional info in regards to how to use **eSCAPE** a series of examples and tutorials is provided in the docker container (`escape-docker`) and is also available for download from the [eSCAPE-demo](https://github.com/Geodels/eSCAPE-demo) repository.

## Contributions

We welcome all kinds of contributions! Please get in touch if you would like to help out. Everything from code to notebooks to examples and documentation are all equally valuable so please don't feel you can't contribute. To contribute please fork the project make your changes and submit a pull request. We will do our best to work through any issues with you and get your code merged into the main branch.

If you found a bug, have questions, or are just having trouble with **eSCAPE**, you can:
+ join the _eSCAPE User Group on Slack_ by sending an email request to: <a href="MAILTO:tristan.salles@sydney.edu.au?subject=eSCAPE User Group&body=Please send me an invite to join the eSCAPE User Group">Join eSCAPE User Group</a>
+ open an issue in our [issue tracker](https://github.com/Geodels/eSCAPE/issues/new) and we'll try to help resolve the concern.
