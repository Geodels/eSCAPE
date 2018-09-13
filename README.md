<div align="center">
    <img width=1000 src="https://github.com/Geodels/eSCAPE/blob/master/images/escape.jpg" alt="eSCAPE" title="eSCAPE Model"</img>
</div>

## Overview

**eSCAPE** is a parallel TIN-based landscape evolution model, built to simulate topography dynamic at various space and time scales. The model accounts for hillslope processes (soil creep using **linear** diffusion), fluvial incision (**stream power law**), spatially and temporally varying tectonics (vertical displacements) and climatic forces (temporal and spatial precipitation changes and/or sea-level fluctuations). 

## Getting started

For installation information and documentation visit our github [**wiki page**](https://github.com/Geodels/eSCAPE/wiki) which provides a quick guide on the installation dependencies.

A set of examples are available in the [eSCAPE-demo](https://github.com/Geodels/eSCAPE-demo) repository.

The easiest way to get started is with the Docker container
[https://hub.docker.com/u/geodels/](https://hub.docker.com/u/geodels/) using [Kitematic](https://docs.docker.com/kitematic/userguide/). Once **Kitematic** is installed on your computer, open it and look for **escape-docker** via the *search* menu.

If you want to install it yourself, you can follow the steps provided in the [wiki](https://github.com/Geodels/eSCAPE/wiki/Installation-on-HPC) page.

## The specs...

The model is based on the following approaches:
* an adaptation of the implicit, parallelizable method for calculating drainage area from [**Richardson & Perron (2014)**](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326),
* the extension of the parallel priority-flood algorithm from [**Barnes (2016)**](http://www.sciencedirect.com/science/article/pii/S0169555X12004618) to unstructured mesh,
* the methods developped in [**pyBadlands**](https://github.com/badlands-model/pyBadlands_serial) ([**Salles et al. (2018)**](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0195557)). 

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

To use/test **eSCAPE** quickly, it is recommended to use the `escape-docker` image that is shipped with all the required libraries.

[https://hub.docker.com/u/geodels/](https://hub.docker.com/u/geodels/)

## Installation

Once the libraries mentioned above have been installed, **eSCAPE** will need to be cloned and compiled using the following:

```bash
git clone https://github.com/Geodels/eSCAPE.git
python setup.py install
```

An example on how to install it on a _HPC server_ is provided in the [**wiki**](https://github.com/Geodels/eSCAPE/wiki/Installation-on-HPC) page.

## Usage

Either via _jupyter notebooks_ or _python_ files.

```bash
python run_eSCAPE.py -i input.yml -v
```

where the `run_eSCAPE.py` script takes one required argument the input filename and an optional verbose command (`-v`).  To run the script in parallel simply use the `mpirun` command. As an example with N processors it will look like:

```bash
mpirun -np N python run_eSCAPE.py -i input.yml
```

`run_eSCAPE.py` consists of a limited number of calls to **eSCAPE**

```python
import eSCAPE
model = eSCAPE.LandscapeEvolutionModel(***)
model.runProcesses()
model.destroy()
```

as shown below:

```python
import argparse
import eSCAPE as sim

# Parsing command line arguments
parser = argparse.ArgumentParser(description='This is a simple entry to run eSCAPE model.',add_help=True)
parser.add_argument('-i','--input', help='Input file name (YAML file)',required=True)
parser.add_argument('-v','--verbose',help='True/false option for verbose', required=False,action="store_true",default=False)
parser.add_argument('-l','--log',help='True/false option for PETSC log', required=False,action="store_true",default=False)

args = parser.parse_args()
if args.verbose:
  print("Input file: {}".format(args.input))
  print(" Verbose is on? {}".format(args.verbose))
  print(" PETSC log is on? {}".format(args.log))

# Reading input file
model = sim.LandscapeEvolutionModel(args.input,args.verbose,args.log)

# Running model
model.runProcesses()

# Cleaning model
model.destroy()
```

### Input file

Input files for **eSCAPE** are based on [YAML](https://circleci.com/blog/what-is-yaml-a-beginner-s-guide/) syntax.

A typical file will look like this:

```YAML
name: Description of the what is going to be done in this simulation...

domain:
    filename: ['data/inputfileparameters.vtu','Z']
    flowdir: 1

time:
    start: 0.
    end: 1000000.
    tout: 1000.
    dt: 100.

sea:
    position: 0.
    curve: 'data/sealevel.csv'

climate:
    - start: 0.
      uniform: 1.0
    - start: 500000.
      map: ['data/inputfileparameters.vtu','R']
    - start: 500000.
      uniform: 2.0

tectonic:
    - start: 0.
      map: ['data/inputfileparameters.vtu','T1']
    - start: 100000.
      uniform: 0.
    - start: 50000.
      map: ['data/inputfileparameters.vtu','T2']

spl:
    m: 0.5
    n: 1.0
    Ke: 1.e-5

diffusion:
    hillslopeK: 5.e-2
    streamK: 300.
    oceanK: 100.
    maxIT: 2000

output:
    dir: 'outputDir'
    makedir: False

```

### Tutorials

To get some additional info in regards to how to use **eSCAPE** a series of examples and tutorials is provided in the docker container (`escape-docker`) and is also available for  download from the [eSCAPE-demo](https://github.com/Geodels/eSCAPE-demo) repository.
