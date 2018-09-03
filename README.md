# gSCAPE

_Global Landscape Evolution Model_

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

To use/test **gSCAPE** quickly, it is recommended to use the `gscape-docker` image that is shipped with all the required libraries.

[https://hub.docker.com/u/geodels/](https://hub.docker.com/u/geodels/)

## Installation

Once the libraries mentioned above have been installed, **gSCAPE** will need to be cloned and compiled using the following:

```bash
git clone https://github.com/Geodels/gSCAPE.git
python setup.py install
```

An example on how to install it on HPC server is provided in the [wiki](https://github.com/Geodels/gSCAPE/wiki/Installation-on-HPC) page.

## Usage

Either via _jupyter notebooks_ or _python_ files.

```bash
python run_gSCAPE.py -i input.yml -v
```

where the `run_gSCAPE.py` script takes one required argument the input filename and an optional verbose command (`-v`).  To run the script in parallel simply use the `mpirun` command. As an example with N processors it will look like:

```bash
mpirun -np N python run_gSCAPE.py -i input.yml
```

`run_gSCAPE.py` consists of a limited number of calls to **gSCAPE**

```python
import gSCAPE
model = gSCAPE.LandscapeEvolutionModel(***)
model.runProcesses()
model.destroy()
```

as shown below:

```python
import argparse
import gSCAPE as sim

# Parsing command line arguments
parser = argparse.ArgumentParser(description='This is a simple entry to run gSCAPE model.',add_help=True)
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

Input files for **gSCAPE** are based on [YAML](https://circleci.com/blog/what-is-yaml-a-beginner-s-guide/) syntax.

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
    step: 50

output:
    dir: 'outputDir'
    makedir: False

```

### Tutorials

To get some additional info in regards to how to use **gSCAPE** a series of examples and tutorials is provided in the docker container (`gscape-docker`) and is also available for  download from the [gSCAPE-demo](https://github.com/Geodels/gSCAPE-demo) repository.
