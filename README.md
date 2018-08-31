# gSCAPE

_Global Landscape Evolution Model_

## Dependencies

### FillIT

+ [numpy](http://numpy.org)
+ [meshplex](https://github.com/nschloe/meshplex)
+ fortran compiler, preferably [gfortran](https://gcc.gnu.org/wiki/GFortran)

### Notebooks

To run the Examples (jupyter notebooks) it is recommended to use the `gscape-docker` image which has all the libraries already installed.

[https://hub.docker.com/u/geodels/dashboard/](https://hub.docker.com/u/geodels/dashboard/)

## Installation

```bash
python setup.py install
```

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
args = parser.parse_args()
print("Input file: {}".format(args.input))
print(" Verbose is on? {}".format(args.verbose))

# Reading input file
model = sim.LandscapeEvolutionModel(args.input,args.verbose)

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
