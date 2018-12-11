###
# Copyright 2017-2018 Tristan Salles
#
# This file is part of eSCAPE.
#
# eSCAPE is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or any later version.
#
# eSCAPE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with eSCAPE.  If not, see <http://www.gnu.org/licenses/>.
###

"""
**LandscapeEvolutionModel** super class is the main entry point for eSCAPE and contains an inherited class (LandscapeEvolutionModelClass) which defines the following 3 main functions:

.. list-table::
    :header-rows: 1
    :widths: 10 60
    :stub-columns: 1
    :align: left

    *  -  Main functions
       -  Summary
    *  -  LandscapeEvolutionModel()
       -  Instantiates eSCAPE model object and performs surface processes evolution.
    *  -  runProcesses()
       -  Run eSCAPE Earth surface processes.
    *  -  destroy()
       -  Destroy PETSc DMPlex objects and associated Petsc local/global Vectors and Matrices.

A typical call to eSCAPE will be like this

.. code-block:: python
    :linenos:

    import eSCAPE as sim

    # Reading input file
    model = sim.LandscapeEvolutionModel('input_globe.yml',False,False)

    # Running model
    model.runProcesses()

    # Cleaning model
    model.destroy()

Submodules
==========

.. autosummary::
    :toctree: _autosummary

    mesher
    flow
    pit
    tools

"""

import petsc4py,sys
petsc4py.init(sys.argv)
from time import clock
from mpi4py import MPI
from .mesher import UnstMesh as _UnstMesh
from .pit import UnstPit as _UnstPit
from .tools import ReadYaml as _ReadYaml
from .tools import WriteMesh as _WriteMesh
from petsc4py import PETSc as _PETSc
from .flow import SPMesh as _SPMesh

import tools

MPIrank = MPI.COMM_WORLD.Get_rank()

def LandscapeEvolutionModel(filename, *args, **kwargs):
    """
    Instantiates eSCAPE model object and performs surface processes evolution.

    This object contains methods for the following operations:
     - initialisation of eSCAPE mesh based on input file options.
     - computation of surface processes
     - cleaning/destruction of PETSC objects

    Args
        filename : YAML input file
        verbose : True/False
            Output option for model main functions
        showlog : True/False
            Output option for PETSC logging file

    Returns:
        LandscapeEvolutionModel : object
    """

    class LandscapeEvolutionModelClass(_ReadYaml, _WriteMesh, _UnstMesh, _UnstPit, _SPMesh):

        def __init__(self, filename, verbose=True, showlog=False, *args, **kwargs):

            self.showlog = showlog
            if self.showlog:
                self.log = _PETSc.Log()
                self.log.begin()

            self.modelRunTime = clock()
            t_init = clock()
            self.verbose = verbose
            _ReadYaml.__init__(self, filename)

            _UnstMesh.__init__(self, self.meshFile, *args, **kwargs)

            _WriteMesh.__init__(self)

            # Surface processes initialisation
            _SPMesh.__init__(self,*args, **kwargs)

            # Pit filling algorithm initialisation
            _UnstPit.__init__(self,*args, **kwargs)

            # Get external forces
            _UnstMesh.applyForces(self)
            if MPIrank == 0:
                print('--- Initialisation Phase (%0.02f seconds)'% (clock() - t_init))

            return

        def runProcesses(self):
            """
            Run eSCAPE Earth surface processes.

            This function contains methods for the following operations:
             - calculating flow accumulation
             - erosion/deposition induced by stream power law
             - depression identification and pit filling
             - stream induced deposition diffusion
             - hillslope diffusion
            """

            while(self.tNow<=self.tEnd):
                tstep = clock()

                # Compute Flow Accumulation
                _SPMesh.FlowAccumulation(self)

                # Output time step for first step
                if self.tNow == self.tStart:
                    _WriteMesh.outputMesh(self, remesh=False)
                    self.saveTime += self.tout

                # Compute Stream Power Law
                _SPMesh.StreamPowerLaw(self)

                # Find depressions chracteristics (volume and spill-over nodes)
                if self.frac_fine < 1.:
                    _UnstPit.computeDepression(self)

                # Apply diffusion to deposited sediments
                if self.frac_fine < 1.:
                    _SPMesh.SedimentDiffusion(self)

                # _WriteMesh.outputMesh(self, remesh=False)
                # dedede
                # Compute Hillslope Diffusion Law
                _SPMesh.HillSlope(self)

                # Update Boundaries
                _UnstMesh.updateBoundaries(self)

                # Output time step
                if self.tNow >= self.saveTime:
                    _WriteMesh.outputMesh(self, remesh=False)
                    self.saveTime += self.tout

                # Update Tectonic, Sea-level & Climatic conditions
                if self.tNow<self.tEnd:
                    _UnstMesh.applyForces(self)

                # Advance time
                self.tNow += self.dt

                if MPIrank == 0:
                    print('--- Computational Step (%0.02f seconds)'% (clock() - tstep))

            return

        def destroy(self):
            """
            Destroy PETSc DMPlex objects and associated Petsc local/global Vectors and Matrices.
            Safely quit eSCAPE model.
            """

            _UnstMesh.destroy_DMPlex(self)

            if self.showlog:
                self.log.view()

            if MPIrank == 0:
                print('\n+++\n+++ Total run time (%0.02f seconds)\n+++'% (clock() - self.modelRunTime))

            return

    return LandscapeEvolutionModelClass(filename, *args, **kwargs)
