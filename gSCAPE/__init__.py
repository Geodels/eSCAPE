"""
Copyright 2017-2018 Tristan Salles

This file is part of gSCAPE.

gSCAPE is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

gSCAPE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with gSCAPE.  If not, see <http://www.gnu.org/licenses/>.
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
    Instantiates gSCAPE model object and performs surface processes evolution.

    This object contains methods for the following operations:
     - initialisation of gSCAPE mesh based on input file options.
     - computation of surface processes
     - cleaning/destruction of PETSC objects

    Parameters
    ----------
     filename : YAML input file
     verbose : True/False
        Output option for model main functions
    showlog : True/False
        Output option for PETSC logging file

    Returns
    -------
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
            Run gSCAPE Earth surface processes.

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
                _UnstPit.depressionDefinition(self)

                # Apply diffusion to deposited sediments
                _SPMesh.SedimentDiffusion(self)

                # Compute Hillslope Diffusion Law
                _SPMesh.HillSlope(self)

                # Output time step
                if self.tNow >= self.saveTime:
                    _WriteMesh.outputMesh(self, remesh=False)
                    self.saveTime += self.tout

                # Advance time
                self.tNow += self.dt

                # Update Tectonic, Sea-level & Climatic conditions
                _UnstMesh.applyForces(self)

                if MPIrank == 0:
                    print('--- Computational Step (%0.02f seconds)'% (clock() - tstep))

            return

        def destroy(self):
            """
            Destroy PETSc DMPlex objects and associated Petsc local/global Vectors and Matrices.
            Safely quit gSCAPE model.
            """

            _UnstMesh.destroy_DMPlex(self)

            if self.showlog:
                self.log.view()

            if MPIrank == 0:
                print('\n+++\n+++ Total run time (%0.02f seconds)\n+++'% (clock() - self.modelRunTime))

            return

    return LandscapeEvolutionModelClass(filename, *args, **kwargs)
