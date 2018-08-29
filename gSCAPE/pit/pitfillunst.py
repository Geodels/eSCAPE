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

import numpy as np
import pandas as pd
from mpi4py import MPI
import sys,petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from time import clock
import warnings;warnings.simplefilter('ignore')

import fillpy as fillAlgo
from gSCAPE._fortran import fillDepression
from gSCAPE._fortran import pitVolume
from gSCAPE._fortran import pitHeight
from gSCAPE._fortran import spillPoints

MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = PETSc.COMM_WORLD

try: range = xrange
except: pass

class UnstPit(object):
    """
    Building the priority flooding algorithm for depression determination
    """
    def __init__(self, *args, **kwargs):

        t0 = clock()
        self.first = 2

        self.pitData = None

        self.fillGlobal = self.dm.createGlobalVector()
        self.fillLocal = self.dm.createLocalVector()

        self.pitGlobal = self.dm.createGlobalVector()
        self.pitLocal = self.dm.createLocalVector()

        self.watershedGlobal = self.dm.createGlobalVector()
        self.watershedLocal = self.dm.createLocalVector()

        # Construct pit filling algorithm vertex indices
        vIS = self.dm.getVertexNumbering()

        # Local mesh points used in the pit filling algo
        self.idLocal =  np.where(vIS.indices>=0)[0]
        self.inIDs = np.zeros(self.npoints,dtype=int)
        self.inIDs[self.idLocal] = 1
        masknodes = np.isin(self.lcells, self.idLocal)
        tmp = np.sum(masknodes.astype(int),axis=1)
        out = np.where(np.logical_and(tmp>0,tmp<3))[0]
        ids = np.invert(masknodes[out]).flatten()
        vIS.destroy()

        # Local points that will be updated by the neighboring partition
        idComm = np.unique(self.lcells[out].flatten()[ids])
        self.idComm = np.zeros(self.npoints,dtype=int)
        self.idComm[idComm] = 1

        # Local points that are part of the global mesh boundary
        self.idGBounds = np.where(np.isin(self.idLocal,self.localboundIDs))[0]
        ids = masknodes[out].flatten()
        self.gbounds = np.zeros(self.npoints,dtype=int)
        self.gbounds[self.idGBounds] = 1

        # Border global ID
        tmpL = self.hLocal.duplicate()
        tmpG = self.hGlobal.duplicate()
        tmpL.setArray(self.gbounds)
        self.dm.localToGlobal(tmpL, tmpG, 1)
        GIDs =  tmpG.getArray()
        self.boundGIDs = np.where(GIDs>0)
        del GIDs
        tmpG.destroy()
        tmpL.destroy()

        # Local points that will be used to update watershed spill over points
        self.idLocalComm = np.unique(self.lcells[out].flatten()[ids])

        # Local points that represent the local pit filling mesh boundary
        self.idLBounds = np.unique(np.concatenate((self.idLocalComm,self.idGBounds)))

        if MPIrank == 0 and self.verbose:
            print('Priority-flood algorithm initialisation (%0.02f seconds)' % (clock() - t0))

        return

    def _performZhouAlgo(self):

        t0 = clock()
        fillZ = self.hLocal.getArray()
        self.fillLocal.setArray(fillZ)
        self.dm.localToGlobal(self.fillLocal, self.fillGlobal, 1)
        self.dm.globalToLocal(self.fillGlobal, self.fillLocal, 1)

        # Perform priority-flood on local domain
        watershed = np.zeros((fillZ.shape),dtype=int)
        fillZ = self.fillLocal.getArray()
        if self.first == 2:
            lcoords = self.lcoords
            lcoords[:,2] = fillZ
            self.zhouPit = fillAlgo.depressionFillingZhou(coords=lcoords, ngbIDs=self.FVmesh_ngbID,
                                                                                ngbNb=self.FVmesh_ngbNbs, meshIDs=self.inIDs,
                                                                                boundary=self.idLBounds, cartesian=False,
                                                                                sealevel=self.sealevel, extent=self.gbounds, first=self.first)
            del lcoords
        else:
            self.zhouPit = fillAlgo.depressionFillingZhou(Z=fillZ, cartesian=False,
                                                                                sealevel=self.sealevel, first=self.first)
        fillZ, locLabel, watershed, graph = self.zhouPit.performPitFillingUnstruct(simple=False)

        # Define globally unique watershed index
        label_offset = -np.ones(MPIsize+1, dtype=int)
        label_offset[MPIrank+1] = max(1,len(graph))+1
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, label_offset, op=MPI.MAX)
        offset = np.cumsum(label_offset)
        watershed += offset[MPIrank]
        locLabel[locLabel>0] += offset[-1]
        graph[:,0] += offset[MPIrank]
        ids = np.where(graph[:,1]>0)[0]
        graph[ids,1] += offset[MPIrank]

        # Transfer watershed values along local borders
        self.watershedLocal.setArray(watershed.astype(int))
        self.dm.localToGlobal(self.watershedLocal, self.watershedGlobal, 1)
        self.dm.globalToLocal(self.watershedGlobal, self.watershedLocal, 1)
        watershed = self.watershedLocal.getArray()

        # Transfer filled values along the local borders
        self.fillLocal.setArray(fillZ)
        self.dm.localToGlobal(self.fillLocal, self.fillGlobal, 1)
        self.dm.globalToLocal(self.fillGlobal, self.fillLocal, 1)
        fillZ = self.fillLocal.getArray()
        if self.first == 2:
            cgraph = self.zhouPit.combineUnstructGrids(fillZ, watershed, self.idLocalComm, self.idComm)
        else:
            cgraph = self.zhouPit.combineUnstructGrids(fillZ, watershed, None, self.idLocalComm)
        self.first = 0
        cgraph = np.concatenate((graph,cgraph))

        # Define global spillover graph on master
        label_offset = np.zeros(MPIsize+1, dtype=int)
        label_offset[MPIrank+1] = max(1,len(cgraph))
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, label_offset, op=MPI.MAX)
        offset = np.cumsum(label_offset)

        graph = -np.ones((np.sum(label_offset),3),dtype=float)
        graph[offset[MPIrank]:offset[MPIrank]+len(cgraph),:] =  cgraph
        if MPIrank == 0:
            mgraph = -np.ones((np.sum(label_offset),3),dtype=float)
        else:
            mgraph = None
        MPI.COMM_WORLD.Reduce(graph, mgraph, op=MPI.MAX, root=0)

        if MPIrank == 0:
            # Build bidrectional edges connections
            cgraph = pd.DataFrame(mgraph,columns=['source','target','weight'])
            cgraph = cgraph.sort_values('weight')
            cgraph = cgraph.drop_duplicates(['source', 'target'], keep='first')
            c12 = np.concatenate((cgraph['source'].values,cgraph['target'].values))
            cmax = np.max(np.bincount(c12.astype(int)))
            # Applying Barnes priority-flood algorithm on the bidirectional graph
            graph = self.zhouPit.fillGraph(cgraph.values,cmax)
        else:
            graph = None

        # Send filled graph dataset to each processors and perform pit filling
        graph = MPI.COMM_WORLD.bcast(graph, root=0)
        nn = locLabel.max() + int(watershed.max()) + 2
        self.z0 = self.hLocal.getArray()
        fillZ, self.pitIDs = fillDepression(self.z0 , fillZ, locLabel, watershed.astype(int), graph, nn)
        self.fillLocal.setArray(fillZ)
        self.dm.localToGlobal(self.fillLocal, self.fillGlobal, 1)
        self.dm.globalToLocal(self.fillGlobal, self.fillLocal, 1)
        del graph, cgraph, locLabel, watershed
        del fillZ, offset, label_offset, mgraph

        if MPIrank == 0 and self.verbose:
            print('Pit filling algorithm (%0.02f seconds)'% (clock() - t0))

        return

    def _definePitParams(self):

        t0 = clock()

        fillZ = self.fillLocal.getArray()
        label_offset = np.zeros(MPIsize+1, dtype=int)
        label_offset[MPIrank+1] = max(1,self.pitIDs.max())+2
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, label_offset, op=MPI.MAX)
        offset = np.cumsum(label_offset)
        self.pitIDs[self.pitIDs>0] += offset[MPIrank]
        pitVols = np.zeros(offset[-1]+1)
        equal = False

        while(not equal):
            self.pitLocal.setArray(self.pitIDs.astype(int))
            self.dm.localToGlobal(self.pitLocal, self.pitGlobal)
            self.dm.globalToLocal(self.pitGlobal, self.pitLocal)
            self.pitIDs = self.pitLocal.getArray().astype(int)

            self.pitIDs, newVols, spillPts = self.zhouPit.getPitData_unst(self.z0, fillZ, self.FVmesh_area,
                                                                                self.pitIDs, int(offset[-1]+1))
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, newVols, op=MPI.SUM)

            spillComb = None
            if MPIrank == 0:
                spillComb = np.empty([(offset[-1]+1)*MPIsize*2])
            MPI.COMM_WORLD.Gather(spillPts.flatten(), spillComb, root=0)

            if MPIrank == 0:
                spillPts = np.zeros(((offset[-1]+1)*MPIsize,5))
                spillPts[:,:2] = spillComb.reshape(((offset[-1]+1)*MPIsize,2))
                pitN = np.arange(1,offset[-1]+2)
                for k in range(MPIsize):
                    p = k*(offset[-1]+1)
                    spillPts[p:p+offset[-1]+1,2] = pitN
                    spillPts[p:p+offset[-1]+1,3] = k
                    spillPts[p:p+offset[-1]+1,4] = newVols
                spill = pd.DataFrame(spillPts,columns=['nid','z','pit','proc','vol'])
                spill = spill[spill['nid']>-1]
                spill = spill.sort_values(['pit', 'z'], ascending=[True, True])
                spill = spill.drop_duplicates(['pit'], keep='first').values
                spill[:,0] -= 1
            else:
                spill = None

            self.pitData = MPI.COMM_WORLD.bcast(spill, root=0)
            equal = np.array_equal(pitVols,newVols)
            pitVols = newVols.copy()

        self.fillLocal.setArray(fillZ)
        self.dm.localToGlobal(self.fillLocal, self.fillGlobal, 1)
        self.dm.globalToLocal(self.fillGlobal, self.fillLocal, 1)

        self.pitLocal.setArray(self.pitIDs)
        self.dm.localToGlobal(self.pitLocal, self.pitGlobal, 1)
        self.dm.globalToLocal(self.pitGlobal, self.pitLocal, 1)

        if MPIrank == 0 and self.verbose:
            print('Pit parameters definition (%0.02f seconds)'% (clock() - t0))

        return

    def depressionDefinition(self):

        self._performZhouAlgo()

        self._definePitParams()

        return

    def depositDepression(self):

        if len(self.pitData) == 0:
            self.pitData = np.zeros((1,5))
            self.pitData[0,2] = -1

        elev = self.hLocal.getArray()
        fillZ = self.fillLocal.getArray()
        cumed = self.cumEDLocal.getArray()

        # Find the deposited volume in each depression
        vol = self.vSedLocal.getArray()
        pitDep = np.zeros(self.npoints)
        pitDep[self.pitID] = vol[self.pitID]

        # Remove deposits on the domain boundary
        pitDep[self.idGBounds] = 0.

        # Compute cumulative deposition on each depression globally
        depLocal = np.zeros(self.npoints)
        depLocal[self.idLocal] = pitDep[self.idLocal]
        pitID = self.pitLocal.getArray()
        pitID[elev<=self.sealevel] = -1
        nb = max(int(self.pitData[:,2].max())+1,1)
        pitsedVol,diffDepLocal = pitVolume(depLocal,pitID,nb)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, pitsedVol, op=MPI.SUM)
        pitVol = np.zeros(len(pitsedVol))
        if self.pitData[:,2].max() > 0:
            pitVol[self.pitData[:,2].astype(int)-1] = self.pitData[:,4]

        # Get the percentage that will be deposited in each depression
        newZ,remainSed = pitHeight(elev,fillZ,pitID,pitVol,pitsedVol)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, remainSed, op=MPI.SUM)

        # Update elevation and cumulative erosion/deposition
        cumed += newZ-elev
        self.hLocal.setArray(newZ)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        self.cumEDLocal.setArray(cumed.reshape(-1))
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)

        # Distribute remaining sediment that will need to be diffused
        distID = np.where(remainSed>0.)[0]
        spillPts = spillPoints(len(distID),remainSed,self.pitData)
        distID = np.where(np.logical_and(spillPts[:,-1]>0.,spillPts[:,1]==MPIrank))[0]
        spillPts = spillPts[distID,:]

        # Add spill over points sediment volumes
        diffDepLocal[spillPts[:,0].astype(int)] += spillPts[:,-1]
        self.diffDepLocal.setArray(diffDepLocal)
        self.dm.localToGlobal(self.diffDepLocal, self.diffDep, 1)
        self.dm.globalToLocal(self.diffDep, self.diffDepLocal, 1)

        return
