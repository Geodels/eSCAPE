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

import numpy as np
import pandas as pd
from mpi4py import MPI
import sys,petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from time import clock
import warnings;warnings.simplefilter('ignore')

import fillit as fillAlgo
from eSCAPE._fortran import fillDepression
from eSCAPE._fortran import combinePit
from eSCAPE._fortran import pitVolume
from eSCAPE._fortran import pitHeight
from eSCAPE._fortran import addExcess

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
        self.sealimit = 2000.

        self.fillGlobal = self.dm.createGlobalVector()
        self.fillLocal = self.dm.createLocalVector()

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
        self.commID = idComm

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

    def _defineDepressionParameters(self):
        """
        Perform the pit filling algorithm proposed by Zhou et al., 2016 -- (https://github.com/cageo/Zhou-2016)
        """

        t0 = clock()
        self.sl_limit = self.sealevel-self.sealimit
        fillZ = self.hLocal.getArray()+self.sealimit
        self.fillLocal.setArray(fillZ)
        self.dm.localToGlobal(self.fillLocal, self.fillGlobal, 1)
        self.dm.globalToLocal(self.fillGlobal, self.fillLocal, 1)
        elev = self.hLocal.getArray().copy()
        seaIDs = np.where(elev<self.sl_limit)

        # Perform priority-flood on local domain
        watershed = np.zeros((fillZ.shape),dtype=int)
        fillZ = self.fillLocal.getArray()

        if self.first == 2:
            lcoords = self.lcoords.copy()
            lcoords[:,2] = fillZ
            self.eScapePit = fillAlgo.depressionFillingScape(coords=lcoords, ngbIDs=self.FVmesh_ngbID,
                                                          ngbNb=self.FVmesh_ngbNbs, meshIDs=self.inIDs,
                                                          boundary=self.idLBounds, seaIDs=seaIDs,
                                                          extent=self.gbounds, first=self.first)
            del lcoords
        else:
            self.eScapePit = fillAlgo.depressionFillingScape(Z=fillZ, seaIDs=seaIDs, first=self.first)
        fillZ, watershed, graph = self.eScapePit.performPitFillingUnstruct(simple=False)

        del seaIDs

        # Define globally unique watershed index
        label_offset = -np.ones(MPIsize+1, dtype=int)
        label_offset[MPIrank+1] = max(1,len(graph))+1
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, label_offset, op=MPI.MAX)
        label_offset[0] = 0
        offset = np.cumsum(label_offset)
        watershed += offset[MPIrank]
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
            cgraph = self.eScapePit.combineUnstructGrids(fillZ, watershed, self.idLocalComm, self.idComm)
        else:
            cgraph = self.eScapePit.combineUnstructGrids(fillZ, watershed, None, self.idLocalComm)
        self.first = 0
        cgraph = np.concatenate((graph,cgraph))

        # Define global spillover graph on master
        label_offset = np.zeros(MPIsize+1, dtype=int)
        label_offset[MPIrank+1] = max(1,len(cgraph))
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, label_offset, op=MPI.MAX)
        offset = np.cumsum(label_offset)
        graph = -np.ones((np.sum(label_offset),5),dtype=float)
        graph[offset[MPIrank]:offset[MPIrank]+len(cgraph),:4] =  cgraph
        graph[offset[MPIrank]:offset[MPIrank]+len(cgraph),4] =  MPIrank
        if MPIrank == 0:
            mgraph = -np.ones((np.sum(label_offset),5),dtype=float)
        else:
            mgraph = None
        MPI.COMM_WORLD.Reduce(graph, mgraph, op=MPI.MAX, root=0)

        if MPIrank == 0:
            # Build bidrectional edges connections
            cgraph = pd.DataFrame(mgraph,columns=['source','target','weight','spill','rank'])
            cgraph = cgraph.sort_values('weight')
            cgraph = cgraph.drop_duplicates(['source', 'target'], keep='first')
            c12 = np.concatenate((cgraph['source'].values,cgraph['target'].values))
            cmax = np.max(np.bincount(c12.astype(int)))+1
            # Applying Barnes priority-flood algorithm on the bidirectional graph
            ggraph = self.eScapePit.fillGraph(cgraph.values,cmax)
            ggraph[:,0] -= self.sealimit
        else:
            ggraph = None

        # Send filled graph dataset to each processors and perform pit filling
        # self.graph array contains:
        # - pit filling elevation,
        # - spill over node,
        # - rank,
        # - pit to which the considered pit is spilling towards
        # - order of pit filling
        self.graph = MPI.COMM_WORLD.bcast(ggraph, root=0)

        # Drain pit on the boundary towards the edges
        keep = self.graph[:,2].astype(int)==MPIrank
        proc = -np.ones(len(self.graph))
        proc[keep] = self.graph[keep,1]
        keep = proc>-1
        proc[keep] = self.gbounds[proc[keep].astype(int)]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, proc, op=MPI.MAX)
        ids = np.where(proc==1)[0]
        ids2 = np.where(self.graph[ids,0]==self.graph[self.graph[ids,3].astype(int),0])[0]
        self.graph[ids[ids2],3] = 0.
        self.graph[self.graph[:,0] < -1.e8,0] = -1.e8
        self.graph[self.graph[:,0] > 1.e6,0] = -1.e8
        del ids2

        # Define global solution by combining depressions/flat together
        fillZ, self.pitIDs, pitVol, combPit = fillDepression(self.hLocal.getArray(), fillZ-self.sealimit,
                                                             watershed.astype(int), self.graph[:,0],
                                                             self.graph[:,3].astype(int), self.FVmesh_area,
                                                             self.inIDs)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, pitVol, op=MPI.MAX)

        gPit, self.pitIDG, gVol, gOver = combinePit(len(combPit), len(self.pitIDs), combPit, pitVol,
                                            self.pitIDs, self.graph[:,-1].astype(int))

        self.pitDef = -np.ones((len(gVol),5))
        ids = np.where(self.graph[:,2]>-1)[0]
        ids2 = np.where(gPit[ids]==-1)[0]
        # Combined pit IDs
        self.pitDef[ids,0] = gPit[ids]
        self.pitDef[ids[ids2],0] = ids[ids2]
        # Spill over point ID
        self.pitDef[ids,1] = self.graph[gOver[ids],1]
        self.pitDef[ids[ids2],1] = self.graph[ids[ids2],1]
        # Spill over point processor ID
        self.pitDef[ids,2] = self.graph[gOver[ids],2]
        self.pitDef[ids[ids2],2] = self.graph[ids[ids2],2]
        # Volume
        self.pitDef[ids,3] = gVol[ids]
        # Elevation
        self.pitDef[ids,4] = self.graph[ids,0]
        del ids, ids2, gPit, gVol, gOver

        self.fillLocal.setArray(fillZ)
        self.dm.localToGlobal(self.fillLocal, self.fillGlobal, 1)
        self.dm.globalToLocal(self.fillGlobal, self.fillLocal, 1)
        del graph, cgraph, watershed, ggraph, mgraph
        del fillZ, offset, label_offset

        if MPIrank == 0 and self.verbose:
            print('Define depression parameters (%0.02f seconds)'% (clock() - t0))

        return

    def _depositDepression(self):
        """
        Perform pit filling based on eroded sediment volume.
        """

        # Find the deposited volume in each depression
        pitDep = self.vLoc.getArray()

        elev = self.hLocal.getArray()
        fillZ = self.fillLocal.getArray()
        cumed = self.cumEDLocal.getArray()

        # Remove deposits on the domain boundary
        pitDep[self.idGBounds] = 0.

        # Compute cumulative deposition on each depression globally
        depLocal = np.zeros(self.npoints)
        hAdd = depLocal.copy()
        depLocal[self.idLocal] = pitDep[self.idLocal]

        nonzeroDep = np.where(depLocal>0)[0]
        definedPit = np.where(self.pitIDG>0)[0]
        isnt = np.isin(nonzeroDep,definedPit,invert=True)

        if np.sum(isnt) > 0:
            ids = nonzeroDep[isnt]
            nids = self.FVmesh_ngbID[ids,:]
            hAdd[nonzeroDep[isnt]] = np.divide(depLocal[nonzeroDep[isnt]], self.FVmesh_area[nonzeroDep[isnt]],
                                out=np.zeros_like(hAdd[nonzeroDep[isnt]]), where=self.FVmesh_area[nonzeroDep[isnt]]!=0)
            hAdd[nonzeroDep[isnt]] = 0.
            if hAdd.max() > 1.e-4:
                print('Warning: Some pits are not defined as being part of a depression!')
            depLocal[nonzeroDep[isnt]] = 0.
        del isnt,nonzeroDep,definedPit

        pitsedVol = pitVolume(depLocal,self.pitIDG,self.pitDef[:,0].max()+1)

        # Get the percentage that will be deposited in each depression
        # and the excess sediment to distribute if any...
        keep = self.pitDef[:,0].astype(int)>-1
        percDep = np.zeros(len(self.pitDef))
        percDep[keep] = pitsedVol[self.pitDef[keep,0].astype(int)]
        remainSed = percDep - self.pitDef[:,3]
        remainSed[remainSed<1.e-8] = 0.
        remainSed[self.pitDef[:,3]<0] = 0.

        percDep = np.divide(percDep, self.pitDef[:,3], out=np.zeros_like(percDep),
                            where=self.pitDef[:,3]!=0)
        percDep[percDep>1.] = 1.
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, percDep, op=MPI.MAX)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, remainSed, op=MPI.MAX)

        # Update depression volume
        self.pitDef[:,3] = np.multiply(self.pitDef[:,3],1.-percDep)

        # Update height from deposited sediment volume in each depression
        perc = np.zeros(len(pitsedVol))
        perc[self.pitDef[keep,0].astype(int)] = percDep[keep]
        newZ = pitHeight(elev, fillZ, self.pitIDG, perc)

        # Update elevation and cumulative erosion/deposition
        newZ += hAdd
        cumed += newZ-elev
        self.hLocal.setArray(newZ)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        self.cumEDLocal.setArray(cumed.reshape(-1))
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)

        # Move excess sediment volume to the downstream node of the overspill node...
        ids = np.where(remainSed>0.)[0]

        if len(ids) == 0:
            del pitDep, elev, fillZ, cumed, newZ, ids
            del perc, percDep, remainSed, keep
            if MPIrank == 0 and self.verbose:
                print('Fill Pit Depression (%0.02f seconds)'% (clock() - t0))
            return False

        downSed = np.zeros((len(ids),3))
        downSed[:,0] = self.pitDef[ids,2]
        downSed[:,1] = self.pitDef[ids,1]
        downSed[:,2] = remainSed[ids]
        keep = downSed[:,0].astype(int) == MPIrank
        if np.sum(keep)>0:
            depLocal = addExcess(newZ, downSed[keep,1:])
        else:
            depLocal = np.zeros(self.npoints)
        self.diffDepLocal.setArray(np.divide(depLocal,self.dt))
        self.dm.localToGlobal(self.diffDepLocal, self.diffDep, 1)
        self.dm.globalToLocal(self.diffDep, self.diffDepLocal, 1)
        del pitDep, elev, fillZ, cumed, newZ, ids, downSed
        del perc, percDep, remainSed, keep, depLocal, hAdd

        if MPIrank == 0 and self.verbose:
            print('Fill Pit Depression (%0.02f seconds)'% (clock() - t0))

        return True

    def computeDepression(self):
        """
        Successively fill depression until all deposits in depressions are distributed.
        """

        # Find the deposited volume in each depression
        if self.vGlob.max()[1] == 0.:
            return

        iters = 0
        self._defineDepressionParameters()
        keepFilling = self._depositDepression()
        while(keepFilling):
            self.downstreamDistribute(depoForce=False)
            self._defineDepressionParameters()
            keepFilling = self._depositDepression()
            iters += 1
            if iters>100:
                keepFilling = False

        if self.diffDep.max()>0.:
            self.downstreamDistribute(depoForce=True)

        return
