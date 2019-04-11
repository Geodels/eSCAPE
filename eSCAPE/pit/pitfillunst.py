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
from mpi4py import MPI
import sys,petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from time import clock
import warnings;warnings.simplefilter('ignore')

import fillit as fillAlgo
from eSCAPE._fortran import pitData

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
        self.first = 1
        self.sealimit = 1.e4

        area = np.zeros(self.gpoints)
        area[self.natural2local] = self.FVmesh_area

        if MPIrank == 0:
            garea = np.zeros(self.gpoints)
        else:
            garea = None
        MPI.COMM_WORLD.Reduce(area, garea, op=MPI.MAX, root=0)
        del area

        if MPIrank == 0:
            self.eScapeGPit = fillAlgo.depressionFillingScape(ngbIDs=self.Gmesh_ngbID, ngbNb=self.Gmesh_ngbNbs,
                                                              boundary=self.boundGlob, first=-1, area=garea)
        else:
            self.eScapeGPit = None
        del garea

        self.fillGlobal = self.dm.createGlobalVector()
        self.fillLocal = self.dm.createLocalVector()

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
        self.idComm = np.unique(self.lcells[out].flatten()[ids])

        # Local points that are part of the global mesh boundary
        self.idGBounds = np.where(np.isin(self.idLocal,self.localboundIDs))[0]
        self.gbounds = np.zeros(self.npoints,dtype=int)
        self.gbounds[self.idGBounds] = 1

        if MPIrank == 0 and self.verbose:
            print('Priority-flood algorithm initialisation (%0.02f seconds)' % (clock() - t0))

        return

    def getDepressions(self):
        """
        Perform pit filling based on vertex elevation.
        """

        tot = np.zeros(1, dtype=int)
        tot[0] = len(self.seaID)+self.nbPit
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, tot, op=MPI.SUM)
        if tot[0] == 0 or self.frac_fine == 1.:
            self.shedID.set(-1.)
            self.pHeight = None
            self.pVol = None
            self.pitNode = None
            self.pitProc = None
            return

        t0 = clock()
        elev = np.zeros(self.gpoints)
        elev.fill(-1.e8)
        h = self.hLocal.getArray().copy()
        elev[self.natural2local] = h
        if MPIrank == 0:
            melev = np.zeros(self.gpoints)
        else:
            melev = None
        MPI.COMM_WORLD.Reduce(elev, melev, op=MPI.MAX, root=0)
        del elev
        self.sl_limit = self.sealevel-self.sealimit

        if MPIrank == 0:
            eps = 1.e-8
            seaIDs = np.where(melev<self.sl_limit)[0]
            fill,wshed,pitvol,pith,pitNode = self.eScapeGPit.performPitFillingEpsilon(melev,seaIDs,eps,type=1)
            tmp1 = np.concatenate((fill,wshed))
            tmp2 = np.concatenate((pitvol,pith,pitNode))
            nn = len(pitNode)
            del fill,wshed,pitvol,pith,pitNode
        else:
            tmp1 = None
            tmp2 = None
            nn = None
        del melev

        MPI.COMM_WORLD.Barrier()
        nb = MPI.COMM_WORLD.bcast(nn, root=0)
        gtmp2 = MPI.COMM_WORLD.bcast(tmp2, root=0)
        self.pVol = gtmp2[:nb]
        self.pHeight = gtmp2[nb:2*nb]
        pNode = gtmp2[2*nb:]

        # Make the for loop in fortran!!!!
        tmp3 = pitData(MPIrank,pNode,self.inIDs,self.natural2local)

        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, tmp3, op=MPI.MAX)
        self.pitProc = tmp3[:nb]
        self.pitNode = tmp3[nb:]
        gtmp1 = MPI.COMM_WORLD.bcast(tmp1, root=0)
        gfill = gtmp1[:self.gpoints]
        gwshed = gtmp1[self.gpoints:]

        fillLocal = gfill[self.natural2local]
        wshedLocal = gwshed[self.natural2local]
        del gfill, gwshed, h
        del tmp3, tmp2, tmp1, gtmp1, pNode

        self.fillLocal.setArray(fillLocal)
        del fillLocal
        self.dm.localToGlobal(self.fillLocal, self.fillGlobal, 1)
        self.dm.globalToLocal(self.fillGlobal, self.fillLocal, 1)

        self.shedIDLocal.setArray(wshedLocal)
        del wshedLocal

        self.dm.localToGlobal(self.shedIDLocal, self.shedID, 1)
        self.dm.globalToLocal(self.shedID, self.shedIDLocal, 1)

        if MPIrank == 0 and self.verbose:
            print('Fill Pit Depression (%0.02f seconds)'% (clock() - t0))

        return
