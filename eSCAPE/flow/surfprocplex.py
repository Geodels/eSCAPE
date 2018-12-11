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
from scipy import sparse
# from scipy.spatial import cKDTree

import sys,petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from time import clock
import warnings;warnings.simplefilter('ignore')

from eSCAPE._fortran import setHillslopeCoeff
from eSCAPE._fortran import setDiffusionCoeff
from eSCAPE._fortran import MFDreceivers
from eSCAPE._fortran import singlePit
from eSCAPE._fortran import minHeight
from eSCAPE._fortran import diffusionDT
from eSCAPE._fortran import distributeHeight
from eSCAPE._fortran import distributeVolume

MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = PETSc.COMM_WORLD

try: range = xrange
except: pass

class SPMesh(object):
    """
    Building the surface processes based on different neighbour conditions
    """
    def __init__(self, *args, **kwargs):

        # KSP solver parameters
        self.rtol = 1.0e-8

        data = np.zeros(1,dtype=int)
        data[0] = self.FVmesh_ngbNbs.max()
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, data, op=MPI.MAX)
        self.maxNgbhs = data[0]

        # Diffusion matrix construction
        if self.Cd > 0.:
            diffCoeffs = setHillslopeCoeff(self.npoints,self.Cd*self.dt)
            self.Diff = self._matrix_build_diag(diffCoeffs[:,0])

            for k in range(0, self.maxNgbhs):
                tmpMat = self._matrix_build()
                indptr = np.arange(0, self.npoints+1, dtype=PETSc.IntType)
                indices = self.FVmesh_ngbID[:,k].copy()
                data = np.zeros(self.npoints)
                ids = np.nonzero(indices<0)
                indices[ids] = ids
                data = diffCoeffs[:,k+1]
                ids = np.nonzero(data==0.)
                indices[ids] = ids
                tmpMat.assemblyBegin()
                tmpMat.setValuesLocalCSR(indptr, indices.astype(PETSc.IntType), data,
                                                                PETSc.InsertMode.INSERT_VALUES)
                tmpMat.assemblyEnd()
                self.Diff += tmpMat
                tmpMat.destroy()
            del ids, indices, indptr

        # Identity matrix construction
        self.iMat = self._matrix_build_diag(np.ones(self.npoints))

        # Petsc vectors
        self.vecG = self.hGlobal.duplicate()
        self.vecL = self.hLocal.duplicate()
        self.vGlob = self.hGlobal.duplicate()
        self.vLoc = self.hLocal.duplicate()
        self.seaG = self.hGlobal.duplicate()
        self.seaL = self.hLocal.duplicate()
        self.tmpG = self.hGlobal.duplicate()
        self.tmpL = self.hLocal.duplicate()
        del data

        return

    def _matrix_build(self, nnz=(1,1)):
        """
        Define PETSC Matrix
        """
        matrix = PETSc.Mat().create(comm=MPIcomm)
        matrix.setType('aij')
        matrix.setSizes(self.sizes)
        matrix.setLGMap(self.lgmap_row, self.lgmap_col)
        matrix.setFromOptions()
        matrix.setPreallocationNNZ(nnz)

        return matrix

    def _matrix_build_diag(self, V, nnz=(1,1)):
        """
        Define PETSC Diagonal Matrix
        """

        matrix = self._matrix_build()

        # Define diagonal matrix
        I = np.arange(0, self.npoints+1, dtype=PETSc.IntType)
        J = np.arange(0, self.npoints, dtype=PETSc.IntType)
        matrix.assemblyBegin()
        matrix.setValuesLocalCSR(I, J, V, PETSc.InsertMode.INSERT_VALUES)
        matrix.assemblyEnd()

        return matrix

    def _make_reasons(self,reasons):
        return dict([(getattr(reasons, r), r)
             for r in dir(reasons) if not r.startswith('_')])

    def _solve_KSP(self, guess, matrix, vector1, vector2):
        """
        Set PETSC KSP solver.

        Args
            guess: Boolean specifying if the iterative KSP solver initial guess is nonzero
            matrix: PETSC matrix used by the KSP solver
            vector1: PETSC vector corresponding to the initial values
            vector2: PETSC vector corresponding to the new values

        Returns:
            vector2: PETSC vector of the new values
        """

        # guess = False
        ksp = PETSc.KSP().create(PETSc.COMM_WORLD)
        if guess:
            ksp.setInitialGuessNonzero(guess)
        ksp.setOperators(matrix,matrix)
        ksp.setType('richardson')
        pc = ksp.getPC()
        pc.setType('bjacobi')
        ksp.setTolerances(rtol=self.rtol)
        ksp.solve(vector1, vector2)
        r = ksp.getConvergedReason()
        if r < 0:
            KSPReasons = self._make_reasons(PETSc.KSP.ConvergedReason())
            print('LinearSolver failed to converge after %d iterations',ksp.getIterationNumber())
            print('with reason: %s',KSPReasons[r])
            raise RuntimeError("LinearSolver failed to converge!")
        ksp.destroy()

        return vector2

    def _buildFlowDirection(self):
        """
        Build multiple flow direction based on neighbouring slopes.
        """

        t0 = clock()

        # Attempt to remove single pit...
        hArrayLocal = self.hLocal.getArray()-self.sealevel
        h, shore = singlePit(self.inIDs, self.gbounds, hArrayLocal)
        self.hLocal.setArray(h+self.sealevel)
        self.shoreID = np.where(shore==1)[0]
        del h, shore
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        # Account for marine regions
        hArrayLocal = self.hLocal.getArray().copy()
        hArrayLocal -= self.sealevel
        self.seaID = np.where(hArrayLocal<=0.)[0]

        # Define shoreline
        # mask = shore>0
        # shorePtsXY = self.lcoords[mask,:2]
        # self.seaDep = np.zeros(hArrayLocal.shape)
        # pts_offset = np.zeros(MPIsize+1, dtype=int)
        # pts_offset[MPIrank+1] = max(1,len(shorePtsXY))
        # MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, pts_offset, op=MPI.MAX)
        # offset = np.cumsum(pts_offset)
        # ptsXY = np.full((np.sum(pts_offset),2),-1.e15)
        # ptsXY[offset[MPIrank]:offset[MPIrank]+len(shorePtsXY),:] =  shorePtsXY
        # MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, ptsXY, op=MPI.MAX)
        # _, uniqueX = np.unique(ptsXY[:,0].round(decimals=6), return_index=True)
        # _, uniqueY = np.unique(ptsXY[:,0].round(decimals=6), return_index=True)
        # unique = np.intersect1d(uniqueX,uniqueY)
        # shore = ptsXY[unique]
        # shore = shore[~np.all(shore == -1.e15, axis=1)]
        # del pts_offset, offset, mask
        # del ptsXY, shorePtsXY, unique, uniqueX, uniqueY

        # Get distance to shoreline for marine vertex
        # if len(shore)>0:
        #     tree = cKDTree(shore)
        #     seaPts = self.lcoords[self.seaID,:2]
        #     dist, _ = tree.query(seaPts,k=1)
        #     distShore = np.empty(hArrayLocal.shape)
        #     if len(dist)>0:
        #         distShore.fill(np.max(dist)+1.e4)
        #         distShore[self.seaID] = dist
        #         self.seaSort = np.argsort(distShore)
        #         mask = hArrayLocal[self.seaSort]>0
        #         self.seaSort[mask] = -1
        #         del mask
        #     del seaPts, dist, tree, distShore
        # else:
        #     self.seaSort = -np.ones(self.npoints)
        # del shore

        # Define multiple flow directions
        self.rcvID, self.slpRcv, self.distRcv, self.wghtVal = MFDreceivers(self.flowDir, self.inIDs, hArrayLocal)
        del hArrayLocal

        # Account for pit regions
        self.pitID = np.where(self.slpRcv[:,0]<=0.)[0]
        self.rcvID[self.pitID,:] = np.tile(self.pitID, (self.flowDir,1)).T
        self.distRcv[self.pitID,:] = 0.
        self.wghtVal[self.pitID,:] = 0.

        # Account for shore vertices
        # self.rcvID[self.shoreID,:] = np.tile(self.shoreID, (self.flowDir,1)).T
        # self.distRcv[self.shoreID,:] = 0.
        # self.wghtVal[self.shoreID,:] = 0.

    	if MPIrank == 0 and self.verbose:
            print('Flow Direction declaration (%0.02f seconds)'% (clock() - t0))

        return

    def FlowAccumulation(self):
        """
        Compute multiple flow accumulation.
        """

        self._buildFlowDirection()

        t0 = clock()
        # Build drainage area matrix
        if self.rainFlag:
            WAMat = self.iMat.copy()
            WeightMat = self._matrix_build_diag(np.zeros(self.npoints))
            for k in range(0, self.flowDir):
                tmpMat = self._matrix_build()
                data = -self.wghtVal[:,k].copy()
                indptr = np.arange(0, self.npoints+1, dtype=PETSc.IntType)
                nodes = indptr[:-1]
                data[self.rcvID[:,k].astype(PETSc.IntType)==nodes] = 0.0
                tmpMat.assemblyBegin()
                tmpMat.setValuesLocalCSR(indptr, self.rcvID[:,k].astype(PETSc.IntType),
                                         data, PETSc.InsertMode.INSERT_VALUES)
                tmpMat.assemblyEnd()
                WeightMat += tmpMat
                WAMat += tmpMat
                tmpMat.destroy()

            # Solve flow accumulation
            WAtrans = WAMat.transpose()
            self.WeightMat = WAtrans.copy()
            if self.tNow == self.tStart:
                self._solve_KSP(False, WAtrans, self.bG, self.drainArea)
            else:
                self._solve_KSP(True, WAtrans, self.bG, self.drainArea)
            WAMat.destroy()
            WAtrans.destroy()
            WeightMat.destroy()
        else:
            self.drainArea.set(0.)
        self.dm.globalToLocal(self.drainArea, self.drainAreaLocal, 1)

        if MPIrank == 0 and self.verbose:
            print('Compute Flow Accumulation (%0.02f seconds)'% (clock() - t0))

        return

    def _getErosionRate(self,Hsoil):
        """
        Compute sediment and bedrock erosion rates.
        """

        Kcoeff = self.drainAreaLocal.getArray()
        Kbr = np.power(Kcoeff,self.mbr)*self.Kbr*self.dt
        Kbr[self.seaID] = 0.
        Ksed = np.power(Kcoeff,self.msed)*self.Ksed*self.dt
        Ksed[self.seaID] = 0.

        # Initialise system matrices...
        if self.Kbr > 0.:
            EbedMat = self.iMat.copy()
        if self.Ksed > 0. and self.frac_fine < 1.:
            EsedMat = self.iMat.copy()
        wght = self.wghtVal.copy()

        # Define erosion coefficients
        for k in range(0, self.flowDir):

            indptr = np.arange(0, self.npoints+1, dtype=PETSc.IntType)
            nodes = indptr[:-1]
            # Define erosion limiter to prevent formation of flat
            dh = self.hOldArray-self.hOldArray[self.rcvID[:,k]]
            limiter = np.divide(dh, dh+1.e-3, out=np.zeros_like(dh), where=dh!=0)

            # Bedrock erosion processes SPL computation (maximum bedrock incision)
            if self.Kbr > 0.:
                data = np.divide(Kbr*limiter, self.distRcv[:,k], out=np.zeros_like(Kcoeff),
                                            where=self.distRcv[:,k]!=0)
                tmpMat = self._matrix_build()
                wght[self.seaID,k] = 0.
                data = np.multiply(data,-wght[:,k])
                data[self.rcvID[:,k].astype(PETSc.IntType)==nodes] = 0.0
                tmpMat.assemblyBegin()
                tmpMat.setValuesLocalCSR(indptr, self.rcvID[:,k].astype(PETSc.IntType),
                                         data, PETSc.InsertMode.INSERT_VALUES)
                tmpMat.assemblyEnd()
                EbedMat += tmpMat
                Mdiag = self._matrix_build_diag(data)
                EbedMat -= Mdiag
                tmpMat.destroy()
                Mdiag.destroy()
                del data

            # Sediment erosion processes SPL computation (maximum alluvial sediment incision)
            if self.Ksed > 0. and self.frac_fine < 1.:
                data = np.divide(Ksed*limiter, self.distRcv[:,k], out=np.zeros_like(Kcoeff),
                                            where=self.distRcv[:,k]!=0)
                tmpMat = self._matrix_build()
                if self.Kbr == 0.:
                    wght[self.seaID,k] = 0.
                data = np.multiply(data,-wght[:,k])
                data[self.rcvID[:,k].astype(PETSc.IntType)==nodes] = 0.0
                tmpMat.assemblyBegin()
                tmpMat.setValuesLocalCSR(indptr, self.rcvID[:,k].astype(PETSc.IntType),
                                         data, PETSc.InsertMode.INSERT_VALUES)
                tmpMat.assemblyEnd()
                EsedMat += tmpMat
                Mdiag = self._matrix_build_diag(data)
                EsedMat -= Mdiag
                tmpMat.destroy()
                Mdiag.destroy()
                del data

        if self.flowDir > 0:
            del dh, limiter, wght

        # Solve bedrock erosion thickness
        if self.Kbr > 0.:
            self._solve_KSP(True, EbedMat, self.hOld, self.vGlob)
            EbedMat.destroy()
            self.stepED.waxpy(-1.0,self.hOld,self.vGlob)
            # Define erosion rate (positive for incision)
            E = -self.stepED.getArray().copy()
            E = np.divide(E,self.dt)
            Ecrit = np.divide(E, self.crit_br, out=np.zeros_like(E),
                                   where=self.crit_br!=0)
            E -= self.crit_br*(1.0-np.exp(-Ecrit))
            E[E<0.] = 0.
            # We use soil thickness from previous time step !
            E = np.multiply(E,np.exp(-Hsoil/self.Hstar))
            self.Eb.setArray(E)
            self.dm.globalToLocal(self.Eb, self.EbLocal, 1)
            E = self.EbLocal.getArray().copy()
            E[self.seaID] = 0.
            self.EbLocal.setArray(E)
            del E, Ecrit
        else:
            self.EbLocal.set(0.)
        self.dm.localToGlobal(self.EbLocal, self.Eb, 1)

        # Solve sediment erosion rate
        if self.Ksed > 0. and self.frac_fine < 1.:
            self._solve_KSP(True, EsedMat, self.hOld, self.vGlob)
            EsedMat.destroy()
            self.stepED.waxpy(-1.0,self.hOld,self.vGlob)
            # Define erosion rate (positive for incision)
            E = -self.stepED.getArray().copy()
            E = np.divide(E,self.dt)
            Ecrit = np.divide(E, self.crit_sed, out=np.zeros_like(E),
                                   where=self.crit_sed!=0)
            E -= self.crit_sed*(1.0-np.exp(-Ecrit))
            E[E<0.] = 0.
            # Limit sediment erosion based on soil thickness
            E = np.multiply(E,self.dt*(1.0-np.exp(-Hsoil/self.Hstar)))
            ids = E>Hsoil
            E[ids] = Hsoil[ids]
            E = np.divide(E,self.dt)
            self.Es.setArray(E)
            self.dm.globalToLocal(self.Es, self.EsLocal, 1)
            E = self.EsLocal.getArray().copy()
            E[self.seaID] = 0.
            self.EsLocal.setArray(E)
            del ids, E, Ecrit
        else:
            self.EsLocal.set(0.)
        self.dm.localToGlobal(self.EsLocal, self.Es, 1)

        del Kcoeff, Kbr, Ksed

        return

    def _getSedFlux(self, Ds=None, Qw=None):
        """
        Compute sediment flux.
        """

        Qw = self.drainAreaLocal.getArray().copy()
        Qs = self.vland*self.FVmesh_area
        Qs[self.seaID] = 0.
        Kcoeff = np.divide(Qs, Qw, out=np.zeros_like(Qw), where=Qw!=0)
        SLMat = self._matrix_build_diag(Kcoeff)
        SLMat += self.WeightMat
        del Kcoeff, Qs

        # Define combined volumic erosion rate accounting for
        # fraction of fine and porosity
        Eb = self.Eb.getArray().copy()
        Eb = np.multiply(Eb,1.0-self.frac_fine)
        Es = self.Es.getArray().copy()
        Es = np.multiply(Es,1.0-self.phi)
        self.stepED.setArray(Es+Eb)
        self.stepED.pointwiseMult(self.stepED,self.areaGlobal)

        if self.tNow == self.tStart:
            self._solve_KSP(False, SLMat, self.stepED, self.vSed)
        else :
            self._solve_KSP(True, SLMat, self.stepED, self.vSed)
        SLMat.destroy()

        # Update local vector
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)
        self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)

        del Es, Eb

        return Qw

    def _getMaxDepositionRate(self, Qw=None, type=0):
        """
        Limit sediment deposition rate based on updated donor elevation.
        """

        # Transported sediment volume to deposit during time step
        if type == 0:
            tmpQ = self.vSedLocal.getArray().copy()
            tmpQ[tmpQ<0] = 0.
            Qs = self.vland*tmpQ
            depo = np.divide(Qs, Qw, out=np.zeros_like(Qw), where=Qw!=0)
            del Qs
        else:
            tmpQ = self.tmpL.getArray().copy()
            tmpQ[tmpQ<0] = 0.
            depo = np.zeros(tmpQ.shape)

        # Condition for pits and marine vertices
        depo[self.pitID] = np.divide(tmpQ[self.pitID], self.FVmesh_area[self.pitID],
                                     out=np.zeros_like(tmpQ[self.pitID]),
                                     where=self.FVmesh_area[self.pitID]!=0)
        depo[self.seaID] = np.divide(tmpQ[self.seaID], self.FVmesh_area[self.seaID],
                                     out=np.zeros_like(tmpQ[self.seaID]),
                                     where=self.FVmesh_area[self.seaID]!=0)
        depo[self.idGBounds] = 0.
        depo[depo<0.] = 0.
        del tmpQ

        # Define deposition thickness...
        depo = np.multiply(depo,self.dt)
        depo = np.divide(depo,1.0-self.phi)

        dhmax = np.zeros(depo.shape)
        if type == 0 and self.vland > 0.:
            # Erosion thicknesses during time step
            Eb = self.EbLocal.getArray().copy()
            Es = self.EsLocal.getArray().copy()
            # Limit deposition if required
            hArrayLocal = self.hLocal.getArray().copy()
            hArrayLocal -= np.multiply(Eb+Es,self.dt)
            dhmax = minHeight(self.inIDs, hArrayLocal)
            dhmax = np.multiply(0.9,dhmax)
            dhmax[dhmax<0] = 0.
            dhmax[self.idGBounds] = 0.
            del Eb, Es, hArrayLocal

        # Limit the maximum aerial deposition
        Ds = np.minimum(dhmax,depo)
        Ds[Ds<0] = 0.
        Ds[self.seaID] = depo[self.seaID]
        Ds[self.pitID] = depo[self.pitID]
        Ds = np.multiply(Ds,1.0-self.phi)
        Ds[self.idGBounds] = 0.
        del dhmax, depo

        # Define the deposition rate
        Ds = np.divide(Ds,self.dt)

        return Ds

    def _getDepositionVolume(self, Ds, force=False):
        """
        Compute sediment deposition volume.
        """

        # From deposition rate to volume
        Ds = np.multiply(Ds,self.FVmesh_area)
        tmpQ = np.divide(Ds*self.dt,1.0-self.phi)
        del Ds

        # Define excess sediment volume to deposit in depressions
        excessDep = np.zeros(self.npoints)
        excessDep[self.pitID] = tmpQ[self.pitID]
        excessDep[self.seaID] = 0.
        excessDep[self.idGBounds] = 0.
        excessDep[excessDep<0.] = 0.
        directAdd = np.zeros(self.npoints)
        if excessDep.max()<1. or force:
            directAdd = excessDep.copy()
            excessDep = np.zeros(self.npoints)
        self.vLoc.setArray(excessDep)
        self.dm.localToGlobal(self.vLoc, self.vGlob, 1)
        self.dm.globalToLocal(self.vGlob, self.vLoc, 1)

        # Update land deposition thicknesses
        tmpQ -= excessDep
        tmpQ[self.idGBounds] = 0.
        tmpQ += directAdd
        depo = np.divide(tmpQ, self.FVmesh_area, out=np.zeros_like(tmpQ),
                         where=self.FVmesh_area!=0)
        del tmpQ, excessDep, directAdd

        # Define marine deposition that will need to be diffused...
        diff = np.zeros(self.npoints)
        diff[self.seaID] = depo[self.seaID]
        depo[self.seaID] = 0.
        self.tmpL.setArray(diff)
        self.seaL.axpy(1.,self.tmpL)
        self.dm.localToGlobal(self.seaL, self.seaG, 1)
        self.dm.globalToLocal(self.seaG, self.seaL, 1)
        del diff

        # Update land deposition thickness
        self.vecL.setArray(depo)
        self.dm.localToGlobal(self.vecL, self.vecG, 1)
        del depo

        return

    def StreamPowerLaw(self):
        """
        Perform stream power law.
        """

        t0 = clock()

        # Constant local & global vectors/arrays
        self.Es.set(0.)
        self.Eb.set(0.)
        self.seaL.set(0.)
        self.seaG.set(0.)
        tmp = self.Eb.getArray()
        shape = tmp.shape
        del tmp
        self.seaL.set(0.)

        self.diffDep.set(0.)
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal, 1)
        self.hOldArray = self.hOldLocal.getArray().copy()

        # Get erosion values for considered time step
        Hsoil = self.Hsoil.getArray().copy()
        if self.rainFlag:
            self._getErosionRate(Hsoil)

        # Get sediment flux rate
        if self.rainFlag and self.frac_fine < 1.:
            Qw = self._getSedFlux()
        else:
            self.vSed.set(0.)
            self.vSedLocal.set(0.)

        # Get deposition volume
        if self.rainFlag and self.frac_fine < 1.:
            Ds = self._getMaxDepositionRate(Qw, type=0)
            self._getDepositionVolume(Ds, force=False)
            del Qw, Ds
        else:
            self.vGlob.set(0.) # Pit deposition
            self.vecG.set(0.) # Land deposition
            self.seaG.set(0.) # Marine diffusion

        # Update bedrock thicknesses due to erosion
        if  self.Kbr > 0.:
            Eb = self.Eb.getArray()
            Eb *= self.dt
            Eb[Eb<0] = 0.
        else:
            Eb = np.zeros(shape)

        # Update sediment thicknesses due to erosion
        if self.Ksed > 0. and self.frac_fine < 1.:
            Es = self.Es.getArray()
            Es *= self.dt
            ids = Es>Hsoil
            Es[ids] = Hsoil[ids]
            Es[Es<0] = 0.
            del ids
        else:
            Es = np.zeros(shape)
        del Hsoil

        # Update parameters
        self.stepED.setArray(-Es-Eb)
        self.stepED.axpy(1.,self.vecG)

        self.cumED.axpy(1.,self.stepED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        self.hGlobal.axpy(1.,self.stepED)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

        if self.Ksed > 0. and self.frac_fine < 1.:
            self.stepED.setArray(-Es)
            self.stepED.axpy(1.,self.vecG)
            self.Hsoil.axpy(1.,self.stepED)
            self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
            self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)

        del Es, Eb
        if MPIrank == 0 and self.verbose:
            print('Compute Stream Power Law (%0.02f seconds)'% (clock() - t0))

        return

    def downstreamDistribute(self, depoForce=False):
        """
        Distribute excess sediment volume downstream
        """

        t0 = clock()
        hOldArray = self.hLocal.getArray().copy()

        # Compute sediment fluxes for excess sediments
        self.diffDep.pointwiseMult(self.diffDep,self.areaGlobal)
        self._solve_KSP(False, self.WeightMat, self.diffDep, self.tmpG)
        self.dm.globalToLocal(self.tmpG, self.tmpL, 1)
        self.dm.localToGlobal(self.tmpL, self.tmpG, 1)

        # Get maximum possible deposition rate
        Ds = self._getMaxDepositionRate(None, type=1)

        # Get deposition volume
        self._getDepositionVolume(Ds, force=depoForce)
        del Ds

        self.cumED.axpy(1.,self.vecG)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        self.hGlobal.axpy(1.,self.vecG)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

        if self.Ksed > 0. and self.frac_fine < 1.:
            self.Hsoil.axpy(1.,self.vecG)
            self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
            self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)

        if MPIrank == 0 and self.verbose:
            print('Distribute excess sediment downstream (%0.02f seconds)'% (clock() - t0))

        return

    def _diffuse_TopSed(self, limit, elev, elev0, dh, dt):
        """
        Compute top sediment diffusion
        """

        # Diffusion matrix construction
        sedCoeffs = setDiffusionCoeff(self.sedimentK*dt, limit, elev, elev0, dh)

        sedDiff = self._matrix_build_diag(sedCoeffs[:,0])
        for k in range(0, self.maxNgbhs):
            tmpMat = self._matrix_build()
            indptr = np.arange(0, self.npoints+1, dtype=PETSc.IntType)
            indices = self.FVmesh_ngbID[:,k].copy()
            data = np.zeros(self.npoints)
            ids = np.nonzero(indices<0)
            indices[ids] = ids
            data = sedCoeffs[:,k+1]
            ids = np.nonzero(data==0.)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(indptr, indices.astype(PETSc.IntType), data,
                                     PETSc.InsertMode.INSERT_VALUES)
            tmpMat.assemblyEnd()
            sedDiff += tmpMat
            tmpMat.destroy()

        del sedCoeffs, data, ids, indices, indptr

        return sedDiff

    def _distributeSediment(self):
        """
        Distribute newly deposited sediments
        """

        t1 = clock()

        # Nothing to diffuse...
        if self.seaG.sum() <= 0. or self.maxIters==0:
            return

        self.hLocal.copy(result=self.hOldLocal)
        hL0 = self.hOldLocal.getArray().copy()

        # Get sediment volume in the marine environment
        self.tmpG.pointwiseMult(self.seaG,self.areaGlobal)
        self.dm.globalToLocal(self.tmpG, self.tmpL, 1)
        self.dm.localToGlobal(self.tmpL, self.tmpG, 1)
        vSed = self.tmpL.getArray().copy()
        hL = hL0.copy()
        # it = 0
        hmax = self.tmpG.max()[1]
        while hmax > 0.:
            nhL,nvSed = distributeHeight(self.inIDs, self.sealevel, hL, hL0, vSed)
            nhL[self.idGBounds] = hL0[self.idGBounds]
            nvSed[self.idGBounds] = 0.
            self.vecL.setArray(nhL)
            self.dm.localToGlobal(self.vecL, self.vecG, 1)
            self.dm.globalToLocal(self.vecG, self.vecL, 1)
            hL = self.vecL.getArray().copy()
            self.vecL.setArray(nvSed)
            self.dm.localToGlobal(self.vecL, self.vecG, 1)
            self.dm.globalToLocal(self.vecG, self.vecL, 1)
            vSed = self.vecL.getArray().copy()
            vSed[self.idGBounds] = 0.
            nvSed = distributeVolume(self.inIDs, self.sealevel, hL, hL0, vSed)
            nvSed[self.idGBounds] = 0.
            self.vecL.setArray(nvSed)
            self.dm.localToGlobal(self.vecL, self.vecG, 1)
            self.dm.globalToLocal(self.vecG, self.vecL, 1)
            vSed = self.vecL.getArray().copy()
            vSed[self.idGBounds] = 0.
            self.tmpG.pointwiseDivide(self.vecG,self.areaGlobal)
            hmax = self.tmpG.max()[1]
            # if hmax == 0.:
            #     break
            # it += 1

        self.hLocal.setArray(hL)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        limit = 1.e-6
        dt = self.dt/float(self.maxIters)
        step = int(self.dt/dt) + 1
        h0 = self.hOldLocal.getArray().copy()

        it = 0
        hL0 = self.hLocal.getArray().copy()
        stop = np.zeros(1)
        stop[0] = 0
        while(it<=step):
            elev = hL0
            self.hLocal.setArray(elev)
            self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
            self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
            dh = elev-h0
            dh[dh<0] = 0.
            sedDiff = self._diffuse_TopSed(limit, elev, h0-self.sealevel, dh, dt)
            self._solve_KSP(True, sedDiff, self.hGlobal, self.tmpG)
            self.dm.globalToLocal(self.tmpG, self.tmpL, 1)
            self.dm.localToGlobal(self.tmpL, self.tmpG, 1)
            self.tmpL.copy(result=self.hLocal)
            self.tmpG.copy(result=self.hGlobal)
            sedDiff.destroy()
            hL0 = self.hLocal.getArray().copy()
            stop[0] = 0
            bid = np.where(hL0>h0)[0]
            if len(bid)>0:
                if hL0[bid].max()>self.sealevel:
                    stop[0] = 1
                if hL0[bid].max()-self.sealevel < 1.e-3:
                    idd = np.where(hL0[bid]>self.sealevel)[0]
                    hL0[bid[idd]] = self.sealevel-1.e-3
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, stop, op=MPI.MAX)
            if stop[0] == 0 and it == step:
                break
            if stop[0] == 1 and it == step:
                it = step-1
            it += 1

        # Update elevation change vector
        self.tmpL.setArray(dh)
        self.dm.localToGlobal(self.tmpL, self.stepED, 1)
        self.dm.globalToLocal(self.stepED, self.tmpL, 1)
        del dh

        self.hLocal.waxpy(1.0,self.tmpL,self.hOldLocal)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        self.cumED.axpy(1.,self.stepED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)
        if self.Ksed > 0. and self.frac_fine < 1.:
            self.Hsoil.axpy(1.,self.stepED)
            Hsoil = self.Hsoil.getArray().copy()
            Hsoil[Hsoil<0.] = 0.
            self.Hsoil.setArray(Hsoil)
            self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
            del Hsoil

        if MPIrank == 0 and self.verbose:
            print('Deposited sediment distribution (%0.02f seconds)'% (clock() - t1))

        return

    def SedimentDiffusion(self):
        """
        Initialise sediment diffusion from pit and marine deposition.
        """

        t0 = clock()
        if self.sedimentK > 0:
            self._distributeSediment()
            if MPIrank == 0 and self.verbose:
                print('Compute Sediment Diffusion (%0.02f seconds)'% (clock() - t0))
        else:
            self.hLocal.axpy(1.,self.seaL)
            self.hGlobal.axpy(1.,self.seaG)
            self.cumED.axpy(1.,self.seaG)
            if self.Ksed > 0. and self.frac_fine < 1.:
                self.Hsoil.axpy(1.,self.seaG)
                Hsoil = self.Hsoil.getArray().copy()
                Hsoil[Hsoil<0.] = 0.
                self.Hsoil.setArray(Hsoil)
                self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
                del Hsoil

        return

    def HillSlope(self):
        """
        Perform hillslope diffusion.
        """

        t0 = clock()
        if self.Cd > 0.:
            # Get erosion values for considered time step
            self.hGlobal.copy(result=self.hOld)
            self._solve_KSP(True, self.Diff, self.hOld, self.hGlobal)

            # Update cumulative erosion/deposition and soil/bedrock elevation
            self.stepED.waxpy(-1.0,self.hOld,self.hGlobal)
            self.cumED.axpy(1.,self.stepED)
            if self.Ksed > 0. and self.frac_fine < 1.:
                self.Hsoil.axpy(1.,self.stepED)
                Hsoil = self.Hsoil.getArray().copy()
                Hsoil[Hsoil<0.] = 0.
                self.Hsoil.setArray(Hsoil)
                self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
                self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)
            self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
            self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)
            self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
            self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

        # Remove erosion/deposition on boundary nodes
        hArray = self.hLocal.getArray()
        hArray[self.idGBounds] = self.hOldArray[self.idGBounds]
        self.hLocal.setArray(hArray)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

        if MPIrank == 0 and self.verbose:
            print('Compute Hillslope Processes (%0.02f seconds)'% (clock() - t0))

        return
