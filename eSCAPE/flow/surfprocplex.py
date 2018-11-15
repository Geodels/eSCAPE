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
from scipy.spatial import cKDTree

import sys,petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from time import clock
import warnings;warnings.simplefilter('ignore')

from eSCAPE._fortran import setHillslopeCoeff
from eSCAPE._fortran import initDiffCoeff
from eSCAPE._fortran import MFDreceivers
from eSCAPE._fortran import minHeight
from eSCAPE._fortran import diffusionDT

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

        # Identity matrix construction
        self.iMat = self._matrix_build_diag(np.ones(self.npoints))

        # Sediment pile diffusion coefficients
        dt = np.zeros(1,dtype=float)
        self.streamKd, self.oceanKd, dt[0] = initDiffCoeff(self.npoints,self.dt,self.streamCd,self.oceanCd)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, dt, op=MPI.MIN)
        self.diffDT = dt[0]
        self.streamKd =  self.streamKd*self.diffDT
        self.oceanKd =  self.oceanKd*self.diffDT
        iters = int(self.dt/self.diffDT)+1
        self.maxIters = max(self.minIters,iters)

        # Smoothing matrix construction
        smthcoeff = max(self.streamCd,self.oceanCd)
        diffCoeffs = setHillslopeCoeff(self.npoints,smthcoeff*self.diffDT)
        self.Smooth = self._matrix_build_diag(diffCoeffs[:,0])
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
            self.Smooth += tmpMat
            tmpMat.destroy()

        # Petsc vectors
        self.hG0 = self.hGlobal.duplicate()
        self.hL0 = self.hLocal.duplicate()
        self.vecG = self.hGlobal.duplicate()
        self.vecL = self.hLocal.duplicate()
        self.vGlob = self.hGlobal.duplicate()
        self.vLoc = self.hLocal.duplicate()

        del diffCoeffs, data, ids, indices, indptr

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

    def _buidFlowDirection(self):
        """
        Build multiple flow direction based on neighbouring slopes.
        """

        t0 = clock()
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        hArrayLocal = self.hLocal.getArray()-self.sealevel
        self.rcvID, shore, self.slpRcv, self.distRcv, self.wghtVal = MFDreceivers(self.flowDir, self.inIDs, hArrayLocal)
        # Account for pit regions
        pitID = np.where(self.slpRcv[:,0]<=0.)[0]
        self.pitID = pitID[pitID >= 0]
        self.rcvID[self.pitID,:] = np.tile(self.pitID, (self.flowDir,1)).T
        self.distRcv[self.pitID,:] = 0.
        self.wghtVal[self.pitID,:] = 0.
        # Account for marine regions
        self.seaID = np.where(hArrayLocal<=0)[0]

        # Define shoreline
        mask = shore>0
        shorePtsXY = self.lcoords[mask,:2]
        self.seaDep = np.zeros(hArrayLocal.shape)
        pts_offset = np.zeros(MPIsize+1, dtype=int)
        pts_offset[MPIrank+1] = max(1,len(shorePtsXY))
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, pts_offset, op=MPI.MAX)
        offset = np.cumsum(pts_offset)
        ptsXY = np.full((np.sum(pts_offset),2),-1.e15)
        ptsXY[offset[MPIrank]:offset[MPIrank]+len(shorePtsXY),:] =  shorePtsXY
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, ptsXY, op=MPI.MAX)
        _, uniqueX = np.unique(ptsXY[:,0].round(decimals=6), return_index=True)
        _, uniqueY = np.unique(ptsXY[:,0].round(decimals=6), return_index=True)
        unique = np.intersect1d(uniqueX,uniqueY)
        shore = ptsXY[unique]
        shore = shore[~np.all(shore == -1.e15, axis=1)]
        del pts_offset, offset, mask
        del ptsXY, shorePtsXY, unique, uniqueX, uniqueY

        # Get distance to shoreline for marine vertex

        if len(shore)>0:
            tree = cKDTree(shore)
            seaPts = self.lcoords[self.seaID,:2]
            dist, _ = tree.query(seaPts,k=1)
            distShore = np.zeros(hArrayLocal.shape)
            distShore[self.seaID] = dist
            del seaPts, dist, tree

            # Get maximum elevation based on delta slope 0.1% (0.1 m per 100 m)
            self.seaDep[self.seaID] = -1.e-4*distShore[self.seaID]-hArrayLocal[self.seaID]
            self.seaDep[self.seaDep<0] = 0.
            del distShore, hArrayLocal
        else:
            del shore,hArrayLocal

    	if MPIrank == 0 and self.verbose:
            print('Flow Direction declaration (%0.02f seconds)'% (clock() - t0))

        # import pandas as pd
        # hArrayLocal = self.hLocal.getArray()-self.sealevel
        # depomar = np.zeros((len(hArrayLocal),3))
        # depomar[:,:2] = self.lcoords[:,:2]
        # depomar[:,2] = self.seaDep
        # df = pd.DataFrame(depomar,columns=['x','y','z'])
        # df.to_csv('out.csv')
        #
        # grtfh
        return

    def FlowAccumulation(self):
        """
        Compute multiple flow accumulation.
        """

        self._buidFlowDirection()

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
                tmpMat.setValuesLocalCSR(indptr, self.rcvID[:,k].astype(PETSc.IntType), data,
                                                                PETSc.InsertMode.INSERT_VALUES)
                tmpMat.assemblyEnd()
                WeightMat += tmpMat
                WAMat += tmpMat
                tmpMat.destroy()

            # Solve flow accumulation
            WAtrans = WAMat.transpose()
            self.WeightMat = WAtrans.copy()
            self.dm.localToGlobal(self.bL, self.bG, 1)
            self.dm.globalToLocal(self.bG, self.bL, 1)
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
        EbedMat = self.iMat.copy()
        EsedMat = self.iMat.copy()
        wght = self.wghtVal.copy()

        # Define erosion coefficients
        for k in range(0, self.flowDir):

            indptr = np.arange(0, self.npoints+1, dtype=PETSc.IntType)
            nodes = indptr[:-1]

            # Define erosion limiter to prevent formation of flat
            dh = self.hOldArray-self.hOldArray[self.rcvID[:,k]]
            limiter = np.divide(dh, dh+1.e-6, out=np.zeros_like(dh), where=dh!=0)

            # Bedrock erosion processes SPL computation (maximum bedrock incision)
            if self.Kbr > 0.:
                data = np.divide(Kbr*limiter, self.distRcv[:,k], out=np.zeros_like(Kcoeff),
                                            where=self.distRcv[:,k]!=0)
                tmpMat = self._matrix_build()
                wght[self.seaID,k] = 0.
                data = np.multiply(data,-wght[:,k])

                data[self.rcvID[:,k].astype(PETSc.IntType)==nodes] = 0.0
                tmpMat.assemblyBegin()
                tmpMat.setValuesLocalCSR(indptr, self.rcvID[:,k].astype(PETSc.IntType), data,
                                                                PETSc.InsertMode.INSERT_VALUES)
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
                tmpMat.setValuesLocalCSR(indptr, self.rcvID[:,k].astype(PETSc.IntType), data,
                                                                PETSc.InsertMode.INSERT_VALUES)
                tmpMat.assemblyEnd()
                EsedMat += tmpMat
                Mdiag = self._matrix_build_diag(data)
                EsedMat -= Mdiag
                tmpMat.destroy()
                Mdiag.destroy()
                del data
            del dh,limiter

        # Solve bedrock erosion thickness
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
        # Erosion rate (semi-implicit as we account for Hsoil from previous time step)
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
            del ids
        else:
            self.EsLocal.set(0.)
        self.dm.localToGlobal(self.EsLocal, self.Es, 1)

        del E, Ecrit, wght
        del Kcoeff, Kbr, Ksed

        return

    def _getSedFlux(self, first, Ds=None, Qw=None):
        """
        Compute sediment flux.
        """

        if first:
            Qw = self.drainAreaLocal.getArray().copy()
            if self.vland == 0. and self.vsea == 0.:
                return Qw
            Qs = self.vland*self.FVmesh_area
            #Qs[self.seaID] = self.vsea*self.FVmesh_area[self.seaID]
            Kcoeff = np.divide(Qs, Qw, out=np.zeros_like(Qw), where=Qw!=0)
            Kcoeff[self.seaID] = 0.

            # Maximum volumic rate of marine deposition
            seaDmax = np.zeros(Qw.shape)
            hArrayLocal = self.hLocal.getArray().copy()
            seaDmax[self.seaID] = np.maximum(0.,self.seaDep[self.seaID])
            seaDmax = np.multiply(seaDmax, 1.0-self.phi)
            seaDmax = np.divide(seaDmax,self.dt)
            self.vSedLocal.setArray(seaDmax)
            self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)
            seaDmax = self.vSed.getArray().copy()
            SLMat = self._matrix_build_diag(Kcoeff)
            SLMat += self.WeightMat
            del Kcoeff, Qs, hArrayLocal
        else:
            if self.vland == 0. and self.vsea == 0.:
                return
            Qs = np.zeros(Ds.shape)
            Qs[self.deepSed] = self.vsea*self.FVmesh_area[self.deepSed]
            Kcoeff = np.divide(Qs, Qw, out=np.zeros_like(Qw), where=Qw!=0)
            SLMat = self._matrix_build_diag(Kcoeff)
            SLMat += self.WeightMat

        # Define combined volumic erosion rate accounting for fraction of fine and porosity
        Eb = self.Eb.getArray().copy()
        Eb = np.multiply(Eb,1.0-self.frac_fine)
        Es = self.Es.getArray().copy()
        Es = np.multiply(Es,1.0-self.phi)

        if first:
            self.stepED.setArray(Es+Eb-seaDmax)
            self.stepED.pointwiseMult(self.stepED,self.areaGlobal)
            del seaDmax
        else:
            self.stepED.setArray(Es+Eb)
            self.stepED.axpy(-1.,self.vSed) # Remove deposition rate
            self.stepED.pointwiseMult(self.stepED,self.areaGlobal)

        if self.tNow == self.tStart:
            self._solve_KSP(False, SLMat, self.stepED, self.vSed)
        else :
            self._solve_KSP(True, SLMat, self.stepED, self.vSed)

        SLMat.destroy()

        del Es, Eb

        if first:
            return Qw

        return

    def _getMaxDepositionRate(self, Qw):
        """
        Limit sediment deposition rate based on updated donor elevation.
        """

        # First get maximum deposition during time step
        tmpQ = self.vSedLocal.getArray().copy()
        self.deepSed = np.where(tmpQ<0.)[0]
        tmpQ[tmpQ<0] = 0.
        Qs = self.vland*tmpQ
        Qs[self.seaID] = np.multiply(tmpQ[self.seaID],Qw[self.seaID]) #self.vsea*tmpQ[self.seaID]
        Qs[self.seaID] = np.divide(Qs[self.seaID], self.FVmesh_area[self.seaID],
                            out=np.zeros_like(Qs[self.seaID]), where=self.FVmesh_area[self.seaID]!=0)
        depo = np.divide(Qs, Qw, out=np.zeros_like(Qw), where=Qw!=0)
        depo = np.multiply(depo,self.dt)
        depo = np.divide(depo,1.0-self.phi)

        # Get combined erosion thickness during time step
        Eb = self.EbLocal.getArray().copy()
        Es = self.EsLocal.getArray().copy()

        # Limit deposition if required
        hArrayLocal = self.hLocal.getArray().copy()
        hArrayLocal -= np.multiply(Eb+Es,self.dt)
        dhmax = minHeight(self.inIDs, hArrayLocal)
        dhmax = np.multiply(0.9,dhmax)
        dhmax[dhmax<0] = 0.
        dhmax[self.seaID] = np.maximum(0.,self.seaDep[self.seaID]) #0.9*(self.sealevel-hArrayLocal[self.seaID])-0.1)
        Ds = np.minimum(dhmax,depo)
        # print 'depo limit',Ds.max(),Ds.min()
        Ds[Ds<0] = 0.
        Ds = np.multiply(Ds,1.0-self.phi)
        Ds[self.idGBounds] = 0.
        del Eb, Es, hArrayLocal, dhmax, Qs, depo
        Ds = np.divide(Ds,self.dt)
        self.vSedLocal.setArray(Ds)
        self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)

        return Ds

    def _getDepositionVolume(self, Ds):
        """
        Compute sediment deposition volume.
        """

        # Remove sediment flux on domain boundaries
        tmpQ = self.vSedLocal.getArray().copy()
        tmpQ[self.idGBounds] = 0.
        self.vSedLocal.setArray(tmpQ)
        self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)

        # From deposition rate to volume
        Ds = np.multiply(Ds,self.FVmesh_area)
        tmpQ = np.divide(Ds*self.dt,1.0-self.phi)

        # Set marine deposition volume
        depo = np.zeros(tmpQ.shape)
        depo[self.seaID] = tmpQ[self.seaID]
        self.vLoc.setArray(depo)
        tmp = np.divide(depo, self.FVmesh_area, out=np.zeros_like(tmpQ), where=self.FVmesh_area!=0)
        del tmp

        # Update land deposits thicknesses
        tmpQ[self.seaID] = 0.
        depo = np.divide(tmpQ, self.FVmesh_area, out=np.zeros_like(tmpQ), where=self.FVmesh_area!=0)
        self.vecL.setArray(depo)
        self.dm.localToGlobal(self.vecL, self.vecG, 1)

        del tmpQ, depo, Ds

        return

    def StreamPowerLaw(self):
        """
        Perform stream power law.
        """

        t0 = clock()

        # Get erosion values for considered time step
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal, 1)
        self.hOldArray = self.hOldLocal.getArray()
        self.Es.set(0.)
        self.Eb.set(0.)
        Hsoil = self.Hsoil.getArray().copy()
        if self.rainFlag:
            self._getErosionRate(Hsoil)

        # Get sediment flux rate
        if self.rainFlag and self.frac_fine < 1.:
            Qw = self._getSedFlux(True, None, None)
        else:
            self.vSed.set(0.)
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)
        self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)

        # Get deposition volume
        if self.rainFlag and self.frac_fine < 1.:
            Ds = self._getMaxDepositionRate(Qw)
            self._getSedFlux(False, Ds, Qw)
            self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)
            self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)
            self._getDepositionVolume(Ds)
            del Qw
        else:
            self.vLoc.set(0.)
            self.vecG.set(0.)

        # Update bedrock thicknesses due to erosion
        Eb = self.Eb.getArray()
        Eb *= self.dt
        Eb[Eb<0] = 0.

        # Update sediment thicknesses due to erosion
        if self.Ksed > 0. and self.frac_fine < 1.:
            Es = self.Es.getArray()
            Es *= self.dt
            ids = Es>Hsoil
            Es[ids] = Hsoil[ids]
            Es[Es<0] = 0.
            del ids
        else:
            Es = np.zeros(Eb.shape)
        del Hsoil

        # Update parameters
        self.stepED.setArray(-Es-Eb)
        self.stepED.axpy(1.,self.vecG)
        self.cumED.axpy(1.,self.stepED)
        # self.dm.globalToLocal(self.stepED, self.vecL, 1)

        self.hGlobal.axpy(1.,self.stepED)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        self.stepED.setArray(-Es)
        self.stepED.axpy(1.,self.vecG)

        if self.Ksed > 0. and self.frac_fine < 1.:
            self.Hsoil.axpy(1.,self.stepED)
            self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
            self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)

        del Es, Eb
        if MPIrank == 0 and self.verbose:
            print('Compute Stream Power Law (%0.02f seconds)'% (clock() - t0))

        return

    def SmoothingDeposit(self):
        """
        Perform smoothing of newly deposited sediments.
        """

        t0 = clock()
        self.vGlob.copy(result=self.vecG)
        self._solve_KSP(True, self.Smooth, self.vecG, self.vGlob)

        self.dm.globalToLocal(self.vGlob, self.vLoc, 1)
        tmpQ = self.vLoc.getArray().copy()
        tmpQ[self.idGBounds] = 0.
        self.vLoc.setArray(tmpQ)
        self.dm.localToGlobal(self.vLoc, self.vGlob, 1)
        del tmpQ

        if MPIrank == 0 and self.verbose:
            print('Smoothing of Deposited Sediments (%0.02f seconds)'% (clock() - t0))

        return

    def SedimentDiffusion(self):
        """
        Initialise sediment diffusion from pit and marine deposition.
        """

        t0 = clock()

        # Deposition in depressions
        self.depositDepression()
        if MPIrank == 0 and self.verbose:
            print('Fill Pit Depression (%0.02f seconds)'% (clock() - t0))

        t0 = clock()
        if self.streamCd > 0 or self.oceanCd > 0:
            self._diffuseSediment()
            if MPIrank == 0 and self.verbose:
                print('Compute Sediment Diffusion (%0.02f seconds)'% (clock() - t0))
        else:
            # self.vLoc.set(0.)
            # self.diffDep.set(0.)
            self.dm.localToGlobal(self.vLoc, self.vGlob, 1)
            self.vGlob.axpy(1.,self.diffDep)
            self.vGlob.pointwiseDivide(self.vGlob,self.areaGlobal)
            self.hGlobal.axpy(1.,self.vGlob)
            self.cumED.axpy(1.,self.vGlob)

            self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
            self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
            self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
            self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

            if self.Ksed > 0. and self.frac_fine < 1.:
                self.Hsoil.axpy(1.,self.vGlob)
                Hsoil = self.Hsoil.getArray().copy()
                Hsoil[Hsoil<0.] = 0.
                self.Hsoil.setArray(Hsoil)
                self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)

        return

    def _diffuseSediment(self):
        """
        Marine and aerial diffusion for deposited sediments.
        """

        t1 = clock()
        # Constant local & global vectors/arrays
        self.hGlobal.copy(result=self.hG0)
        self.hLocal.copy(result=self.hL0)
        hG0 = self.hG0.getArray().copy()

        # Add sediment volume from pits to diffuse
        self.dm.localToGlobal(self.vLoc, self.vGlob, 1)
        self.vGlob.axpy(1.,self.diffDep)
        self.dm.globalToLocal(self.vGlob, self.vLoc, 1)

        # From volume to sediment thickness
        self.vGlob.pointwiseDivide(self.vGlob,self.areaGlobal)

        # Get maximum diffusion iteration number
        iters = int((self.vGlob.max()[1]+1.)*2.0)
        iters = max(10,iters)
        if iters < self.maxIters:
            iters = self.maxIters
        itflux =  int(iters*0.5)
        # if iters > self.maxIters:
        #     # Smoothing of newly deposited sediments
        #     while iters > self.maxIters:
        #         self.SmoothingDeposit()
        #         iters = int((self.vGlob.max()[1]+1.)*2.0)
        #         iters = max(10,iters)
        #         if iters < self.maxIters:
        #             iters = self.maxIters
        #         itflux =  int(iters*0.5)
        #     iters *= 2

        # Prepare diffusion arrays
        self.vGlob.scale(2.0/float(iters))
        sFlux = self.vGlob.getArray().copy()
        self.hGlobal.setArray(hG0+sFlux)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        # Nothing to diffuse...
        if self.vGlob.sum() <= 0. :
            del hG0, sFlux
            return
        del hG0

        # Solve temporal diffusion equation
        self.dm.globalToLocal(self.vGlob, self.vLoc, 1)
        sFlux = self.vLoc.getArray().copy()
        ierr = diffusionDT(self.dm.fortran, self.hLocal.fortran, self.hL0.fortran,
                           self.gbounds, iters, itflux, self.inIDs, sFlux, self.streamKd,
                           self.oceanKd, self.sealevel)
        if ierr: raise PETSc.Error(ierr)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        # Cleaning
        del sFlux

        # Update erosion/deposition local/global vectors
        self.stepED.waxpy(-1.0,self.hG0,self.hGlobal)
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
            print('Deposited sediment diffusion (%0.02f seconds)'% (clock() - t1))

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
