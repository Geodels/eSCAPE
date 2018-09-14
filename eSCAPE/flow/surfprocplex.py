"""
Copyright 2017-2018 Tristan Salles

This file is part of eSCAPE.

eSCAPE is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

eSCAPE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with eSCAPE.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
from mpi4py import MPI
from scipy import sparse
import sys,petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from time import clock
import warnings;warnings.simplefilter('ignore')

from eSCAPE._fortran import setHillslopeCoeff
from eSCAPE._fortran import initDiffCoeff
from eSCAPE._fortran import MFDreceivers
from eSCAPE._fortran import distributeHeight
from eSCAPE._fortran import distributeVolume
from eSCAPE._fortran import getMaxEro
from eSCAPE._fortran import getDiffElev

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

        # Petsc vectors
        self.hG0 = self.hGlobal.duplicate()
        self.hL0 = self.hLocal.duplicate()
        self.vecG = self.hGlobal.duplicate()
        self.vecL = self.hLocal.duplicate()
        self.vGlob = self.hGlobal.duplicate()
        self.vLoc = self.hLocal.duplicate()

        return

    def _matrix_build(self, nnz=(1,1)):

        matrix = PETSc.Mat().create(comm=MPIcomm)
        matrix.setType('aij')
        matrix.setSizes(self.sizes)
        matrix.setLGMap(self.lgmap_row, self.lgmap_col)
        matrix.setFromOptions()
        matrix.setPreallocationNNZ(nnz)

        return matrix

    def _matrix_build_diag(self, V, nnz=(1,1)):

        matrix = self._matrix_build()

        # Define diagonal matrix
        I = np.arange(0, self.npoints+1, dtype=PETSc.IntType)
        J = np.arange(0, self.npoints, dtype=PETSc.IntType)
        matrix.assemblyBegin()
        matrix.setValuesLocalCSR(I, J, V, PETSc.InsertMode.INSERT_VALUES)
        matrix.assemblyEnd()

        return matrix

    def _solve_KSP(self, guess, matrix, vector1, vector2):

        ksp = PETSc.KSP().create(PETSc.COMM_WORLD)
        if guess:
            ksp.setInitialGuessNonzero(guess)
        ksp.setOperators(matrix,matrix)
        ksp.setType('richardson')
        pc = ksp.getPC()
        pc.setType('bjacobi')
        ksp.setTolerances(rtol=self.rtol)
        ksp.solve(vector1, vector2)
        ksp.destroy()

        return vector2

    def _buidFlowDirection(self):

        t0 = clock()
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        hArrayLocal = self.hLocal.getArray()
        self.rcvID, self.slpRcv, self.distRcv, self.wghtVal = MFDreceivers(self.flowDir, self.inIDs, hArrayLocal)
        # Account for pit regions
        self.pitID = np.where(self.slpRcv[:,0]<=0.)[0]
        self.rcvID[self.pitID,:] = np.tile(self.pitID,  (self.flowDir,1)).T
        self.distRcv[self.pitID,:] = 0.
        self.wghtVal [self.pitID,:] = 0.
        # Account for marine regions
        self.seaID = np.where(hArrayLocal<self.sealevel)[0]
        self.rcvID[self.seaID,:] = np.tile(self.seaID,(self.flowDir,1)).T
        self.distRcv[self.seaID,:] = 0.
        self.wghtVal [self.seaID,:] = 0.
        del hArrayLocal

    	if MPIrank == 0 and self.verbose:
            print('Flow Direction declaration (%0.02f seconds)'% (clock() - t0))

        return

    def FlowAccumulation(self):

        self._buidFlowDirection()

        t0 = clock()
        # Build drainage area matrix
        if self.rainFlag:
            WAMat = self.iMat.copy()
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
                WAMat += tmpMat
                tmpMat.destroy()

            # Solve flow accumulation
            WAtrans = WAMat.transpose()
            self.WAtrans = WAtrans.copy()
            self.dm.localToGlobal(self.bL, self.bG, 1)
            self.dm.globalToLocal(self.bG, self.bL, 1)
            if self.tNow == self.tStart:
                self._solve_KSP(False,WAtrans, self.bG, self.drainArea)
            else:
                self._solve_KSP(True,WAtrans, self.bG, self.drainArea)
            WAMat.destroy()
            WAtrans.destroy()
        else:
            self.drainArea.set(0.)
        self.dm.globalToLocal(self.drainArea, self.drainAreaLocal, 1)

        if MPIrank == 0 and self.verbose:
            print('Compute Flow Accumulation (%0.02f seconds)'% (clock() - t0))

        return

    def StreamPowerLaw(self):

        t0 = clock()
        # Get erosion values for considered time step
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal, 1)
        self.hOldArray = self.hOldLocal.getArray()

        if self.rainFlag:
            # Define linear solver coefficients
            Kcoeff = self.drainAreaLocal.getArray()
            Kcoeff = -np.power(Kcoeff,self.m)*self.dt*self.Ke
            WHMat = self.iMat.copy()
            for k in range(0, self.flowDir):
                # Erosion processes computation
                data = np.divide(Kcoeff, self.distRcv[:,k], out=np.zeros_like(Kcoeff),
                                            where=self.distRcv[:,k]!=0)
                tmpMat = self._matrix_build()
                data = np.multiply(data,self.wghtVal[:,k])
                indptr = np.arange(0, self.npoints+1, dtype=PETSc.IntType)
                nodes = indptr[:-1]
                data[self.rcvID[:,k].astype(PETSc.IntType)==nodes] = 0.0
                tmpMat.assemblyBegin()
                tmpMat.setValuesLocalCSR(indptr, self.rcvID[:,k].astype(PETSc.IntType), data,
                                                                PETSc.InsertMode.INSERT_VALUES)
                tmpMat.assemblyEnd()
                WHMat += tmpMat
                Mdiag = self._matrix_build_diag(data)
                WHMat -= Mdiag
                tmpMat.destroy()
                Mdiag.destroy()
            # Solve flow accumulation
            self._solve_KSP(True, WHMat, self.hOld, self.hGlobal)
            WHMat.destroy()

        self.stepED.waxpy(-1.0,self.hOld,self.hGlobal)
        self.cumED.axpy(1.,self.stepED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        # Update elevation locally due to aerial erosion
        hL0 = self.hLocal.getArray().copy()
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        # Get sediment load value for given time step
        if self.rainFlag:
            # Define erosion deposition volume
            self.stepED.scale(-1.)
            self.stepED.pointwiseMult(self.stepED,self.areaGlobal)
            if self.tNow == self.tStart:
                self._solve_KSP(False,self.WAtrans, self.stepED, self.vSed)
            else :
                self._solve_KSP(True,self.WAtrans, self.stepED, self.vSed)
            # Solution for sediment load
            self.WAtrans.destroy()
        else:
            self.vSed.set(0.)
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)
        self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)

        # Update elevation locally due to erosion
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

        if MPIrank == 0 and self.verbose:
            print('Compute Stream Power Law (%0.02f seconds)'% (clock() - t0))

        return

    def SedimentDiffusion(self):

        t0 = clock()

        # Deposition in depressions
        self.depositDepression()
        if MPIrank == 0 and self.verbose:
            print('Fill Pit Depression (%0.02f seconds)'% (clock() - t0))

        t0 = clock()
        if self.streamCd > 0 and self.oceanCd > 0:
            self._diffuseSediment()
            if MPIrank == 0 and self.verbose:
                print('Compute Sediment Diffusion (%0.02f seconds)'% (clock() - t0))

        return

    def _diffusionTimeStep(self, hArr, hG0, hL0, sFlux):

        shape = hG0.shape
        for tStep in range(self.iters):

            # Solve PDE for diffusion system explicitly
            self.hGlobal.copy(result=self.hOld)
            eroCoeffs = getMaxEro( self.sealevel, self.inIDs, hArr, hL0, self.streamKd, self.oceanKd)
            self.vLoc.setArray(eroCoeffs)
            self.dm.localToGlobal(self.vLoc, self.vGlob, 1)
            self.dm.globalToLocal(self.vGlob, self.vLoc, 1)
            eroCoeffs = self.vLoc.getArray()
            dh = getDiffElev( self.sealevel, self.inIDs, hArr, eroCoeffs, self.streamKd, self.oceanKd)
            hArr += dh
            self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

             # Remove portion of boundary deposits
            self.hGlobal.copy(result=self.vGlob)
            self.vGlob.axpy(-1.,self.hOld)
            tmpArr = np.zeros(shape)
            tmpArr[self.boundGIDs] = 0.
            self.vecG.setArray(tmpArr)
            self.vecG.pointwiseMult(self.vecG,self.vGlob)
            self.hGlobal.axpy(-1.,self.vecG)

            # Update elevation vectors
            if tStep < int(self.iters*0.5)-1:
                tmpArr = self.hGlobal.getArray()
                tmpArr += sFlux

            self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

            # Update upper layer fraction
            hArr = self.hLocal.getArray()

        # Cleaning
        del eroCoeffs, tmpArr

        return

    def _diffuseSediment(self):

        # Constant local & global vectors/arrays
        self.hGlobal.copy(result=self.hG0)
        hG0 = self.hG0.getArray().copy()
        hL0 = self.hLocal.getArray().copy()

        # Get sediment volume in the marine environment
        vSed = self.vSedLocal.getArray().copy()
        hL = self.hLocal.getArray().copy()
        vSea = np.zeros(vSed.shape)
        vSea[self.seaID] = vSed[self.seaID]
        vSea[self.idGBounds] = 0.

        # Add inland sediment volume to diffuse
        self.vLoc.setArray(vSea)
        self.dm.localToGlobal(self.vLoc, self.vGlob, 1)
        splitNb = 1
        self.vGlob.axpy(1.,self.diffDep)
        self.vGlob.scale(1./float(splitNb))
        self.dm.globalToLocal(self.vGlob, self.vLoc, 1)
        vSed = self.vLoc.getArray().copy()
        inFlux = vSed.copy()
        it2 = 0
        while it2 < splitNb:
            it = 0
            while it < self.maxIters:
                nhL,nvSed = distributeHeight(self.inIDs, self.sealevel, hL, hL0, vSed)
                # Update local elevation globally
                self.vecL.setArray(nhL)
                self.dm.localToGlobal(self.vecL, self.vecG, 1)
                self.dm.globalToLocal(self.vecG, self.vecL, 1)
                hL = self.vecL.getArray().copy()
                self.vecL.setArray(nvSed)
                self.dm.localToGlobal(self.vecL, self.vecG, 1)
                self.dm.globalToLocal(self.vecG, self.vecL, 1)
                vSed = self.vecL.getArray().copy()
                vSed[self.idGBounds] = 0.
                # Update sediment to distribute
                nvSed = distributeVolume(self.inIDs, self.sealevel, hL, hL0, vSed)
                self.vecL.setArray(nvSed)
                self.dm.localToGlobal(self.vecL, self.vecG, 1)
                self.dm.globalToLocal(self.vecG, self.vecL, 1)
                vSed = self.vecL.getArray().copy()
                vSed[self.idGBounds] = 0.
                self.vGlob.pointwiseDivide(self.vecG,self.areaGlobal)
                hmax = self.vGlob.max()[1]
                if hmax <= 0.1:
                    break
                it += 1
            if it2 < splitNb-1:
                vSed += inFlux
            it2 += 1

        vTot = np.multiply(hL-hL0,self.FVmesh_area)
        vTot += vSed
        self.vLoc.setArray(vTot)

        # Get the height to diffuse in the aerial domain
        nhL = hL-hL0
        nhL[self.seaID] = 0.
        self.vecL.setArray(nhL)
        self.dm.localToGlobal(self.vecL, self.vecG, 1)
        self.iters = int(self.vecG.max()[1])
        del vTot, vSed, hL, nvSed, nhL, inFlux

        # From volume to sediment thickness to distribute at each interval
        self.dm.localToGlobal(self.vLoc, self.vGlob, 1)
        self.vGlob.pointwiseDivide(self.vGlob,self.areaGlobal)

        # Distribute sediment downstream to speed up the diffusion algorithm
        self.iters = max(100,int(self.iters*2.0))
        self.iters = min(self.maxIters,self.iters)
        cumdt = self.iters*self.diffDT
        if self.iters*self.diffDT < self.dt:
            self.iters = int(self.dt/self.diffDT)+1
        self.vGlob.scale(1.0/float(self.iters))
        sFlux = self.vGlob.getArray().copy()

        # Prepare arrays
        self.hGlobal.setArray(hG0+sFlux)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        hArr = self.hLocal.getArray()

        # Nothing to diffuse...
        if self.vGlob.sum() <= 0. :
            del hArr,hL0,hG0
            return

        # Solve temporal diffusion equation
        self._diffusionTimeStep(hArr, hG0, hL0, sFlux)

        # Cleaning
        del hArr,sFlux
        del hG0,hL0

        # Update erosion/deposition local/global vectors
        self.stepED.waxpy(-1.0,self.hG0,self.hGlobal)
        self.cumED.axpy(1.,self.stepED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)
        
        return

    def HillSlope(self):

        t0 = clock()
        if self.Cd > 0.:
            # Get erosion values for considered time step
            self.hGlobal.copy(result=self.hOld)
            self._solve_KSP(True, self.Diff, self.hOld, self.hGlobal)

            # Update cumulative erosion/deposition and elevation
            self.stepED.waxpy(-1.0,self.hOld,self.hGlobal)
            self.cumED.axpy(1.,self.stepED)
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