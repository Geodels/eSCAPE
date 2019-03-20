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

import sys,petsc4py
petsc4py.init(sys.argv)
from petsc4py import PETSc
from time import clock
import warnings;warnings.simplefilter('ignore')

from eSCAPE._fortran import setHillslopeCoeff
from eSCAPE._fortran import setDiffusionCoeff
from eSCAPE._fortran import explicitDiff
from eSCAPE._fortran import MFDreceivers
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

        val = np.zeros(1,dtype=int)
        X = np.ma.masked_equal(self.FVmesh_edgeLgt,0)
        val[0] = X.min()
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, val, op=MPI.MIN)
        minlgth = val[0]
        dt = np.divide(np.square(val[0]),2.*self.sedimentK)
        self.diff_step = int(self.dt/dt) + 1
        self.diff_dt = self.dt/float(self.diff_step)
        del X, val

        # KSP solver parameters
        self.rtol = 1.0e-8

        data = np.zeros(1,dtype=int)
        data[0] = self.FVmesh_ngbNbs.max()
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, data, op=MPI.MAX)
        self.maxNgbhs = data[0]

        self.excess = False

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
        self.sedLoadLocal = self.hLocal.duplicate()
        self.FAG = self.hGlobal.duplicate()
        self.FAL = self.hLocal.duplicate()
        self.vecG = self.hGlobal.duplicate()
        self.vecL = self.hLocal.duplicate()
        self.vGlob = self.hGlobal.duplicate()
        self.vLoc = self.hLocal.duplicate()
        self.seaG = self.hGlobal.duplicate()
        self.seaL = self.hLocal.duplicate()
        self.tmpG = self.hGlobal.duplicate()
        self.tmpL = self.hLocal.duplicate()

        tmp = self.Eb.getArray()
        self.shape = tmp.shape
        del tmp

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

    def _buildFlowDirection(self, filled=False):
        """
        Build multiple flow direction based on neighbouring slopes.
        """

        t0 = clock()

        # Account for marine regions
        h1 = self.hLocal.getArray().copy()
        self.seaID = np.where(h1<=self.sealevel)[0]
        if not filled:
            h2 = h1.copy()
        else:
            h2 = self.fillLocal.getArray().copy()

        # Define multiple flow directions
        self.rcvID, self.slpRcv, self.distRcv, self.wghtVal = MFDreceivers(self.flowDir, self.inIDs, h2)

        # Set depression nodes
        if not filled:
            # Account for pit regions
            self.pitID = np.where(self.slpRcv[:,0]<=0.)[0]
            self.rcvID[self.pitID,:] = np.tile(self.pitID,(self.flowDir,1)).T
            self.rcvID[self.seaID,:] = np.tile(self.seaID,(self.flowDir,1)).T
            self.distRcv[self.pitID,:] = 0.
            self.distRcv[self.seaID,:] = 0.
            self.wghtVal[self.pitID,:] = 0.
            self.wghtVal[self.seaID,:] = 0.
            self.nbPit = 0
            self.depID = np.where(np.isin(self.idLocal,self.pitID))[0]
            self.depID = np.setdiff1d(self.depID,self.idGBounds)
            if self.depID is not None:
                self.nbPit = len(self.depID)
        else:
            self.pitID = np.where(h1<h2)[0]
        del h1, h2

    	if MPIrank == 0 and self.verbose:
            print('Flow Direction declaration (%0.02f seconds)'% (clock() - t0))

        return

    def FlowAccumulation(self, filled=False):
        """
        Compute multiple flow accumulation.
        """

        self._buildFlowDirection(filled)

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
            if self.tNow == self.tStart and not filled:
                self._solve_KSP(False, WAtrans, self.bG, self.drainArea)
            else:
                if filled:
                    self._solve_KSP(True, WAtrans, self.bG, self.FAG)
                else:
                    self._solve_KSP(True, WAtrans, self.bG, self.drainArea)
            WAMat.destroy()
            WAtrans.destroy()
            WeightMat.destroy()
        else:
            self.drainArea.set(0.)

        if filled:
            self.dm.globalToLocal(self.FAG, self.FAL, 1)
        else:
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
        Kbr[self.pitID] = 0.
        Ksed = np.power(Kcoeff,self.msed)*self.Ksed*self.dt
        Ksed[self.seaID] = 0.
        Ksed[self.pitID] = 0.

        # Initialise identity matrices...
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
            E[self.pitID] = 0.
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
            E[self.pitID] = 0.
            self.EsLocal.setArray(E)
            del ids, E, Ecrit
        else:
            self.EsLocal.set(0.)
        self.dm.localToGlobal(self.EsLocal, self.Es, 1)

        del Kcoeff, Kbr, Ksed

        return

    def cptErosion(self):
        """
        Compute erosion using stream power law.
        """

        t0 = clock()

        # Constant local & global vectors/arrays
        self.Es.set(0.)
        self.Eb.set(0.)
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal, 1)
        self.hOldArray = self.hOldLocal.getArray().copy()

        # Get erosion rate values for considered time step
        Hsoil = self.Hsoil.getArray().copy()
        if self.rainFlag:
            self._getErosionRate(Hsoil)

        # Update bedrock thicknesses due to erosion
        if  self.Kbr > 0.:
            Eb = self.Eb.getArray().copy()
            Eb *= self.dt
            Eb[Eb<0] = 0.
        else:
            Eb = np.zeros(self.shape)

        # Update sediment thicknesses due to erosion
        if self.Ksed > 0. and self.frac_fine < 1.:
            Es = self.Es.getArray().copy()
            Es *= self.dt
            ids = Es>Hsoil
            Es[ids] = Hsoil[ids]
            Es[Es<0] = 0.
            del ids
        else:
            Es = np.zeros(self.shape)
        del Hsoil

        # Update Elevation due to Erosion
        self.stepED.setArray(-Es-Eb)

        self.cumED.axpy(1.,self.stepED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        self.hGlobal.axpy(1.,self.stepED)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

        if self.Ksed > 0. and self.frac_fine < 1.:
            self.stepED.setArray(-Es)
            self.Hsoil.axpy(1.,self.stepED)
            self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
            self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)

        del Es, Eb

        if MPIrank == 0 and self.verbose:
            print('Get Erosion Thicknesses (%0.02f seconds)'% (clock() - t0))

        return

    def _getSedFlux(self):
        """
        Compute sediment flux.
        """

        Qw = self.drainAreaLocal.getArray().copy()

        if self.vland > 0.:
            Qs = self.vland*self.FVmesh_area
            Qs[self.seaID] = 0.
            Qs[self.pitID] = 0.

            # Define deposition thickness...
            depo = np.divide(Qs, Qw, out=np.zeros_like(Qw), where=Qw!=0)
            depo = np.multiply(depo,self.dt)
            depo = np.divide(depo,1.0-self.phi)

            # Limit deposition if required
            dhmax = np.zeros(depo.shape)
            hArrayLocal = self.hLocal.getArray().copy()
            dhmax = minHeight(self.inIDs, hArrayLocal)
            dhmax = np.multiply(0.9,dhmax)
            dhmax[dhmax<0] = 0.
            dhmax[self.idGBounds] = 0.
            del hArrayLocal

            # Limit the maximum aerial deposition
            Ds = np.minimum(dhmax,depo)
            Ds[Ds<0] = 0.
            Ds[self.seaID] = 0.
            Ds[self.pitID] = 0.
            Ds[self.idGBounds] = 0.
            del dhmax, depo

            # Define the deposition rate
            Qs = np.multiply(Ds,1.0-self.phi)
            Qs = np.divide(Ds,self.dt)
            Qs = np.multiply(Ds,Qw)

            # Update Elevation due to Deposition
            self.stepEDL.setArray(Ds)
            self.dm.localToGlobal(self.stepEDL, self.stepED, 1)
            self.cumED.axpy(1.,self.stepED)
            self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
            self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

            self.hGlobal.axpy(1.,self.stepED)
            self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
            self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)

            if self.Ksed > 0. and self.frac_fine < 1.:
                self.Hsoil.axpy(1.,self.stepED)
                self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
                self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)
            del Ds
        else:
            Qs = np.zeros(Qw.shape)

        Kcoeff = np.divide(Qs, Qw, out=np.zeros_like(Qw), where=Qw!=0)
        SLMat = self._matrix_build_diag(Kcoeff)
        SLMat += self.WeightMat
        del Kcoeff, Qs, Qw

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

        return

    def cptSedFlux(self):
        """
        Compute sediment flux.
        """

        t0 = clock()

        # Get erosion rate values for considered time step
        if self.rainFlag:
            # Get sediment flux rate
            self._getSedFlux()
            self.vSedLocal.copy(result=self.sedLoadLocal)
        else:
            self.vSed.set(0.)
            self.vSedLocal.set(0.)
            self.sedLoadLocal.set(0.)

        if MPIrank == 0 and self.verbose:
            print('Update Sediment Load (%0.02f seconds)'% (clock() - t0))

        return

    def depositDepressions(self):
        """
        Compute deposition in depressions.
        """

        t0 = clock()

        self.excess = False

        # In case there is no depression
        if self.pHeight is None:
            return

        # In case there is only suspended sediments
        if self.frac_fine == 1.:
            return

        # Sediment flux and volume
        Qs = self.vSedLocal.getArray().copy()
        depo = np.zeros(Qs.shape)
        depV = np.zeros(Qs.shape)
        zsea = np.zeros(Qs.shape)
        zsea.fill(self.sealevel)

        depV[self.depID] = Qs[self.depID]*self.dt
        depV[self.seaID] = Qs[self.seaID]*self.dt

        # Watershed ID
        tmp = self.shedIDLocal.getArray().copy()
        wshed = -np.ones(tmp.shape)
        wshed[self.idLocal] = tmp[self.idLocal]

        hBase = self.hLocal.getArray().copy()
        hTop = self.fillLocal.getArray().copy()
        hDiff = hTop-hBase

        # Find inland depressions
        inlandIDs = self.pHeight>self.sealevel
        excess = np.zeros(1)
        for k in range(len(self.pVol)):
            if inlandIDs[k]:
                # pts = wshed==k
                pts = np.where(wshed==k)[0]
                zsea[pts] = -1.e6
                locVol = np.zeros(1)
                locVol[0] = np.sum(depV[pts])
                MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, locVol, op=MPI.SUM)
                # Case 1: incoming sediment volume lower than pit volume
                if locVol[0] < self.pVol[k]:
                    frac = locVol[0]/self.pVol[k]
                    depo[pts] += frac*hDiff[pts]
                    Qs[pts] = 0.
                    depV[pts] = 0.
                    self.pVol[k] -= locVol[0]
                # Case 2: incoming sediment volume greater than pit volume
                elif locVol[0] > self.pVol[k] and wshed[self.pitNode[k]]>=0:
                    depo[pts] = hDiff[pts]
                    Qs[pts] = 0.
                    excess[0] = 1
                    self.pVol[k] = 0.
                    depV[pts] = 0.
                    if MPIrank==self.pitProc[k]:
                        Qs[self.pitNode[k]] = (locVol[0] - self.pVol[k])/self.dt
                        depV[self.pitNode[k]] = (locVol[0] - self.pVol[k])
                # Case 3: incoming sediment volume equal pit volume
                else:
                    depo[pts] = hDiff[pts]
                    Qs[pts] = 0.
                    depV[pts] = 0.
                    self.pVol[k] = 0.

        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, excess, op=MPI.MAX)
        if excess[0] > 0:
            self.excess = True

        self.tmpL.setArray(zsea)
        self.dm.localToGlobal(self.tmpL, self.tmpG, 1)
        self.dm.globalToLocal(self.tmpG, self.tmpL, 1)

        # Define deposition from local to global
        self.vecL.setArray(depo)
        self.dm.localToGlobal(self.vecL, self.vecG, 1)
        self.dm.globalToLocal(self.vecG, self.vecL, 1)

        # Cumulative erosion deposition
        self.cumED.axpy(1.,self.vecG)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        # Update elevation and soil thickness if any
        self.hGlobal.axpy(1.,self.vecG)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        if self.Ksed > 0. and self.frac_fine < 1.:
            self.Hsoil.axpy(1.,self.vecG)
            self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
            self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)

        # Update local vector
        self.vSedLocal.setArray(np.divide(depV,self.dt))
        self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)

        if MPIrank == 0 and self.verbose:
            print('Compute Deposition in Land Depressions (%0.02f seconds)'% (clock() - t0))

        del hBase, hTop, hDiff, wshed, depo, Qs, tmp, zsea

        return

    def _sedFlux(self):
        """
        Move excesss sediment volume across filled topography.
        """

        self.vSed.copy(result=self.stepED)
        self._solve_KSP(False, self.WeightMat, self.stepED, self.vSed)

        # Update local vector
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)
        self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)

        return

    def downSediment(self, Qw=None, type=0):
        """
        Transport excesss sediment volume to downstream regions.
        """

        t0 = clock()
        while self.excess:
            self.FlowAccumulation(filled=False)
            self._sedFlux()
            self.depositDepressions()

        if MPIrank == 0 and self.verbose:
            print('Distribute sediment downstream (%0.02f seconds)'% (clock() - t0))

        return

    def marineDeposition(self):
        """
        Deposit sediment in marine environment.
        """

        t0 = clock()

        # In case there is only suspended sediments
        if self.frac_fine == 1.:
            return

        # Get underwater directions
        self.FlowAccumulation(filled=True)

        # Define maximum marine volumic deposition rate
        hArray = self.hGlobal.getArray().copy()
        zsea = self.tmpG.getArray().copy()

        dRate = 0.9*(zsea-hArray)/self.dt
        ids = dRate<0.
        dRate[ids] = 0.

        # Get sediment flux
        self.stepED.setArray(-dRate)
        self.stepED.pointwiseMult(self.stepED,self.areaGlobal)
        self.stepED.axpy(1.,self.vSed)
        self._solve_KSP(False, self.WeightMat, self.stepED, self.vSed)

        # Update local vector
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)
        self.dm.localToGlobal(self.vSedLocal, self.vSed, 1)

        # Get marine deposition
        self.vSed.pointwiseDivide(self.vSed,self.areaGlobal)
        Qs = self.vSed.getArray().copy()
        ids = Qs>0.
        Ds = np.zeros(hArray.shape)
        Ds[ids] = dRate[ids]*self.dt
        Qs[ids] = 0.
        dRate[ids] = 0.

        dH = np.multiply(dRate+Qs,self.dt)
        ids = dH>0.
        Ds[ids] = dH[ids]

        # Add deposition from global to local
        self.vecG.setArray(Ds)
        self.dm.globalToLocal(self.vecG, self.vecL, 1)
        self.dm.localToGlobal(self.vecL, self.vecG, 1)

        # Cumulative erosion deposition
        self.cumED.axpy(1.,self.vecG)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        # Update elevation and soil thickness if any
        self.hGlobal.axpy(1.,self.vecG)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        if self.Ksed > 0. and self.frac_fine < 1.:
            self.Hsoil.axpy(1.,self.vecG)
            self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
            self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)

        del dRate, hArray, ids, Qs, zsea, Ds, dH

        if MPIrank == 0 and self.verbose:
            print('Distribute marine sediment (%0.02f seconds)'% (clock() - t0))

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
        if self.seaG.sum() <= 0.:
            return

        it = 0
        flag = False
        limit = 1.e-12
        dt = self.dt/float(self.maxIters)
        step = int(self.dt/dt) + 1

        # Elevation without freshly deposited sediment thickness
        self.hLocal.copy(result=self.hOldLocal)
        h0 = self.hOldLocal.getArray().copy()

        # Local elevation adding newly deposited thickness
        self.hLocal.axpy(1.,self.seaL)
        hL0 = self.hLocal.getArray().copy()

        while(it<step):
            elev = hL0
            self.hLocal.setArray(elev)
            self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
            self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
            elev = self.hLocal.getArray().copy()
            dh = elev-h0
            dh[dh<0] = 0.
            sedDiff = self._diffuse_TopSed(limit, elev, h0-self.sealevel, dh, dt)
            self._solve_KSP(flag, sedDiff, self.hGlobal, self.tmpG)
            self.dm.globalToLocal(self.tmpG, self.tmpL, 1)
            self.dm.localToGlobal(self.tmpL, self.tmpG, 1)
            self.tmpL.copy(result=self.hLocal)
            self.tmpG.copy(result=self.hGlobal)
            sedDiff.destroy()
            hL0 = self.hLocal.getArray().copy()
            it += 1
            flag = True

        # Update elevation change vector
        self.tmpL.setArray(dh)
        self.dm.localToGlobal(self.tmpL, self.stepED, 1)
        self.dm.globalToLocal(self.stepED, self.tmpL, 1)
        del dh, elev, hL0, h0

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

    def _distSedExplicit(self):
        """
        Distribute newly deposited sediments explicitly
        """

        t1 = clock()

        # Nothing to diffuse...
        if self.seaG.sum() <= 0.:
            return

        it = 0
        limit = 1.e-12

        # Elevation without freshly deposited sediment thickness
        self.hLocal.copy(result=self.hOldLocal)
        h0 = self.hOldLocal.getArray().copy()

        # Local elevation adding newly deposited thickness
        self.hLocal.axpy(1.,self.seaL)
        hL0 = self.hLocal.getArray().copy()

        while(it<self.diff_step):
            elev = hL0
            self.hLocal.setArray(elev)
            self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
            self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
            elev = self.hLocal.getArray().copy()
            dh = elev-h0
            dh[dh<0] = 0.
            # Diffusion matrix construction
            newz = explicitDiff(self.sedimentK*self.diff_dt, limit, elev,
                                      h0-self.sealevel, dh)
            self.tmpL.setArray(newz)
            self.dm.localToGlobal(self.tmpL, self.tmpG, 1)
            self.dm.globalToLocal(self.tmpG, self.tmpL, 1)
            self.tmpL.copy(result=self.hLocal)
            self.tmpG.copy(result=self.hGlobal)
            hL0 = self.hLocal.getArray().copy()
            it += 1

        # Update elevation change vector
        self.tmpL.setArray(dh)
        self.dm.localToGlobal(self.tmpL, self.stepED, 1)
        self.dm.globalToLocal(self.stepED, self.tmpL, 1)
        del dh, elev, hL0, h0

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
            self.seaL.waxpy(-1.,self.hOldLocal,self.hLocal)
            dh = self.seaL.getArray().copy()
            depSea = np.zeros(dh.shape)
            depSea[self.seaID] = dh[self.seaID]
            depSea[depSea<0] = 0.
            self.seaL.setArray(depSea)
            self.dm.localToGlobal(self.seaL, self.seaG, 1)
            self.dm.globalToLocal(self.seaG, self.seaL, 1)
            del depSea
            self.hLocal.axpy(-1.,self.seaL)
            self.hGlobal.axpy(-1.,self.seaG)
            self.cumED.axpy(-1.,self.seaG)

            if self.Ksed > 0. and self.frac_fine < 1.:
                 self.Hsoil.axpy(-1.,self.seaG)
                 Hsoil = self.Hsoil.getArray().copy()
                 Hsoil[Hsoil<0.] = 0.
                 self.Hsoil.setArray(Hsoil)
                 self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
                 del Hsoil

            self._distSedExplicit()
            # self._distributeSediment()
            if MPIrank == 0 and self.verbose:
                print('Compute Sediment Diffusion (%0.02f seconds)'% (clock() - t0))

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
