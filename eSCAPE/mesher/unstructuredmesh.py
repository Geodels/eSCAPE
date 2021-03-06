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

import meshio
import meshplex

from eSCAPE._fortran import defineGTIN
from eSCAPE._fortran import defineTIN
from eSCAPE._fortran import slpBounds
from eSCAPE._fortran import flatBounds

from scipy.spatial import cKDTree as _cKDTree
from scipy.interpolate import CloughTocher2DInterpolator

MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD

try: range = xrange
except: pass

class UnstMesh(object):
    """
    Creating a distributed DMPlex and global vector from it based on triangulated mesh
    """
    def __init__(self, filename, dim=2):

        self.natural2local = None

        # Define mesh attributes on root processor
        t0 = clock()
        t = clock()
        gpoints = np.zeros(1)
        if MPIrank == 0:
            cells = np.asarray(self.mdata.cells['triangle'], dtype=np.int32)
            coords = np.asarray(self.mdata.points, dtype=np.double)
            MPIcomm.bcast(cells.shape, root=0)
            MPIcomm.bcast(coords.shape, root=0)
            elev = self.mdata.point_data[filename[1]]
            gpoints[0] = len(coords)
            if self.verbose:
                print('Extract information (%0.02f seconds)'% (clock() - t))
            t = clock()
            Gmesh = meshplex.mesh_tri.MeshTri(coords, cells)
            if self.verbose:
                print('Reading meshplex meshtri (%0.02f seconds)'% (clock() - t))
            t = clock()
            Gmesh.mark_boundary()
            if self.verbose:
                print('Mark boundary TIN (%0.02f seconds)'% (clock() - t))
            ids = np.arange(0, len(Gmesh.node_coords), dtype=int)
            self.boundGlob = ids[Gmesh._is_boundary_node]
            t = clock()
            Gmesh.create_edges()
            if self.verbose:
                print('Defining edges TIN (%0.02f seconds)'% (clock() - t))
            t = clock()
            self.Gmesh_ngbNbs, self.Gmesh_ngbID = defineGTIN(gpoints[0], Gmesh.cells['nodes'], Gmesh.edges['nodes'])
            if self.verbose:
                print('Defining global TIN (%0.02f seconds)'% (clock() - t))
            del Gmesh, ids
        else:
            cell_shape = list(MPIcomm.bcast(None, root=0))
            coord_shape = list(MPIcomm.bcast(None, root=0))
            cell_shape[0] = 0
            coord_shape[0] = 0
            cells = np.zeros(cell_shape, dtype=np.int32)
            coords = np.zeros(coord_shape, dtype=np.double)
            elev = np.zeros(coord_shape[0], dtype=np.double)
            self.Gmesh_ngbNbs = None
            self.Gmesh_ngbID = None
            gcoords = None

        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gpoints, op=MPI.MAX)
        self.gpoints = int(gpoints[0])
        self.gcoords = MPI.COMM_WORLD.bcast(coords, root=0)

        if MPIrank == 0 and self.verbose:
            print('Reading mesh information (%0.02f seconds)' % (clock() - t0))

        # Create DMPlex
        t = clock()
        self._create_DMPlex(dim, coords, cells, elev)
        if MPIrank == 0 and self.verbose:
            print('Creating Petsc DMPlex (%0.02f seconds)'% (clock() - t))

        # Define local vertex & cells
        t = clock()
        cStart, cEnd = self.dm.getHeightStratum(0)
        # Dealing with triangular cells only
        self.lcells = np.zeros((cEnd-cStart,3), dtype=PETSc.IntType)
        for c in range(cStart, cEnd):
            point_closure = self.dm.getTransitiveClosure(c)[0]
            self.lcells[c,:] = point_closure[-3:]-cEnd
        del point_closure
        if MPIrank == 0 and self.verbose:
            print('Defining local DMPlex (%0.02f seconds)'% (clock() - t))

        # Create mesh structure with meshplex
        t = clock()
        # Define mesh characteristics
        Tmesh = meshplex.mesh_tri.MeshTri(self.lcoords, self.lcells)
        self.FVmesh_area = np.abs(Tmesh.control_volumes)
        self.FVmesh_area[np.isnan(self.FVmesh_area)] = 1.
        self.boundary, self.localboundIDs = self._get_boundary()
        self.gbds = self.boundary.astype(int)

        # Voronoi and simplices declaration
        coords = Tmesh.node_coords
        Tmesh.create_edges()
        cc = Tmesh.cell_circumcenters
        edges_nodes = Tmesh.edges['nodes']
        cells_nodes = Tmesh.cells['nodes']
        cells_edges = Tmesh.cells['edges']

        if MPIrank == 0 and self.verbose:
            print('Voronoi creation (%0.02f seconds)'% (clock() - t))

        # Tesselation on unstructured mesh
        t = clock()
        self.FVmesh_ngbNbs, self.FVmesh_ngbID, self.FVmesh_edgeLgt, \
                self.FVmesh_voroDist = defineTIN(coords, cells_nodes, cells_edges,
                                                 edges_nodes, self.FVmesh_area, cc.T)

        if MPIrank == 0 and self.verbose:
            print('Tesselation (%0.02f seconds)'% (clock() - t))

        self.bL = self.hLocal.duplicate()
        self.bG = self.hGlobal.duplicate()
        self.stepED = self.hGlobal.duplicate()
        self.stepEDL = self.hLocal.duplicate()
        self.vSed = self.hGlobal.duplicate()
        self.vSedLocal = self.hLocal.duplicate()
        self.diffDep = self.hGlobal.duplicate()
        self.diffDepLocal = self.hLocal.duplicate()
        self.hOld = self.hGlobal.duplicate()
        self.hOldLocal = self.hLocal.duplicate()

        self.drainArea = self.hGlobal.duplicate()
        self.drainAreaLocal = self.hLocal.duplicate()
        areaLocal = self.hLocal.duplicate()
        self.areaGlobal = self.hGlobal.duplicate()
        areaLocal.setArray(self.FVmesh_area)

        self.dm.localToGlobal(areaLocal, self.areaGlobal)
        self.dm.globalToLocal(self.areaGlobal, areaLocal)
        self.FVmesh_Garea = self.areaGlobal.getArray().copy()
        areaLocal.destroy()

        self.cumED = self.hGlobal.duplicate()
        self.cumED.set(0.0)
        self.cumEDLocal = self.hLocal.duplicate()
        self.cumEDLocal.set(0.0)

        self.shedID = self.hGlobal.duplicate()
        self.shedIDLocal = self.hLocal.duplicate()

        self.Es = self.hGlobal.duplicate()
        self.Es.set(0.0)
        self.Eb = self.hGlobal.duplicate()
        self.Eb.set(0.0)
        self.EsLocal = self.hLocal.duplicate()
        self.EsLocal.set(0.0)
        self.EbLocal = self.hLocal.duplicate()
        self.EbLocal.set(0.0)

        self._getSoilThickness()

        if MPIrank == 0 and self.verbose:
            print('Finite volume mesh declaration (%0.02f seconds)'% (clock() - t0))

        return

    def _naturalNumbering(self,coords):
        """
        Natural numbering based on DMPlex distribution

        Args:
            coords: mesh coordinates
        """

        self.lcoords = self.dm.getCoordinatesLocal().array.reshape(-1,3)
        self.npoints = self.lcoords.shape[0]

        if MPIsize > 1:
            tree = _cKDTree(coords[:,:2])
            distances, nat2loc = tree.query(self.lcoords[:,:2], 1)
            self.natural2local = nat2loc.copy()
            del tree,distances, nat2loc
        else:
            self.natural2local = np.arange(0,self.npoints,dtype=int)

        return

    def _create_DMPlex(self, dim, coords, cells, elev):
        """
        Create a PETSc DMPlex object from the mesh attributes

        Args:
            dim: mesh dimensions
            coords: mesh coordinates
            cells: cell nodes indices
            elev: nodes elevation
        """

        t0 = clock()
        self.dm = PETSc.DMPlex().createFromCellList(dim, cells, coords, comm=PETSc.COMM_WORLD)
        if MPIrank == 0 and self.verbose:
            print('Create DMPlex (%0.02f seconds)'% (clock() - t0))

        # Create boundary labels
        t0 = clock()
        label = "boundary"
        self._set_DMPlex_boundary_points(label)

        # label coarse DM in case it is ever needed again
        self.dm.createLabel("coarse")
        pStart, pEnd = self.dm.getDepthStratum(0)
        for pt in range(pStart, pEnd):
            self.dm.setLabelValue("coarse", pt, 1)

        # Define one DoF on the nodes
        self.dm.setNumFields(1)
        origSect = self.dm.createSection(1, [1,0,0])
        origSect.setFieldName(0, "points")
        origSect.setUp()
        self.dm.setDefaultSection(origSect)
        origVec = self.dm.createGlobalVector()

        # Distribute to other processors if any
        if MPIsize > 1:
            sf = self.dm.distribute(overlap=1)
            newSect, newVec = self.dm.distributeField(sf, origSect, origVec)
            self.dm.setDefaultSection(newSect)
            newSect.destroy()
            newVec.destroy()
            sf.destroy()
        origVec.destroy()
        origSect.destroy()

        self.hGlobal = self.dm.createGlobalVector()
        self.hLocal = self.dm.createLocalVector()
        self.sizes = self.hGlobal.getSizes(), self.hGlobal.getSizes()

        # Local/Global mapping
        self.lgmap_row = self.dm.getLGMap()
        l2g = self.lgmap_row.indices.copy()
        offproc = l2g < 0
        l2g[offproc] = -(l2g[offproc] + 1)
        self.lgmap_col = PETSc.LGMap().create(l2g, comm=PETSc.COMM_WORLD)
        del l2g

        if MPIrank == 0 and self.verbose:
            print('Distribute DMPlex (%0.02f seconds)'% (clock() - t0))

        # Get natural numbering
        t0 = clock()
        coords = MPI.COMM_WORLD.bcast(coords, root=0)
        elev = MPI.COMM_WORLD.bcast(elev, root=0)
        self._naturalNumbering(coords)

        self.hLocal.setArray(elev[self.natural2local])
        self.dm.localToGlobal(self.hLocal, self.hGlobal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        if MPIrank == 0 and self.verbose:
            print('Distribute field to DMPlex (%0.02f seconds)'% (clock() - t0))

        # Forcing event number
        self.rainNb = -1
        self.tecNb = -1

        return

    def _getSoilThickness(self):
        """
        Specify initial soil thickness
        """

        self.Hsoil = self.hGlobal.duplicate()
        self.HsoilLocal = self.hLocal.duplicate()

        if pd.isnull(self.soildata['sUni'][0]):
            mdata = meshio.read(self.soildata.iloc[0,1])
            soilH = mdata.point_data[self.soildata.iloc[0,2]]
            del mdata
        else:
            soilH = np.full(len(self.elev),self.soildata.iloc[0,0])

        self.HsoilLocal.setArray(soilH[self.natural2local])
        self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)
        self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
        del soilH

        return

    def _get_boundary(self, label="boundary"):
        """
        Find the nodes on the boundary from the DM
        """

        bmask = np.ones(self.npoints, dtype=bool)
        pStart, pEnd = self.dm.getDepthStratum(0)

        labels = []
        for i in range(self.dm.getNumLabels()):
            labels.append(self.dm.getLabelName(i))

        if label not in labels:
            raise ValueError("There is no {} label in the DM".format(label))

        stratSize = self.dm.getStratumSize(label, 1)
        if stratSize > 0:
            labelIS = self.dm.getStratumIS(label, 1)
            pt_range = np.logical_and(labelIS.indices >= pStart, labelIS.indices < pEnd)
            indices = labelIS.indices[pt_range] - pStart
            labelIS.destroy()
        else:
            indices = np.zeros((0,), dtype=np.int)

        bmask[indices] = False

        return bmask, indices

    def updateBoundaries(self):
        """
        Apply boundary forcing for slope and flat conditions.
        """

        t0 = clock()
        if self.tNow > self.tStart :

            hArray = self.hLocal.getArray().copy()
            edArray = self.cumEDLocal.getArray().copy()

            if self.boundCond == 'slope' :
                bElev, bDep = slpBounds(hArray, edArray, self.idGBounds, self.gbds)
            elif self.boundCond == 'flat' :
                bElev, bDep = flatBounds(hArray, edArray, self.idGBounds, self.gbds)
            else:
                bElev = hArray.copy()
                bDep = edArray.copy()

            self.hLocal.setArray(bElev)
            self.cumEDLocal.setArray(bDep)
            del bElev, hArray, bDep, edArray
            self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
            self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)

        if MPIrank == 0 and self.verbose:
            print('Update Boundaries (%0.02f seconds)'% (clock() - t0))

    def applyForces(self):
        """
        Find the different values for climatic and tectonic forces that will be applied to the
        considered time interval
        """

        t0 = clock()
        # Sea level
        self.sealevel = self.seafunction(self.tNow+self.dt)
        # Climate
        self._updateRain()
        # Tectonic
        self._updateTectonic()

        if MPIrank == 0 and self.verbose:
            print('Update External Forces (%0.02f seconds)'% (clock() - t0))

        return

    def _updateRain(self):
        """
        Find the current rain values for the considered time interval
        """

        if self.raindata is None :
            self.rainFlag = False
            return

        nb = self.rainNb
        if nb < len(self.raindata)-1 :
            if self.raindata.iloc[nb+1,0] <= self.tNow+self.dt :
                nb += 1

        if nb > self.rainNb or nb == -1:
            if nb == -1:
                nb = 0

            self.rainNb = nb
            if pd.isnull(self.raindata['rUni'][nb]):
                mdata = meshio.read(self.raindata.iloc[nb,2])
                rainArea = mdata.point_data[self.raindata.iloc[nb,3]]
                del mdata
            else:
                rainArea = np.full(len(self.elev),self.raindata.iloc[nb,1])
            self.rainArea = rainArea[self.natural2local]*self.FVmesh_area

            if rainArea.max() > 0.:
                self.rainFlag = True
            else:
                self.rainFlag = False

        localZ = self.hLocal.getArray()
        seaID = localZ<self.sealevel
        rainArea = self.rainArea.copy()
        rainArea[seaID] = 0.
        rainArea[self.idGBounds] = 0.

        self.bL.setArray(rainArea)
        self.dm.localToGlobal(self.bL, self.bG, 1)
        self.dm.globalToLocal(self.bG, self.bL, 1)
        del rainArea

        return

    def _updateTectonic(self):
        """
        Find the current tectonic rates for the considered time interval
        """

        if self.tecdata is None :
            self.tectonic = None
            return

        nb = self.tecNb
        if nb < len(self.tecdata)-1 :
            if self.tecdata.iloc[nb+1,0] <= self.tNow+self.dt :
                nb += 1

        if nb > self.tecNb or nb == -1:
            if nb == -1:
                nb = 0

            self.tecNb = nb
            if pd.isnull(self.tecdata['tUni'][nb]):

                # Horizontal displacements along X-axis
                flagAdvection = False
                if not pd.isnull(self.tecdata['tMapX'][nb]):
                    mdata = meshio.read(self.tecdata.iloc[nb,2])
                    tectonicX = mdata.point_data[self.tecdata.iloc[nb,5]]
                    flagAdvection = True
                    del mdata
                else:
                    tectonicX = np.zeros(self.gpoints)

                # Horizontal displacements along Y-axis
                if not pd.isnull(self.tecdata['tMapY'][nb]):
                    mdata = meshio.read(self.tecdata.iloc[nb,3])
                    tectonicY = mdata.point_data[self.tecdata.iloc[nb,6]]
                    flagAdvection = True
                    del mdata
                else:
                    tectonicY = np.zeros(self.gpoints)

                # Vertical displacements
                if self.sphere == 0:
                    if not pd.isnull(self.tecdata['tMapZ'][nb]):
                        mdata = meshio.read(self.tecdata.iloc[nb,4])
                        tectonic = mdata.point_data[self.tecdata.iloc[nb,7]]
                        self.tectonic = tectonic[self.natural2local]*self.dt
                        del mdata,tectonic
                    else:
                        tectonic = np.zeros(len(self.elev))
                        self.tectonic = tectonic[self.natural2local]
                        del tectonic
                else:
                    if not pd.isnull(self.tecdata['tMapZ'][nb]):
                        mdata = meshio.read(self.tecdata.iloc[nb,4])
                        tectonicZ = mdata.point_data[self.tecdata.iloc[nb,7]]
                    else:
                        tectonicZ = np.zeros(self.gpoints)

                    tectonic = np.zeros(len(self.elev))
                    self.tectonic = tectonic[self.natural2local]
                    del tectonic

                if flagAdvection:
                    if nb < len(self.tecdata.index)-1:
                        timer = self.tecdata.iloc[nb+1,0]-self.tecdata.iloc[nb,0]
                    else:
                        timer = self.tEnd-self.tecdata.iloc[nb,0]

                    if self.sphere == 0:
                        self._meshAdvector(tectonicX, tectonicY, timer)
                        del tectonicX, tectonicY
                    else:
                        self._meshAdvectorSphere(tectonicX, tectonicY, tectonicZ, timer)
                        del tectonicX, tectonicY, tectonicZ
            else:
                tectonic = np.full(len(self.elev),self.tecdata.iloc[self.tecNb,1])
                self.tectonic = tectonic[self.natural2local]*self.dt
                del tectonic
                # Assign null displacements on the domain boundaries
                if self.boundCond == 'fixed' :
                    self.tectonic[self.idGBounds] = 0.

        # Update elevation
        if self.tNow+self.dt > self.tStart:
            localZ = self.hLocal.getArray()
            localZ += self.tectonic
            self.hLocal.setArray(localZ.reshape(-1))
            self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
            self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        return

    def _meshAdvector(self, tectonicX, tectonicY, timer):
        """
        Advect the mesh horizontally and interpolate mesh information
        """

        # Move horizontally coordinates
        XYZ = np.zeros((self.gpoints,3))
        XYZ[:,0] = self.gcoords[:,0] + tectonicX*timer
        XYZ[:,1] = self.gcoords[:,1] + tectonicY*timer

        # Get mesh variables to interpolate
        # Elevation
        elev = np.zeros(self.gpoints)
        elev.fill(-1.e8)
        elev[self.natural2local] = self.hLocal.getArray().copy()
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, elev, op=MPI.MAX)
        interpolator = CloughTocher2DInterpolator(XYZ[:,:2],elev)
        nelev = interpolator(self.lcoords[:,:2])
        id_NaNs = np.isnan(nelev)
        nelev[id_NaNs] = 0.

        # Erosion/deposition
        erodep = np.zeros(self.gpoints)
        erodep.fill(-1.e8)
        erodep[self.natural2local] = self.cumEDLocal.getArray().copy()
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, erodep, op=MPI.MAX)
        interpolator = CloughTocher2DInterpolator(XYZ[:,:2],erodep)
        nerodep = interpolator(self.lcoords[:,:2])
        nerodep[id_NaNs] = 0.

        # Soil thickness
        if self.Ksed > 0.:
            soil = np.zeros(self.gpoints)
            soil.fill(-1.e8)
            soil[self.natural2local] = self.HsoilLocal.getArray().copy()
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, soil, op=MPI.MAX)
            interpolator = CloughTocher2DInterpolator(XYZ[:,:2],soil)
            nsoil = interpolator(self.lcoords[:,:2])
            nsoil[id_NaNs] = 0.

        # Build kd-tree
        tree = _cKDTree(XYZ)
        distances, indices = tree.query(self.lcoords[id_NaNs,:], k=10)

        # Inverse weighting distance...
        weights = 1.0 / distances**2
        onIDs = np.where(distances[:,0] == 0)[0]

        nelev[id_NaNs] = np.sum(weights*elev[indices],axis=1)/np.sum(weights, axis=1)

        self.hLocal.setArray(nelev)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        del elev, nelev

        self.cumEDLocal.setArray(nerodep)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        del erodep, nerodep

        if self.Ksed > 0.:
            self.HsoilLocal.setArray(nsoil)
            self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)
            self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
            del soil, nsoil

        del XYZ, tree, weights, distances, indices

        return

    def _meshAdvectorSphere(self, tectonicX, tectonicY, tectonicZ, timer):
        """
        Advect spherical mesh horizontally and interpolate mesh information
        """

        # Move coordinates
        XYZ = np.zeros((self.gpoints,3))
        XYZ[:,0] = tectonicX #self.gcoords[:,0] + tectonicX*timer
        XYZ[:,1] = tectonicY #self.gcoords[:,1] + tectonicY*timer
        XYZ[:,2] = tectonicZ #self.gcoords[:,2] + tectonicZ*timer
        XYZ = XYZ[1:]

        # Get mesh variables to interpolate
        # Elevation
        elev = np.zeros(self.gpoints)
        elev.fill(-1.e8)
        elev[self.natural2local] = self.hLocal.getArray().copy()
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, elev, op=MPI.MAX)
        elev = elev[1:]

        # Erosion/deposition
        erodep = np.zeros(self.gpoints)
        erodep.fill(-1.e8)
        erodep[self.natural2local] = self.cumEDLocal.getArray().copy()
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, erodep, op=MPI.MAX)
        erodep = erodep[1:]

        # Soil thickness
        if self.Ksed > 0.:
            soil = np.zeros(self.gpoints)
            soil.fill(-1.e8)
            soil[self.natural2local] = self.HsoilLocal.getArray().copy()
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, soil, op=MPI.MAX)
            soil = soil[1:]

        # Build kd-tree
        tree = _cKDTree(XYZ)
        distances, indices = tree.query(self.lcoords, k=10)

        # Inverse weighting distance...
        weights = 1.0 / distances**2
        onIDs = np.where(distances[:,0] == 0)[0]

        nelev = np.sum(weights*elev[indices],axis=1)/np.sum(weights, axis=1)
        nerodep = np.sum(weights*erodep[indices],axis=1)/np.sum(weights, axis=1)
        if self.Ksed > 0.:
            nsoil = np.sum(weights*soil[indices],axis=1)/np.sum(weights, axis=1)

        if len(onIDs)>0:
            nelev[onIDs] = elev[indices[onIDs,0]]
            nerodep[onIDs] = erodep[indices[onIDs,0]]
            if self.Ksed > 0.:
                nsoil[onIDs] = soil[indices[onIDs,0]]

        self.hLocal.setArray(nelev)
        self.dm.localToGlobal(self.hLocal, self.hGlobal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)
        del elev, nelev

        self.cumEDLocal.setArray(nerodep)
        self.dm.localToGlobal(self.cumEDLocal, self.cumED, 1)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        del erodep, nerodep

        if self.Ksed > 0.:
            self.HsoilLocal.setArray(nsoil)
            self.dm.localToGlobal(self.HsoilLocal, self.Hsoil, 1)
            self.dm.globalToLocal(self.Hsoil, self.HsoilLocal, 1)
            del soil, nsoil

        del XYZ, tree, weights, distances, indices

        return

    def _set_DMPlex_boundary_points(self, label):
        """
        Finds the points that join the edges that have been
        marked as "boundary" faces in the DAG then sets them
        as boundaries.
        """

        self.dm.createLabel(label)
        self.dm.markBoundaryFaces(label)

        pStart, pEnd = self.dm.getDepthStratum(0) # points
        eStart, eEnd = self.dm.getDepthStratum(1) # edges
        edgeIS = self.dm.getStratumIS(label, 1)

        if edgeIS and eEnd - eStart > 0:
            edge_mask = np.logical_and(edgeIS.indices >= eStart, edgeIS.indices < eEnd)
            boundary_edges = edgeIS.indices[edge_mask]

            # Query the DAG  (directed acyclic graph) for points that join an edge
            for edge in boundary_edges:
                vertices = self.dm.getCone(edge)
                # mark the boundary points
                for vertex in vertices:
                    self.dm.setLabelValue(label, vertex, 1)
            del vertices

        edgeIS.destroy()

        return

    def destroy_DMPlex(self):
        """
        Destroy PETSc DMPlex objects and associated Petsc local/global Vectors and Matrices.
        """

        t0 = clock()
        self.bL.destroy()
        self.bG.destroy()
        self.hOld.destroy()
        self.hOldLocal.destroy()

        self.shedID.destroy()
        self.shedIDLocal.destroy()
        self.cumED.destroy()
        self.cumEDLocal.destroy()
        self.drainArea.destroy()
        self.drainAreaLocal.destroy()

        self.vSed.destroy()
        self.vSedLocal.destroy()
        self.stepED.destroy()
        self.stepEDL.destroy()
        self.areaGlobal.destroy()
        self.diffDep.destroy()
        self.diffDepLocal.destroy()
        self.hLocal.destroy()
        self.hGlobal.destroy()
        self.fillGlobal.destroy()
        self.fillLocal.destroy()

        self.Es.destroy()
        self.Eb.destroy()
        self.EsLocal.destroy()
        self.EbLocal.destroy()
        self.Hsoil.destroy()
        self.HsoilLocal.destroy()

        self.iMat.destroy()
        self.lgmap_col.destroy()
        self.lgmap_row.destroy()
        if self.Cd > 0.:
            self.Diff.destroy()
        del self.lcoords

        self.vLoc.destroy()
        self.vecG.destroy()
        self.FAL.destroy()
        self.FAG.destroy()
        self.vecL.destroy()
        self.seaG.destroy()
        self.seaL.destroy()
        self.tmpG.destroy()
        self.tmpL.destroy()
        self.vGlob.destroy()
        self.sedLoadLocal.destroy()
        self.dm.destroy()

        if MPIrank == 0 and self.verbose:
            print('Cleaning Model Dataset (%0.02f seconds)'% (clock() - t0))

        return
