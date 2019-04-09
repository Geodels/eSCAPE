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

import ruamel.yaml as yaml
import meshio

import pandas as pd
from operator import itemgetter
from scipy.interpolate import interp1d

MPIrank = PETSc.COMM_WORLD.Get_rank()
MPIsize = PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD

try: range = xrange
except: pass

bc_types = {'flat'  : 'flat',  'fixed' : 'fixed', 'slope' : 'slope'}

class ReadYaml(object):
    """
    Parsing YAML input file

    Args:
        filename: input filename (.yml YAML file)
    """
    def __init__(self, filename):

        # Check input file exists
        try:
            with open(filename) as finput:
                pass
        except IOError as exc:
            print("Unable to open file: ",filename)
            raise IOError('The input file is not found...')

        # Open YAML file
        with open(filename, 'r') as finput:
            self.input = yaml.load(finput)

        if MPIrank == 0 and 'name' in self.input.keys() and self.verbose:
            print("The following model will be run:     {}".format(self.input['name']))

        # Read simulation parameters
        self._readDomain()
        self._readTime()
        self._readSealevel()
        self._readSoil()
        self._readSPsediment()
        self._readSPbedrock()
        self._readSPdeposition()
        self._readHillslope()
        self._readOut()
        self._readRain()
        self._readTectonic()

        self.tNow = self.tStart
        self.saveTime = self.tNow

        return

    def _readDomain(self):
        """
        Read domain definition, boundary conditions and flow direction parameters.
        """

        try:
            domainDict = self.input['domain']
        except KeyError as exc:
            print("Key 'domain' is required and is missing in the input file!")
            raise KeyError('Key domain is required in the input file!')

        try:
            boundCond = domainDict['bc']
        except KeyError as exc:
            boundCond = 'slope'

        if boundCond in bc_types.keys():
            self.boundCond = bc_types[boundCond]
        else:
            raise TypeError("Boundary condition type {:s} unknown\n\
                Known types: {}".format(boundCond, bc_types.keys()))

        try:
            self.flowDir = domainDict['flowdir']
        except KeyError as exc:
            self.flowDir = 1

        try:
            self.sphere = domainDict['sphere']
        except KeyError as exc:
            self.sphere = 0

        try:
            meshFile = domainDict['filename']
        except KeyError as exc:
            print("Key 'filename' is required and is missing in the 'domain' declaration!")
            raise KeyError('Mesh file definition is not defined!')

        try:
            with open(meshFile[0]) as meshfile:
                pass
        except IOError as exc:
            print("Unable to open file: ",meshFile[0])
            raise IOError('The mesh file is not found...')

        self.meshFile = meshFile
        self.mdata = meshio.read(meshFile[0])

        try:
            self.elev = self.mdata.point_data[meshFile[1]]

        except KeyError as exc:
            print("Field name {} is missing from mesh file {}".format(meshFile[1],meshFile[0]))
            print("The following fields are available: ",self.mdata.point_data.keys())
            print("Check your mesh file fields definition...")
            raise KeyError('Field name for elevation is not defined correctly or does not exist!')

        return

    def _readTime(self):
        """
        Read simulation time declaration.
        """

        try:
            timeDict = self.input['time']
        except KeyError as exc:
            print("Key 'time' is required and is missing in the input file!")
            raise KeyError('Key time is required in the input file!')

        try:
            self.tStart = timeDict['start']
        except KeyError as exc:
            print("Key 'start' is required and is missing in the 'time' declaration!")
            raise KeyError('Simulation start time needs to be declared.')

        try:
            self.tEnd = timeDict['end']
        except KeyError as exc:
            print("Key 'end' is required and is missing in the 'time' declaration!")
            raise KeyError('Simulation end time needs to be declared.')

        try:
            self.dt = timeDict['dt']
        except KeyError as exc:
            print("Key 'dt' is required and is missing in the 'time' declaration!")
            raise KeyError('Simulation discretisation time step needs to be declared.')

        try:
            self.tout = timeDict['tout']
        except KeyError as exc:
            self.tout = self.tEnd-self.tStart
            print("Output time interval 'tout' has been set to {} years".format(self.tout))

        if self.tEnd <= self.tStart:
            raise ValueError('Simulation end/start times do not make any sense!')

        if self.tout < self.dt:
            self.tout = self.dt
            print("Output time interval was changed to {} years to match the time step dt".format(self.dt))

        return

    def _readSealevel(self):
        """
        Define sealevel evolution.
        """

        seafile = None
        seafunction = None
        sealevel = 0.
        self.seafunction = None
        try:
            seaDict = self.input['sea']
            try:
                sealevel = seaDict['position']
                try:
                    seafile = seaDict['curve']
                except KeyError as exc:
                    seafile = None
            except KeyError as exc:
                try:
                    seafile = seaDict['curve']
                except KeyError as exc:
                    seafile = None
        except KeyError as exc:
            sealevel = 0.

        if seafile is not None:
            try:
                with open(seafile) as fsea:
                    try:
                        seadata = pd.read_csv(seafile, sep=r',', engine='c',
                                                  header=None, na_filter=False,
                                                  dtype=np.float, low_memory=False)
                        pass
                    except ValueError as exc:
                        try:
                            seadata = pd.read_csv(seafile, sep=r'\s+',
                                                      engine='c', header=None,
                                                      na_filter=False, dtype=np.float,
                                                      low_memory=False)
                            pass
                        except ValueError as exc:
                            print("The sea-level file is not well formed: it should be comma or tab separated")
                            raise ValueError('Wrong formating of sea-level file.')
            except IOError as exc:
                print("Unable to open file: ",seafile)
                raise IOError('The sealevel file is not found...')

            seadata[1] += sealevel
            if seadata[0].min() > self.tStart:
                tmpS = []
                tmpS.insert(0, {0: self.tStart, 1: seadata[1].iloc[0]})
                seadata = pd.concat([pd.DataFrame(tmpS), seadata], ignore_index=True)
            if seadata[0].max() < self.tEnd:
                tmpE = []
                tmpE.insert(0, {0: self.tEnd, 1: seadata[1].iloc[-1]})
                seadata = pd.concat([seadata,pd.DataFrame(tmpE)], ignore_index=True)
            self.seafunction = interp1d(seadata[0], seadata[1]+sealevel, kind='linear')
        else:
            year = np.linspace(self.tStart, self.tEnd+self.dt, num=11, endpoint=True)
            seaval = np.full(len(year),sealevel)
            self.seafunction = interp1d(year, seaval, kind='linear')

        return

    def _readSPbedrock(self):
        """
        Read surface processes bedrock parameters.
        """

        try:
            spbDict = self.input['sp_br']
            self.mbr = 0.5
            self.nbr = 1.0
            try:
                self.Kbr = spbDict['Kbr']
            except KeyError as exc:
                print("When using the Surface Process Model definition of coefficient Kbr is required.")
                raise ValueError('Surface Process Model: Kbr coefficient not found.')
            try:
                self.crit_br = spbDict['sp_crit_br']
            except KeyError as exc:
                self.crit_br = 0.0
        except KeyError as exc:
            self.mbr = 0.5
            self.nbr = 1.0
            self.Kbr = 0.
            self.crit_br = 0.

        return

    def _readSPsediment(self):
        """
        Read surface processes sediment parameters.
        """

        try:
            spsDict = self.input['sp_sed']
            self.msed = 0.5
            self.nsed = 1.0
            try:
                self.Ksed = spsDict['Ksed']
            except KeyError as exc:
                print("When using the Surface Process Model definition of coefficient Ksed is required.")
                raise ValueError('Surface Process Model: Ksed coefficient not found.')
            try:
                self.crit_sed = spsDict['sp_crit_sed']
            except KeyError as exc:
                self.crit_sed = 0.0
        except KeyError as exc:
            self.msed = 0.5
            self.nsed = 1.0
            self.Ksed = 0.
            self.crit_sed = 0.0

        return

    def _readSPdeposition(self):
        """
        Read surface processes deposition parameters.
        """

        try:
            spdDict = self.input['sp_dep']
            try:
                self.vland = spdDict['vsL']
            except KeyError as exc:
                self.vland = 0.0
            try:
                self.frac_fine = spdDict['Ff']
            except KeyError as exc:
                self.frac_fine = 0.0
            try:
                self.phi = spdDict['phi']
            except KeyError as exc:
                self.phi = 0.0
            try:
                self.Hstar = spdDict['Hstar']
                if self.Hstar <= 0.:
                    print("Hstar needs to be a strictly positive constant!")
                    raise ValueError('Hstar needs to be a strictly positive constant!')
            except KeyError as exc:
                self.Hstar = 1.0
        except KeyError as exc:
            self.vland = 0.0
            self.frac_fine = 0.0
            self.phi = 0.0
            self.Hstar = 1.0

        return

    def _readHillslope(self):
        """
        Read hillslope parameters.
        """

        try:
            hillDict = self.input['diffusion']
            try:
                self.Cd = hillDict['hillslopeK']
            except KeyError as exc:
                print("When declaring diffusion processes, the coefficient hillslopeK is required.")
                raise ValueError('Hillslope: Cd coefficient not found.')
            try:
                self.sedimentK = hillDict['sedimentK']
            except KeyError as exc:
                self.sedimentK = 10.
        except KeyError as exc:
            self.Cd = 0.
            self.sedimentK = 10.

        return

    def _readOut(self):
        """
        Parse output directory.
        """

        try:
            outDict = self.input['output']
            try:
                self.outputDir = outDict['dir']
            except KeyError as exc:
                self.outputDir = 'output'
            try:
                self.makedir = outDict['makedir']
            except KeyError as exc:
                self.makedir = True
        except KeyError as exc:
            self.outputDir = 'output'
            self.makedir = True

        return

    def _readSoil(self):
        """
        Parse soil dataset.
        """
        try:
            soilDict = self.input['soil']
            sUniform = None
            sMap = None
            try:
                sUniform = soilDict['uniform']
            except:
                pass
            try:
                sMap = soilDict['map']
            except:
                pass

            if sMap is not None:
                if self.meshFile[0] != sMap[0]:
                    try:
                        with open(sMap[0]) as soilfile:
                            pass
                    except IOError as exc:
                        print("Unable to open soil file: ",sMap[0])
                        raise IOError('The soil file {} is not found.'.format(sMap[0]))

                    sdata = meshio.read(sMap[0])
                    soilSet = sdata.point_data
                else:
                    soilSet = self.mdata.point_data
                try:
                    soilKey = soilSet[sMap[1]]
                except KeyError as exc:
                    print("Field name {} is missing from soil file {}".format(sMap[1],sMap[0]))
                    print("The following fields are available: ",soilSet.keys())
                    print("Check your soil file fields definition...")
                    raise KeyError('Field name for soil is not defined correctly or does not exist!')


            if sMap is None and sUniform is None:
                print("For soil definition a soil thickness value (uniform) or a soil grid (map) is required.")
                raise ValueError('Soil has no value (uniform) or map (map).')

            tmpSoil = []
            if sMap is None:
                tmpSoil.insert(0, {'sUni': sUniform, 'sMap': None, 'sKey': None})
            else:
                tmpSoil.insert(0, {'sUni': None, 'sMap': sMap[0], 'sKey': sMap[1]})
            self.soildata = pd.DataFrame(tmpSoil, columns=['sUni', 'sMap', 'sKey'])

        except KeyError as exc:
            tmpSoil = []
            tmpSoil.insert(0, {'sUni': 0., 'sMap': None, 'sKey': None})
            self.soildata = pd.DataFrame(tmpSoil, columns=['sUni', 'sMap', 'sKey'])
            pass

        return

    def _readRain(self):
        """
        Parse rain forcing conditions.
        """
        try:
            rainDict = self.input['climate']
            rainSort = sorted(rainDict, key=itemgetter('start'))
            for k in range(len(rainSort)):
                rStart = None
                rUniform = None
                rMap = None
                try:
                    rStart = rainSort[k]['start']
                except:
                    print("For each climate event a start time is required.")
                    raise ValueError('Climate event {} has no parameter start'.format(k))
                try:
                    rUniform = rainSort[k]['uniform']
                except:
                    pass
                try:
                    rMap = rainSort[k]['map']
                except:
                    pass

                if rMap is not None:
                    if self.meshFile[0] != rMap[0]:
                        try:
                            with open(rMap[0]) as rainfile:
                                pass
                        except IOError as exc:
                            print("Unable to open rain file: ",rMap[0])
                            raise IOError('The rain file {} is not found for climatic event {}.'.format(rMap[0],k))

                        mdata = meshio.read(rMap[0])
                        rainSet = mdata.point_data
                    else:
                        rainSet = self.mdata.point_data
                    try:
                        rainKey = rainSet[rMap[1]]
                    except KeyError as exc:
                        print("Field name {} is missing from rain file {}".format(rMap[1],rMap[0]))
                        print("The following fields are available: ",rainSet.keys())
                        print("Check your rain file fields definition...")
                        raise KeyError('Field name for rainfall is not defined correctly or does not exist!')


                if rMap is None and rUniform is None:
                    print("For each climate event a rainfall value (uniform) or a rainfall grid (map) is required.")
                    raise ValueError('Climate event {} has no rainfall value (uniform) or a rainfall map (map).'.format(k))

                tmpRain = []
                if rMap is None:
                    tmpRain.insert(0, {'start': rStart, 'rUni': rUniform, 'rMap': None, 'rKey': None})
                else:
                    tmpRain.insert(0, {'start': rStart, 'rUni': None, 'rMap': rMap[0], 'rKey': rMap[1]})

                if k == 0:
                    raindata = pd.DataFrame(tmpRain, columns=['start', 'rUni', 'rMap', 'rKey'])
                else:
                    raindata = pd.concat([raindata,pd.DataFrame(tmpRain, columns=['start', 'rUni', 'rMap', 'rKey'])],
                                                         ignore_index=True)

            if raindata['start'][0] > self.tStart:
                tmpRain = []
                tmpRain.insert(0, {'start': self.tStart, 'rUni': 0., 'rMap': None, 'rKey': None})
                raindata = pd.concat([pd.DataFrame(tmpRain, columns=['start', 'rUni', 'rMap', 'rKey']),raindata],
                                                                              ignore_index=True)
            self.raindata = raindata

        except KeyError as exc:
            self.raindata = None
            pass

        return

    def _readTectonic(self):
        """
        Parse tectonic forcing conditions.
        """

        try:
            tecDict = self.input['tectonic']
            tecSort = sorted(tecDict, key=itemgetter('start'))
            for k in range(len(tecSort)):
                tecStart = None
                tMapX = None
                tMapY = None
                tMapZ = None
                tUniform = None
                tEnd = None
                tStep = None
                try:
                    tecStart = tecSort[k]['start']
                except:
                    print("For each tectonic event a start time is required.")
                    raise ValueError('Tectonic event {} has no parameter start'.format(k))
                try:
                    tUniform = tecSort[k]['uniform']
                except:
                    pass
                try:
                    tMapX = tecSort[k]['mapX']
                except:
                    pass
                try:
                    tMapY = tecSort[k]['mapY']
                except:
                    pass
                try:
                    tMapZ = tecSort[k]['mapZ']
                except:
                    pass
                try:
                    tStep = tecSort[k]['step']
                except:
                    pass
                try:
                    tEnd = tecSort[k]['end']
                except:
                    pass

                if tMapX is not None:
                    if self.meshFile[0] != tMapX[0]:
                        try:
                            with open(tMapX[0]) as tecfile:
                                pass
                        except IOError as exc:
                            print("Unable to open tectonic file: ",tMapX[0])
                            raise IOError('The tectonic file {} is not found for climatic event {}.'.format(tMapX[0],k))

                        mdata = meshio.read(tMapX[0])
                        tecSet = mdata.point_data
                    else:
                        tecSet = self.mdata.point_data
                    try:
                        tecKey = tecSet[tMapX[1]]
                    except KeyError as exc:
                        print("Field name {} is missing from tectonic file {}".format(tMapX[1],tMapX[0]))
                        print("The following fields are available: ",tecSet.keys())
                        print("Check your tectonic file fields definition...")
                        raise KeyError('Field name for tectonics is not defined correctly or does not exist!')


                if tMapY is not None:
                    if self.meshFile[0] != tMapY[0]:
                        try:
                            with open(tMapY[0]) as tecfile:
                                pass
                        except IOError as exc:
                            print("Unable to open tectonic file: ",tMapY[0])
                            raise IOError('The tectonic file {} is not found for climatic event {}.'.format(tMapY[0],k))

                        mdata = meshio.read(tMapY[0])
                        tecSet = mdata.point_data
                    else:
                        tecSet = self.mdata.point_data
                    try:
                        tecKey = tecSet[tMapY[1]]
                    except KeyError as exc:
                        print("Field name {} is missing from tectonic file {}".format(tMapY[1],tMapY[0]))
                        print("The following fields are available: ",tecSet.keys())
                        print("Check your tectonic file fields definition...")
                        raise KeyError('Field name for tectonics is not defined correctly or does not exist!')

                if tMapZ is not None:
                    if self.meshFile[0] != tMapZ[0]:
                        try:
                            with open(tMapZ[0]) as tecfile:
                                pass
                        except IOError as exc:
                            print("Unable to open tectonic file: ",tMapZ[0])
                            raise IOError('The tectonic file {} is not found for climatic event {}.'.format(tMapZ[0],k))

                        mdata = meshio.read(tMapZ[0])
                        tecSet = mdata.point_data
                    else:
                        tecSet = self.mdata.point_data
                    try:
                        tecKey = tecSet[tMapZ[1]]
                    except KeyError as exc:
                        print("Field name {} is missing from tectonic file {}".format(tMapZ[1],tMapZ[0]))
                        print("The following fields are available: ",tecSet.keys())
                        print("Check your tectonic file fields definition...")
                        raise KeyError('Field name for tectonics is not defined correctly or does not exist!')


                if tMapX is None and tMapY is None and tMapZ is None and tUniform is None:
                    print("For each tectonic event a tectonic value (uniform) or a tectonic grid (map) is required.")
                    raise ValueError('Tectonic event {} has no tectonic value (uniform) or a tectonic map (map).'.format(k))

                tmpTec = []
                if tUniform is not None:
                    tmpTec.insert(0, {'start': tecStart, 'tUni': tUniform, 'tMapX': None, 'tMapY': None, 'tMapZ': None,
                                      'tKeyX': None, 'tKeyY': None, 'tKeyZ': None})
                else:
                    if tMapX is None:
                        tMapX = [None] * 2
                    if tMapY is None:
                        tMapY = [None] * 2
                    if tMapZ is None:
                        tMapZ = [None] * 2
                    tmpTec.insert(0, {'start': tecStart, 'tUni': None, 'tMapX': tMapX[0], 'tMapY': tMapY[0], 'tMapZ': tMapZ[0],
                                      'tKeyX': tMapX[1], 'tKeyY': tMapY[1], 'tKeyZ': tMapZ[1]})

                if k == 0:
                    tecdata = pd.DataFrame(tmpTec, columns=['start', 'tUni', 'tMapX', 'tMapY','tMapZ', 'tKeyX', 'tKeyY', 'tKeyZ'])
                else:
                    tecdata = pd.concat([tecdata,pd.DataFrame(tmpTec, columns=['start', 'tUni', 'tMapX', 'tMapY',
                                         'tMapZ', 'tKeyX', 'tKeyY', 'tKeyZ'])], ignore_index=True)

                if tStep is not None:
                    if tEnd is not None:
                        tectime = tecStart+tStep
                        while tectime<tEnd:
                            tmpTec = []
                            tmpTec.insert(0, {'start': tectime, 'tUni': None, 'tMapX': tMapX[0], 'tMapY': tMapY[0], 'tMapZ': tMapZ[0],
                                              'tKeyX': tMapX[1], 'tKeyY': tMapY[1], 'tKeyZ': tMapZ[1]})
                            tecdata = pd.concat([tecdata,pd.DataFrame(tmpTec, columns=['start', 'tUni', 'tMapX', 'tMapY','tMapZ',
                                                'tKeyX', 'tKeyY', 'tKeyZ'])], ignore_index=True)
                            tectime = tectime+tStep

            if tecdata['start'][0] > self.tStart:
                tmpTec = []
                tmpTec.insert(0, {'start': self.tStart, 'tUni': 0., 'tMapX': None, 'tMapY': None, 'tMapZ': None,
                                  'tKeyX': None, 'tKeyY': None, 'tKeyZ': None})
                tecdata = pd.concat([pd.DataFrame(tmpTec, columns=['start', 'tUni', 'tMapX', 'tMapY','tMapZ',
                                    'tKeyX', 'tKeyY', 'tKeyZ']),tecdata], ignore_index=True)
            self.tecdata = tecdata
            
        except KeyError as exc:
            self.tecdata = None
            pass

        return
