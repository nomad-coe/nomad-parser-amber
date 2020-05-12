# Copyright 2017-2018 Berk Onat, Fawzi Mohamed
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from builtins import map
from builtins import range
from builtins import object
import numpy as np
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import mainFunction, ParsingContext
from nomadcore.simple_parser import SimpleMatcher as SM
from .AMBERDictionary import set_excludeList, set_includeList, get_updateDictionary, getList_MetaStrInDict, getDict_MetaStrInDict
from nomadcore.smart_parser.SmartParserDictionary import isMetaStrInDict
from contextlib import contextmanager
from . import AMBERCommon as AmberC
from . import trajectory_reader as TrajRead
import logging
import os
import re
import sys


############################################################
# This is the parser for the main file of AMBER (mdout).
############################################################

LOGGER = logging.getLogger("nomad.AMBERParser")

#PRINTABLE = re.compile(r"\W+")

@contextmanager
def open_section(parser, name):
    gid = parser.openSection(name)
    yield gid
    parser.closeSection(name, gid)

class AMBERParser(AmberC.AMBERParserBase):
    """Context for parsing AMBER main file.

    This class keeps tracks of several AMBER settings to adjust the parsing to them.
    The onClose_ functions allow processing and writing of cached values after a section is closed.
    They take the following arguments:
        backend: Class that takes care of writing and caching of metadata.
        gIndex: Index of the section that is closed.
        section: The cached values and sections that were found in the section that is closed.
    """
    def __init__(self):
        # dictionary of energy values, which are tracked between SCF iterations and written after convergence
        self.totalEnergyList = {
                                'energy_electrostatic': None,
                                'energy_total_T0_per_atom': None,
                                'energy_free_per_atom': None,
                               }
        AmberC.AMBERParserBase.__init__(
            self, re_program_name=re.compile(r"\s*Amber$"))

        self.cachingLevelForMetaName = {
                               'x_amber_trajectory_file_detect': CachingLevel.Cache,
                               'x_amber_geometry_optimization_cdetect': CachingLevel.Cache,
                               'x_amber_mdin_finline': CachingLevel.Ignore,
                               'x_amber_mdin_wt': CachingLevel.Ignore,
                               'x_amber_section_input_output_files': CachingLevel.Ignore,
                               'x_amber_single_configuration_calculation_detect': CachingLevel.Cache,
                              }
        for name in self.metaInfoEnv.infoKinds:
            metaInfo = self.metaInfoEnv.infoKinds[name]
            if (name.startswith('x_amber_mdin_') and
                metaInfo.kindStr == "type_document_content" and
                ("x_amber_mdin_method" in metaInfo.superNames or
                 "x_amber_mdin_run" in metaInfo.superNames or
                 "x_amber_mdin_system" in metaInfo.superNames) or
                name.startswith('x_amber_parm_') and
                metaInfo.kindStr == "type_document_content" and
                ("x_amber_mdin_method" in metaInfo.superNames or
                 "x_amber_mdin_run" in metaInfo.superNames or
                 "x_amber_mdin_system" in metaInfo.superNames) or
                name.startswith('x_amber_mdin_file_') and
                metaInfo.kindStr == "type_document_content" and
                ("x_amber_section_input_output_files" in metaInfo.superNames) or
                name.startswith('x_amber_mdout_') and
                metaInfo.kindStr == "type_document_content" and
                ("x_amber_mdout_method" in metaInfo.superNames or
                 "x_amber_mdout_system" in metaInfo.superNames or
                 "x_amber_mdout_single_configuration_calculation" in metaInfo.superNames)
                or name.startswith('section_single_configuration_calculation')
               ):
                self.cachingLevelForMetaName[name] = CachingLevel.Cache
            if name in self.extraDict.keys():
                self.cachingLevelForMetaName[name] = CachingLevel.Ignore


    def initialize_values(self):
        """Initializes the values of certain variables.

        This allows a consistent setting and resetting of the variables,
        when the parsing starts and when a section_run closes.
        """
        self.secMethodGIndex = None
        self.secSystemGIndex = None
        self.secTopologyGIndex = None
        self.secSamplingGIndex = None
        self.secSingleGIndex = None
        self.secVDWGIndex = None
        self.secAtomType = None
        self.inputMethodIndex = None
        self.mainMethodIndex = None
        self.mainCalcIndex = None
        self.MD = True
        self.topologyDict = None
        self.topologyTable = None
        self.topologyBonds = None
        self.topology = None
        self.topologyFormat = None
        self.topologyFile = None
        self.trajectory = None
        self.trajectoryFormat = None
        self.trajectoryFile = None
        self.readChunk = 300
        self.boxlengths = None
        self.latticevectors = None
        self.atompositions = None
        # start with -1 since zeroth iteration is the initialization
        self.MDiter = -1
        self.singleConfCalcs = []
        self.minConverged = None
        self.parsedLogFile = False
        self.LogSuperContext = None
        self.forces_raw = []
        self.lastCalculationGIndex = None
        self.logFileName = None
        self.lastfInLine = None
        self.lastfInMatcher = None
        self.secOpen = open_section
        self.superP = self.parser.backend.superBackend

    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts.

        Get compiled parser, filename and metadata.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        self.fName = fInName
        # save metadata
        self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv
        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()

    def peekline(self, parser):
        pos = parser.fIn.fIn.tell()
        line = parser.fIn.fIn.readline()
        parser.fIn.fIn.seek(pos)
        return line

    def topologyFileHandler(self, fileItem):
        if 'topology' in self.fileDict[fileItem].infoPurpose:
            topofile = self.fileDict[fileItem].fileName
            self.topologyFile = topofile
            self.trajectory = TrajRead.TrajectoryReader()
            self.trajectory.topofile = topofile
            self.topology = self.trajectory.load_topology()
            if self.topology is None:
                for fileFormat in self.fileDict[fileItem].fileFormat:
                    self.trajectory = TrajRead.TrajectoryReader()
                    self.trajectory.topofile = topofile
                    self.trajectory.topoformat = fileFormat
                    self.topology = self.trajectory.load_topology()
                    if self.topology is not None:
                        return fileFormat

    def trajectoryFileHandler(self, fileItem, topoformat):
        if 'trajectory' in self.fileDict[fileItem].infoPurpose:
            trajfile = self.fileDict[fileItem].fileName
            self.trajectoryFile = trajfile
            self.trajectory = TrajRead.TrajectoryReader()
            self.trajectory.trajfile = trajfile
            self.trajectory.topofile = self.topologyFile
            self.trajectory.topoformat = topoformat
            traj_loaded = self.trajectory.load()
            self.trajectory.trajchunk=self.readChunk
            self.atompositions = self.trajectory.iread()
            if self.atompositions is None:
                for fileFormat in self.fileDict[fileItem].fileFormat:
                    self.trajectory = TrajRead.TrajectoryReader()
                    self.trajectory.trajfile = trajfile
                    self.trajectory.topofile = self.topologyFile
                    self.trajectory.topoformat = topoformat
                    self.trajectory.trajformat = fileFormat
                    self.trajectory.trajchunk=self.readChunk
                    traj_loaded = self.trajectory.load()
                    self.atompositions = self.trajectory.iread()
                    if self.atompositions is not None:
                        self.topology = self.trajectory.get_topology()
                        self.MDiter += 1
                        return fileFormat
            else:
                self.MDiter += 1

    def initializeFileHandlers(self):
        # Files will be loaded using their extensions initially.
        # If this fails, the fileFormat lists will be used in loading process.
        topoformat = None
        trajformat = None
        for fileItem in self.fileDict:
            if (self.fileDict[fileItem].fileSupplied and
                self.fileDict[fileItem].activeInfo):
                # First, check topology file
                topoformat = self.topologyFileHandler(fileItem)
                # Second, check trajectory file
                trajformat = self.trajectoryFileHandler(fileItem, topoformat)
        if self.topology:
            self.topologyTable, self.topologyBonds = self.topology.to_dataframe()
            self.topologyDict = self.topologyTable.to_dict(orient='list')
        self.topologyFormat = topoformat
        self.trajectoryFormat = trajformat
        #if trajformat or topoformat:
        #    return True
        #else:
        #    return False

    def parameter_file_name(self, itemdict):
        """ Function to generate data for parameter files list
        """
        working_dir_name = os.path.dirname(os.path.abspath(self.fName))
        parmmeta = isMetaStrInDict("PARM",self.fileDict)
        filename = []
        if parmmeta is not None:
            if self.fileDict[parmmeta].value is not None:
                filename.append(self.fileDict[parmmeta].value)
        if filename:
            return False, filename, itemdict
        else:
            return False, None, itemdict

    def topology_to_dictionary(self):
        """ This function generates self.topologyDict dictionary
            and converts/stores all the topology information
            according to the meta info definitions
        """
        newDict = {}
        return newDict

    def topology_num_topo_mol(self, itemdict):
        """ Function to generate data for number_of_topology_molecules
        """
        if self.topology:
            residueList = self.topologyDict["resSeq"]
            return False, residueList[len(residueList)-1], itemdict
        else:
            return False, None, itemdict

    def topology_system_name(self, itemdict):
        """ Function to generate data for system_name
        """
        system_name = self.fName.split('.')[-1]
        if self.topology:
            residueList = ','.join(list(set(self.topologyDict["resName"])))
        return False, system_name, itemdict

    def topology_atom_to_mol(self, itemdict):
        """ Function to generate data for atom_to_molecule
        """
        if self.topology:
            residueList = self.topologyDict["resSeq"]
            atomIndex = np.arange(len(residueList))
            atom_to_residue = np.zeros((len(residueList), 2), dtype=int)
            atom_to_residue[:,0] = atomIndex+1
            atom_to_residue[:,1] = np.array(residueList)+1
            return False, atom_to_residue, itemdict
        else:
            return False, None, itemdict

    def topology_atom_type_and_interactions(self, backend, gIndex):
        """ Function to generate data for atom_to_molecule
        """
        sO = open_section
        supB = backend.superBackend
        if self.topology:
            atom_type_list=list(set(self.topologyDict["name"]))
            atom_type_dict = {}
            massesDict = {}
            elementDict = {}
            radiusDict = {}
            for ielm in range(len(atom_type_list)):
                elm = atom_type_list[ielm]
                atom_type_dict.update({elm : ielm+1})
                for atom in self.topology.atoms:
                    if elm == atom.name:
                        massesDict.update({atom.name : atom.element.mass})
                        elementDict.update({atom.name : atom.element.symbol})
                        radiusDict.update({atom.name : atom.element.radius})

            massesList = list(massesDict.values())
            elementList = list(elementDict.values())
            radiusList = list(radiusDict.values())

            interNum = 0
            interDict = {}
            interTypeDict = {}
            for ielm in range(len(atom_type_list)-1):
                for jelm in range(ielm+1, len(atom_type_list)):
                    aelm = atom_type_list[ielm]
                    belm = atom_type_list[jelm]
                    bondList = []
                    typeList = []
                    bondid = 0
                    noninter = True
                    for bond in self.topology.bonds:
                        molname1, molatom1 = str(bond[0]).split('-')
                        molname2, molatom2 = str(bond[1]).split('-')
                        if((aelm == str(molatom1) and belm == str(molatom2)) or
                            (aelm == str(molatom2) and belm == str(molatom1))):
                            bondList.append(list(self.topologyBonds[bondid]))
                            interDict.update({interNum : bondList})
                            if noninter:
                                noninter = False
                                typeList.extend(list([
                                    atom_type_dict[aelm],
                                    atom_type_dict[belm]
                                    ]))
                                interTypeDict.update({interNum : typeList})
                                interNum += 1
                        bondid += 1

            for elm in range(len(atom_type_list)):
                with sO(supB, 'section_atom_type'):
                    ### !!! ----------------------------------------------- !!!
                    ### !!! Need to check the definition of atom type name.
                    ### !!! Which one is better? C1, CA, C:1 or C !!!
                    ### !!! ----------------------------------------------- !!!
                    # Atom name? Atom element?
                    supB.addValue('atom_type_name', str(atom_type_list[elm]) + ':' + str(elm+1))
                    # Atomic mass
                    supB.addValue('atom_type_mass', massesList[elm])
                    # Atomic element/symbol
                    supB.addValue('x_amber_atom_type_element', elementList[elm])
                    # Atomic van der Waals radius
                    supB.addValue('x_amber_atom_type_radius', radiusList[elm])
                    # Atomic charge
                    #self.superP.addValue('atom_type_charge', atom_charge)
                    pass

            for inum in range(interNum):
                with sO(supB, 'section_interaction'):
                    # atom indexes of bound pairs for a specific atom type
                    supB.addArrayValues('interaction_atoms', np.asarray(interDict[inum]))
                    # number of bonds for this type
                    supB.addValue('number_of_interactions', len(interDict[inum]))
                    # number of atoms involved (2 for covalent bonds)
                    supB.addValue('number_of_atoms_per_interaction', len(interDict[inum][0]))

                    #if bondFunctional:
                    #    self.superP.addValue('interaction_kind', bondFunctional)  # functional form of the interaction

                    # this points to the relative section_atom_type
                    supB.addArrayValues('x_amber_interaction_atom_to_atom_type_ref',
                            np.asarray(interTypeDict[inum]))
                    # interaction parameters for the functional
                    #self.superP.addValue('interaction_parameters', bondParameters)


#            for i in range(len(moleculeTypeInfo)):
#
#                with o(p, 'section_molecule_type'):
#                    # gindex = 0
#                    p.addValue('molecule_type_name', 'molecule'+'_'+str(moleculeTypeInfo[i][0]))
#                    p.addValue('number_of_atoms_in_molecule', len(moleculeTypeInfo[i][1]))
#
#                    p.addArrayValues('atom_in_molecule_to_atom_type_ref', np.asarray([x-1 for x in moleculeTypeInfo[i][1]]))
#
#
#            residueList = self.topologyDict["resSeq"]
#            atomIndex = np.arange(len(residueList))
#            atom_to_residue = np.zeros((len(residueList), 2), dtype=int)
#            atom_to_residue[:,0] = atomIndex
#            atom_to_residue[:,1] = residueList
#            return atom_to_residue
#        else:
#            return None

#    mol_to_mol_bond = []
#    mol_to_mol_bond_num = []
#    atom_in_mol_bond = []
#    atom_in_mol_bond_num = []
#    index=0
#    print(len([res for res in mytopology.bonds]))
#    for res in mytopology.bonds:
#        molname1, molatom1 = str(res[0]).split('-')
#        molname2, molatom2 = str(res[1]).split('-')
#        if molname1 in molname2:
#            atom_in_mol_bond.append(res)
#            atom_in_mol_bond_num.append(bonds[index])
#        else:
#            mol_to_mol_bond.append(res)
#            mol_to_mol_bond_num.append(bonds[index])
#        index += 1
#    print(mol_to_mol_bond)
#    print(np.array(mol_to_mol_bond_num))

    def onClose_section_run(self, backend, gIndex, section):
        """Trigger called when section_run is closed.

        Write convergence of geometry optimization.
        Write the keywords from control parametres and the Amber output from the parsed log output, which belong to settings_run.
        Write the last occurrence of a keyword/setting, i.e. [-1], since Amber uses the last occurrence of a keyword.
        Variables are reset to ensure clean start for new run.

        ATTENTION
        backend.superBackend is used here instead of only the backend to write the JSON values,
        since this allows to bybass the caching setting which was used to collect the values for processing.
        However, this also bypasses the checking of validity of the metadata name by the backend.
        The scala part will check the validity nevertheless. (Good to know :)
        """
        # write trajectory
        #valuesDict = section.simpleValues

        frameSequenceGIndex = backend.openSection("section_frame_sequence")
        self.metaStorage.updateBackend(backend,
                startsection=['section_frame_sequence'],
                autoopenclose=False)
        backend.addValue("frame_sequence_to_sampling_ref", self.secSamplingGIndex)
        backend.addArrayValues("frame_sequence_local_frames_ref", np.asarray(self.singleConfCalcs))
        backend.closeSection("section_frame_sequence", frameSequenceGIndex)

        # reset all variables
        self.initialize_values()

    def onClose_x_amber_section_input_output_files(self, backend, gIndex, section):
        """Trigger called when x_amber_section_input_output_files is closed.

        Determine whether topology, trajectory and input coordinate files are
        supplied to the parser

        Initiates topology and trajectory file handles.

        Captures topology, atomic positions, atom labels, lattice vectors and
        stores them before section_system and
        section_single_configuration_calculation are encountered.
        """
        # Checking whether topology, input
        # coordinates and trajectory files exist
        atLeastOneFileExist = False
        working_dir_name = os.path.dirname(os.path.abspath(self.fName))
        for k,v in section.simpleValues.items():
            if k.startswith('x_amber_mdin_file'):
                file_name = os.path.normpath(os.path.join(working_dir_name, v[-1]))
                self.fileDict[k].fileSupplied = os.path.isfile(file_name)
                if self.fileDict[k].fileSupplied:
                    self.fileDict[k].fileName = file_name
                    if self.fileDict[k].activeInfo:
                        self.fileDict[k].value = v[-1]
                        atLeastOneFileExist = True
        if atLeastOneFileExist:
            self.initializeFileHandlers()

    def onOpen_section_method(self, backend, gIndex, section):
        # keep track of the latest method section
        self.secMethodGIndex = gIndex
        if self.inputMethodIndex is None:
            self.inputMethodIndex = gIndex
        else:
            backend.openNonOverlappingSection("section_method_to_method_refs")
            backend.addValue("method_to_method_kind", "core_settings")
            backend.addValue("method_to_method_ref", self.inputMethodIndex)
            backend.closeNonOverlappingSection("section_method_to_method_refs")
        if self.mainMethodIndex is None:
            self.mainMethodIndex = gIndex

    def onClose_section_method(self, backend, gIndex, section):
        """Trigger called when section_method is closed.
        """
        # input method
        pass

    def onOpen_section_sampling_method(self, backend, gIndex, section):
        # keep track of the latest sampling description section
        self.secSamplingGIndex = gIndex

    def onClose_section_sampling_method(self, backend, gIndex, section):
        """Trigger called when section_sampling_method is closed.

        Writes sampling method details for minimization and molecular dynamics.
        """
        # check control keywords were found throguh dictionary support
        section_sampling_Dict = get_updateDictionary(self, 'sampling')
        updateDict = {
            'startSection' : [['section_sampling_method']],
            'dictionary' : section_sampling_Dict
            }
        #self.secSamplingGIndex = backend.openSection("section_sampling_method")
        self.metaStorage.update(updateDict)
        self.metaStorage.updateBackend(backend,
                startsection=['section_sampling_method'],
                autoopenclose=False)
        self.secRestrictGIndex = backend.openSection("section_restricted_uri")
        section_restrictions_Dict = get_updateDictionary(self, 'restrictions')
        updateDict = {
            'startSection' : [['section_restricted_uri']],
            'dictionary' : section_restrictions_Dict
            }
        self.metaStorageRestrict.update(updateDict)
        self.metaStorageRestrict.updateBackend(backend,
                startsection=['section_restricted_uri'],
                autoopenclose=False)
        backend.closeSection("section_restricted_uri", self.secRestrictGIndex)

    def onOpen_section_topology(self, backend, gIndex, section):
        # keep track of the latest topology description section
        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            self.secTopologyGIndex = backend.superBackend.openSection("section_topology")
        else:
            self.secTopologyGIndex = gIndex

    def onClose_section_topology(self, backend, gIndex, section):
        """Trigger called when section_topology is closed.
        """
        section_topology_Dict = get_updateDictionary(self, 'topology')
        updateDict = {
            'startSection' : [
                ['section_topology']],
            'muteSections' : [['section_sampling_method']],
            'dictionary' : section_topology_Dict
            }
        self.metaStorage.update(updateDict)
        self.metaStorage.updateBackend(backend.superBackend,
                startsection=['section_topology'],
                autoopenclose=False)
        self.topology_atom_type_and_interactions(backend, gIndex)

        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            backend.superBackend.closeSection("section_topology", self.secTopologyGIndex)

    def onOpen_section_system(self, backend, gIndex, section):
        # keep track of the latest system description section
        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            self.secSystemGIndex = backend.superBackend.openSection("section_system")
        else:
            self.secSystemGIndex = gIndex

    def onClose_section_system(self, backend, gIndex, section):
        """Trigger called when section_system is closed.

        Writes atomic positions, atom labels and lattice vectors.
        """
        # Our special recipe for the confused backend and more
        # See our other recipes at the back of this parser: MetaInfoStorage! :)
        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            SloppyBackend = backend.superBackend
        else:
            SloppyBackend = backend

        if self.topology:
            if (self.secTopologyGIndex is None or
                (self.secTopologyGIndex == -1 or
                self.secTopologyGIndex == "-1")):
                self.onOpen_section_topology(backend, None, None)
                self.onClose_section_topology(backend, None, None)
            SloppyBackend.addValue("topology_ref", self.secTopologyGIndex)

        if self.atompositions is not None:
            self.trajRefSingleConfigurationCalculation = gIndex
            unit_cell = np.asarray(self.metaStorage.convertUnits(
                self.atompositions.unitcell_vectors[0], "Angstrom", self.unitDict))
            SloppyBackend.addArrayValues('simulation_cell', unit_cell)
            SloppyBackend.addArrayValues('lattice_vectors', unit_cell)
            SloppyBackend.addArrayValues('atom_positions', np.transpose(np.asarray(
                self.metaStorage.convertUnits(self.atompositions.xyz[0], "Angstrom", self.unitDict))))
            if self.topology is not None:
                atom_labels = self.topologyDict["element"]
                SloppyBackend.addArrayValues('atom_labels', np.asarray(atom_labels))
                pass

            # Read the next step at trajectory in advance
            # If iread returns None, it will be the last step
            try:
                self.atompositions = self.trajectory.iread()
                self.MDiter += 1
            except AttributeError:
                self.atompositions = None
                pass

            #if atom_vel:
            # need to transpose array since its shape is [number_of_atoms,3] in the metadata
            #backend.addArrayValues('atom_velocities', np.transpose(np.asarray(atom_vel)))

        if (gIndex is None or gIndex == -1 or gIndex == "-1"):
            SloppyBackend.closeSection("section_system", self.secSystemGIndex)

    def onOpen_section_single_configuration_calculation(self, backend, gIndex, section):
        # write the references to section_method and section_system
        backend.addValue('single_configuration_to_calculation_method_ref', self.secMethodGIndex)
        backend.addValue('single_configuration_calculation_to_system_ref', self.secSystemGIndex)
        self.singleConfCalcs.append(gIndex)
        self.secSingleGIndex = backend.superBackend.openSection("section_single_configuration_calculation")

    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        """Trigger called when section_single_configuration_calculation is closed.

        Write number of steps in MD or Minimization.
        Check for convergence of geometry optimization.
        Write energy values at MD and with error in Minimization.
        Write reference to section_method and section_system
        """

        self.lastCalculationGIndex = gIndex

        section_frameseq_Dict = get_updateDictionary(self, 'frameseq')
        updateFrameDict = {
            'startSection' : [
                ['section_frame_sequence']],
            'muteSections' : [['section_sampling_method']],
            'dictionary' : section_frameseq_Dict
            }
        self.metaStorage.update(updateFrameDict)
        section_singlevdw_Dict = get_updateDictionary(self, 'singlevdw')
        updateDictVDW = {
            'startSection' : [
                ['section_energy_van_der_Waals']],
            'muteSections' : [['section_sampling_method']],
            'dictionary' : section_singlevdw_Dict
            }
        self.secVDWGIndex = backend.superBackend.openSection("section_energy_van_der_Waals")
        self.metaStorage.update(updateDictVDW)
        self.metaStorage.updateBackend(backend.superBackend,
                startsection=['section_energy_van_der_Waals'],
                autoopenclose=False)
        backend.superBackend.closeSection("section_energy_van_der_Waals", self.secVDWGIndex)
        section_singlecalc_Dict = get_updateDictionary(self, 'singleconfcalc')
        updateDict = {
            'startSection' : [
                ['section_single_configuration_calculation']],
            'muteSections' : [['section_sampling_method']],
            'dictionary' : section_singlecalc_Dict
            }
        self.metaStorage.update(updateDict)
        self.metaStorage.updateBackend(backend.superBackend,
                startsection=['section_single_configuration_calculation'],
                autoopenclose=False)
        if(self.topology is not None or
           self.atompositions is not None):
            self.onOpen_section_system(backend, None, None)
            self.onClose_section_system(backend, None, None)
        backend.superBackend.closeSection("section_single_configuration_calculation", self.secSingleGIndex)

        # write number of Minimization steps
        #backend.addValue('number_of_steps', self.minStepNr)
#        # write Minimization convergence and reset
#        backend.addValue('single_configuration_calculation_converged', self.minConvergence)
#        self.minConvergence = False
#        # start with -1 since zeroth iteration is the initialization
#        self.minStepNr = -1
#        # write converged energy/thermodynamic values
#        # write forces
#        forces_free = []
#        for i in ['x', 'y', 'z']:
#            fi = section['x_amber_atom_forces_free_' + i]
#            if fi is not None:
#                forces_free.append(fi)
#        if forces_free:
#            # need to transpose array since its shape is [number_of_atoms,3] in the metadata
#            backend.addArrayValues('atom_forces_free', np.transpose(np.asarray(forces_free)))
#        if self.forces_raw:
#            # need to transpose array since its shape is [number_of_atoms,3] in the metadata
#            backend.addArrayValues('atom_forces_free_raw', np.transpose(np.asarray(self.forces_raw)))

    def setStartingPointCalculation(self, parser):
        backend = parser.backend
        backend.openSection('section_calculation_to_calculation_refs')
        if self.lastCalculationGIndex:
            backend.addValue('calculation_to_calculation_ref', self.lastCalculationGIndex)
#        backend.addValue('calculation_to_calculation_kind', 'pertubative GW')
            backend.closeSection('section_calculation_to_calculation_refs')
        return None

    def check_namelist_store(self, parser, lastLine, stopOnMatchRe, quitOnMatchRe,
            metaNameStart, matchNameList, matchNameDict, onlyCaseSensitive, stopOnFirstLine):
        stopOnMatch = False
        if stopOnMatchRe.findall(lastLine):
            stopOnMatch = True
            if self.firstLine==0:
                if stopOnFirstLine:
                    stopOnMatch = True
                else:
                    stopOnMatch = False
        if quitOnMatchRe is not None:
            if quitOnMatchRe.findall(lastLine):
                stopOnMatch = True
        if stopOnMatch:
            return True
        else:
            # If there is at least one namelist in the line,
            # search for all others in the dictionary.
            if self.MD is not True:
                newLine = parser.fIn.readline()
                lastLine = ' = '.join([ "%s" % str(line) for line in zip(lastLine, newLine)])
            for cName, key in getDict_MetaStrInDict(matchNameDict).items():
                reDict={key:value for value in
                        re.compile(r"(?:\s%s|^%s|,%s)\s*=\s*(?:'|\")?(?P<%s>[\-+0-9.a-zA-Z:]+)(?:'|\")?\s*,?"
                        % (cName, cName, cName, key)).findall(lastLine)}
                if onlyCaseSensitive is not True:
                    reDict.update({key:value for value in
                                   re.compile(r"(?:\s%s|^%s|,%s)\s*=\s*(?:'|\")?(?P<%s>[\-+0-9.a-zA-Z:]+)(?:'|\")?\s*,?"
                                   % (cName.upper(), cName.upper(), cName.upper(), key)).findall(lastLine)})
                if reDict:
                    for k,v in reDict.items():
                        if k == key:
                            if k in list(parser.lastMatch.keys()):
                                parser.lastMatch[k]=v
                            else:
                                matchNameDict[k].value=v
                                matchNameDict[k].activeInfo=True
            return False

    def adHoc_read_namelist_stop_parsing(self, parser, stopOnMatchStr, quitOnMatchStr,
            metaNameStart, matchNameList, matchNameDict, onlyCaseSensitive, stopOnFirstLine):
        lastLine = parser.fIn.fInLine
        self.firstLine = 0
        # Check the captured line has Fortran namelist variables and store them.
        # Continue search and store until the line matches with stopOnMatch.
        stopOnMatchRe = re.compile(stopOnMatchStr)
        quitOnMatchRe = None
        if quitOnMatchStr is not None:
            quitOnMatchRe = re.compile(quitOnMatchStr)
        if self.check_namelist_store(parser, lastLine,
                stopOnMatchRe, quitOnMatchRe,
                metaNameStart, matchNameList,
                matchNameDict, onlyCaseSensitive,
                stopOnFirstLine) is not True:
            while True:
                lastLine = self.peekline(parser)
                self.firstLine += 1
                if not lastLine:
                    break
                else:
                    # Matched with stopOnMatch. Discarding the line and return SimpleMatcher context.
                    # Can this line be discarded since it is the end of line for input control
                    # variables or end of namelist ?
                    if self.check_namelist_store(parser, lastLine,
                            stopOnMatchRe, quitOnMatchRe,
                            metaNameStart, matchNameList,
                            matchNameDict, onlyCaseSensitive,
                            stopOnFirstLine):
                        break
                    else:
                        lastLine = parser.fIn.readline()

    def build_mdinKeywordsSimpleMatchers(self):
        cntrlNameList=getList_MetaStrInDict(self.cntrlDict)
        ewaldNameList=getList_MetaStrInDict(self.ewaldDict)
        qmmmNameList=getList_MetaStrInDict(self.qmmmDict)
        wtNameList=getList_MetaStrInDict(self.wtDict)
        return [
            SM(name="cntrl",
               startReStr=r"\s*&cntrl",
               endReStr=r"(?:(?:^/|\s*/)|&end|&END)\s*$",
               subMatchers=[
                   SM(startReStr=(r"\s*(?:" +
                      '|'.join(["%s" % (cName) for cName in cntrlNameList]) +
                      "AMBER)\s*=\s*(?:'|\")?(?P<x_amber_mdin_finline>[\-+0-9.:a-zA-Z]+)(?:'|\")?\s*,?"),
                      coverageIgnore=True,
                      adHoc=lambda p:
                      self.adHoc_read_namelist_stop_parsing(p,
                      stopOnMatchStr=r"(?:(?:^/|\s*/)|&end|&END)\s*$",
                      quitOnMatchStr=None,
                      metaNameStart="x_amber_mdin_",
                      matchNameList=cntrlNameList,
                      matchNameDict=self.cntrlDict,
                      onlyCaseSensitive=True,
                      stopOnFirstLine=True)
                      )
               ]),
            SM(name="wt",
                startReStr=r"(?:\s*&wt|^&wt)",
                endReStr=r"\s*&wt\s*(?:type|TYPE)\s*=\s*(?:\'|\")(?:end|END)(?:\'|\")\s*\/\s*",
                forwardMatch=True,
                subMatchers=[
                    SM(startReStr=r"(?:\s*&wt|^&wt)\s*(?:type|TYPE)\s*=" +
                       "(?:'|\")?(?P<x_amber_mdin_wt>[-+0-9.:a-zA-Z]+)(?:'|\")?",
                       repeats=True,
                       adHoc=lambda p:
                       self.adHoc_read_namelist_stop_parsing(p,
                       stopOnMatchStr=r"\s*&wt\s*(?:type|TYPE)\s*=\s*(?:\'|\")(?:end|END)(?:\'|\")\s*\/\s*",
                       quitOnMatchStr=None,
                       metaNameStart="x_amber_mdin_",
                       matchNameList=wtNameList,
                       matchNameDict=self.wtDict,
                       onlyCaseSensitive=False,
                       stopOnFirstLine=True)
                       )
                ]),
            SM(name="ewald",
               startReStr=r"\s*&ewald",
               endReStr=r"(?:(?:^/|\s*/)|&end|&END)\s*$",
               subMatchers=[
                   SM(startReStr=(r"\s*(?:" +
                      '|'.join(["%s" % (cName) for cName in ewaldNameList]) +
                      "AMBER)\s*=\s*(?:'|\")?(?P<x_amber_mdin_finline>[\-+0-9.:a-zA-Z]+)(?:'|\")?\s*,?"),
                      coverageIgnore=True,
                      adHoc=lambda p:
                      self.adHoc_read_namelist_stop_parsing(p,
                      stopOnMatchStr=r"(?:(?:^/|\s*/)|&end|&END)\s*$",
                      quitOnMatchStr=None,
                      metaNameStart="x_amber_mdin_",
                      matchNameList=ewaldNameList,
                      matchNameDict=self.ewaldDict,
                      onlyCaseSensitive=True,
                      stopOnFirstLine=True)
                      )
               ]),
            SM(name="qmmm",
               startReStr=r"\s*&qmmm",
               endReStr=r"(?:(?:^/|\s*/)|&end|&END)\s*$",
               subMatchers=[
                   SM(startReStr=(r"\s*(?:" +
                      '|'.join(["%s" % (cName) for cName in qmmmNameList]) +
                      "AMBER)\s*=\s*(?:'|\")?(?P<x_amber_mdin_finline>[\-+0-9.a-zA-Z:]+)(?:'|\")?\s*,?"),
                      coverageIgnore=True,
                      adHoc=lambda p:
                      self.adHoc_read_namelist_stop_parsing(p,
                      stopOnMatchStr=r"(?:(?:^/|\s*/)|&end|&END)\s*$",
                      quitOnMatchStr=None,
                      metaNameStart="x_amber_mdin_",
                      matchNameList=qmmmNameList,
                      matchNameDict=self.qmmmDict,
                      onlyCaseSensitive=True,
                      stopOnFirstLine=True)
                      )
               ])
            ]

    def build_parmKeywordsSimpleMatchers(self):
        parmNameList=getList_MetaStrInDict(self.parmDict)
        return [
            SM(name="parm",
               startReStr=r"\|\s*Version\s*=\s*(?P<x_amber_parm_file_version>[0-9.eEdD]+)" +
                           "\s*Date\s*=\s*(?P<x_amber_parm_file_date>[0-9a-zA-Z\/]+)" +
                           "\s*Time\s*=\s*(?P<x_amber_parm_file_time>[0-9a-zA-Z:]+)",
               endReStr=r"\|\s*Memory\s*Use\s*Allocated\s*",
               subMatchers=[
                   SM(startReStr=(r"\s*(?:" +
                      '|'.join(["%s" % (cName) for cName in parmNameList]) +
                      "AMBER)\s*=\s*(?:'|\")?(?P<x_amber_mdin_finline>[\-+0-9.:a-zA-Z]+)(?:'|\")?\s*,?"),
                      coverageIgnore=True,
                      adHoc=lambda p:
                      self.adHoc_read_namelist_stop_parsing(p,
                      stopOnMatchStr=r"\|\s*Memory\s*Use\s*Allocated\s*",
                      quitOnMatchStr=None,
                      metaNameStart="x_amber_parm_",
                      matchNameList=parmNameList,
                      matchNameDict=self.parmDict,
                      onlyCaseSensitive=True,
                      stopOnFirstLine=True)
                      )
               ])
            ]

    def build_mdoutKeywordsSimpleMatchers(self):
        newDict = self.cntrlDict
        newDict.update(self.qmmmDict)
        mdoutNameList = getList_MetaStrInDict(newDict)
        return [
            SM(name="mdout",
               startReStr=r"\s*General\s*flags:",
               endReStr=r"\s*--{5}--*\s*3\.\s*ATOMIC\s*COORDINATES\s*",
               subMatchers=[
                   SM(startReStr=(r"\s*(?:" +
                      '|'.join(["%s" % (cName) for cName in mdoutNameList]) +
                      "AMBER)\s*=\s*(?:'|\")?(?P<x_amber_mdin_finline>[\-+0-9.:a-zA-Z]+)(?:'|\")?\s*,?"),
                      coverageIgnore=True,
                      adHoc=lambda p:
                      self.adHoc_read_namelist_stop_parsing(p,
                      stopOnMatchStr=r"\s*3\.\s*ATOMIC\s*COORDINATES",
                      quitOnMatchStr=None,
                      metaNameStart="x_amber_mdin_",
                      matchNameList=mdoutNameList,
                      matchNameDict=newDict,
                      onlyCaseSensitive=True,
                      stopOnFirstLine=True)
                      )
               ])
            ]

    def build_fileNameListGroupSubMatcher(self):
        """Builds the Sub Matchers for the Main Parser
        """
        return [SM(r"\|\s*%s:\s*(?P<%s>[0-9a-zA-Z_./\-]+)" %
                (fileNL.matchStr, fileNL.metaHeader + '_' + fileNL.metaNameTag + '_' + fileNL.metaName)) for fileNL in self.fileDict.values()]


    def build_subMatchers(self):
        """Builds the sub matchers to parse the main file of AMBER.
        """
        mdinKeywordsSimpleMatchers = self.build_mdinKeywordsSimpleMatchers()
        parmKeywordsSimpleMatchers = self.build_parmKeywordsSimpleMatchers()
        mdoutKeywordsSimpleMatchers = self.build_mdoutKeywordsSimpleMatchers()
        fileNameListGroupSubMatcher = self.build_fileNameListGroupSubMatcher()

        ########################################
        # submatcher for mdin
        mdinSubMatcher = SM(name='mdinKeywords',
            startReStr=r"\s*Here\s*is\s*the\s*input\s*file:\s*",
            endReStr=r"\s*--{5}--*\s*1\.\s*RESOURCE\s*USE:",
            subFlags=SM.SubFlags.Unordered,
            subMatchers=[
                SM(r"(?:(?P<x_amber_mdin_header>[0-9a-zA-Z]+)?)\s*$",
                    coverageIgnore=True,
                    repeats=False),
                ] + mdinKeywordsSimpleMatchers
            )

        ########################################
        # submatcher for parm
        parmSubMatcher = SM(name='parmKeywords',
            startReStr=r"\s*--{5}--*\s*|\s*Flags:",
            endReStr=r"\s*--{5}--*\s*2\.\s*CONTROL\s*",
            forwardMatch=True,
            #subFlags=SM.SubFlags.Unordered,
            subMatchers=[
                SM(r"\|\s*Flags\s*:\s*(?P<x_amber_parm_flags>[0-9a-zA-Z]+)?\s*$"),
                SM(r"\s*getting\s*new\s*box\s*info\s*from\s*(?:bottom\s*of\s*|netcdf)(?P<x_amber_mdin_finline>[0-9a-zA-Z]+)\s*(file)?\s*",
                   coverageIgnore=True),
                SM(r"\|(\s*NetCDF)?\s*(?P<x_amber_parm_box_info>[_0-9a-zA-Z]+)\s*(?:ntb=[0-9]\s*and\s*igb=[0-9])?:?" +
                    "(?:(?:box|Box)info\s*found|Setting\s*up\s*nonperiodic\s*simulation)\s*"),
                SM(r"\|\s*Largest\s*sphere\s*to\s*fit\s*in\s*unit\s*cell\s*has\s*radius\s*=\s*" +
                    "(?P<x_amber_parm_unitcell_radius>[0-9.eEdD]+)\s*"),
                SM(r"\|\s*(?P<x_amber_parm_file_format>[0-9.eEdD]+)" +
                    "\s*format\s*PARM\s*file\s*being\s*parsed\s*\.")
                ] + parmKeywordsSimpleMatchers + [
                SM(r"\|\s*Real\s*(?P<x_amber_mdin_finline>[0-9.eEdD]+)\s*", coverageIgnore=True),
                SM(r"\|\s*Hollerith\s*(?P<x_amber_mdin_finline>[0-9.eEdD]+)\s*", coverageIgnore=True),
                SM(r"\|\s*Integer\s*(?P<x_amber_mdin_finline>[0-9.eEdD]+)\s*", coverageIgnore=True),
                SM(r"\|\s*Max\s*Pairs\s*(?P<x_amber_mdin_finline>[0-9.eEdD]+)\s*", coverageIgnore=True),
                SM(r"\|\s*nblistReal\s*(?P<x_amber_mdin_finline>[0-9.eEdD]+)\s*", coverageIgnore=True),
                SM(r"\|\s*nblist\s*Int\s*(?P<x_amber_mdin_finline>[0-9.eEdD]+)\s*", coverageIgnore=True),
                SM(r"\|\s*Total\s*(?P<x_amber_parm_total_memory>[0-9.eEdD]+)\s*")
                ]
            )

        ########################################
        # submatcher for mdout
        mdoutSubMatcher = SM(name='mdoutKeywords',
            startReStr=r"\s*General\s*flags:",
            forwardMatch=True,
            #subFlags=SM.SubFlags.Unordered,
            subMatchers=mdoutKeywordsSimpleMatchers
            )

        ########################################
        # submatcher for MD
        mddataNameList=getList_MetaStrInDict(self.mddataDict)
        MDSubMatcher = SM(name='MDStep',
            startReStr=r"\s*(?:NSTEP\s*=|NSTEP\s*ENERGY\s*RMS\s*)",
            #endReStr=r"\s*(?:FINAL\s*RESULTS|A\sV\sE\sR\sA\sG\sE\sS\s*O\sV\sE\sR)",
            forwardMatch=True,
            subMatchers=[
                   SM(startReStr=(r"\s*(?:(?:" +
                   '|'.join(["%s" % (cName) for cName in mddataNameList]) +
                   "AMBER)\s*=\s*|NSTEP\s*ENERGY\s*RMS\s*)(?:'|\")?" +
                   "(?P<x_amber_mdin_finline>[\-+0-9.:a-zA-Z]+)(?:'|\")?\s*,?"),
                   #endReStr=r"\s*(?:FINAL\s*RESULTS|A\sV\sE\sR\sA\sG\sE\sS\s*O\sV\sE\sR)",
                   coverageIgnore=True,
                   adHoc=lambda p:
                   self.adHoc_read_namelist_stop_parsing(p,
                   stopOnMatchStr=r"\s*(?:NSTEP\s*=|NSTEP\s*ENERGY\s*RMS\s*)",
                   quitOnMatchStr=r"\s*(?:FINAL\s*RESULTS|A\sV\sE\sR\sA\sG\sE\sS\s*O\sV\sE\sR)",
                   metaNameStart="x_amber_mdout_",
                   matchNameList=mddataNameList,
                   matchNameDict=self.mddataDict,
                   onlyCaseSensitive=True,
                   stopOnFirstLine=False)
                   )
            ])

        ########################################
        # return main Parser
        return [
            SM(name='NewRun',
                startReStr=r"\s*Amber\s*[0-9]+\s*(?:SANDER|PMEMD)",
                endReStr=r"\s*wallclock\(\) was called\s*[0-9]+\s*times",
                repeats=True,
                required=True,
                forwardMatch=True,
                sections=['section_run'],
                fixedStartValues={'program_name': 'Amber'},
                subMatchers=[
                    # header specifing version, compilation info, task assignment
                    SM(r"\s*Amber\s*(?P<program_version>[0-9a-zA-Z]+)\s*(?P<x_amber_program_module>[0-9a-zA-Z]+)\s*(?P<x_amber_program_version_date>[0-9a-zA-Z]+)"),
                    SM(name='run_date', startReStr=r"\|\s*Run\s*on\s*(?P<x_amber_program_execution_date>[0-9a-zA-Z\/]+)\s*at\s*(?P<x_amber_program_execution_time>[0-9a-zA-Z:]+)"),
                    SM(name='run_path', startReStr=r"\|\s*Executable\s*path\s*:\s*(?P<x_amber_program_execution_path>[0-9a-zA-Z\/]+)"),
                    SM(name='work_path', startReStr=r"\|\s*Working\s*directory\s*:\s*(?P<x_amber_program_working_path>[0-9a-zA-Z\/]+)"),
                    SM(name='run_host', startReStr=r"\|\s*Hostname\s*:\s*(?P<x_amber_program_execution_host>[0-9a-zA-Z\/]+)"),
                    SM(name='FileNameMatch',
                       startReStr=r"\s*File\s*Assignments:",
                       forwardMatch=True,
                       sections=['x_amber_section_input_output_files'],
                       subFlags=SM.SubFlags.Unordered,
                       subMatchers=fileNameListGroupSubMatcher
                       ), # END SectionMethod
                    SM(name='SectionMethod',
                       startReStr=r"\s*Here\s*is\s*the\s*input\s*file:",
                       forwardMatch=True,
                       sections=['section_sampling_method'],
                       subMatchers=[
                           # parse verbatim writeout of mdin file
                           mdinSubMatcher,
                           # parse summary writeout of topology parameters
                           parmSubMatcher,
                           # parse control settings writeout of Amber
                           mdoutSubMatcher
                       ]), # END SectionMethod
                    SM(name='SingleConfigurationCalculationWithSystemDescription',
                       startReStr=r"\s*4\.\s*RESULTS",
                       endReStr=r"\s*(?:FINAL\s*RESULTS|A\sV\sE\sR\sA\sG\sE\sS\s*O\sV\sE\sR)",
                       forwardMatch=True,
                       subMatchers=[
                           # the actual section for a single configuration calculation starts here
                           SM(name='SingleConfigurationCalculation',
                              startReStr=r"\s*(?:NSTEP\s*=|NSTEP\s*ENERGY\s*RMS\s*)",
                              #endReStr=r"\s*(?:FINAL\s*RESULTS|A\sV\sE\sR\sA\sG\sE\sS\s*O\sV\sE\sR)",
                              repeats=True,
                              forwardMatch=True,
                              sections=['section_single_configuration_calculation'],
                              subMatchers=[
                                  # MD
                                  MDSubMatcher
                              ]) # END SingleConfigurationCalculation
                       ]), # END SingleConfigurationCalculationWithSystemDescription
                    SM(r"\s*5\.\s*TIMINGS\s*"),
                    # summary of computation
                    SM(name='ComputationTimings',
                       startReStr=r"\s*Computation timings:",
                       subMatchers=[
                           SM(r"\s*Run\s*done\s*at[0-9:]+"),
                       ]), # END Computation
                    SM(name='end_run', startReStr=r"\s*wallclock\(\)\s*was\s*called\s*[0-9]+\s*times")
                ]) # END NewRun
            ]


class AMBERParserInterface():
   """ A proper class envolop for running this parser from within python. """
   def __init__(self, backend, **kwargs):
       self.backend_factory = backend

   def parse(self, mainfile):
        from unittest.mock import patch
        logging.info('amber parser started')
        logging.getLogger('nomadcore').setLevel(logging.WARNING)
        backend = self.backend_factory("amber.nomadmetainfo.json")
        parserInfo = {'name': 'amber-parser', 'version': '1.0'}
        context = AMBERParser()
        context.coverageIgnore = re.compile(r"^(?:" + r"|".join(context.coverageIgnoreList) + r")$")
        with patch.object(sys, 'argv', ['<exe>', '--uri', 'nmd://uri', mainfile]):
            mainFunction(
                mainFileDescription=context.mainFileDescription(),
                metaInfoEnv=None,
                parserInfo=parserInfo,
                cachingLevelForMetaName=context.cachingLevelForMetaName,
                superContext=context,
                superBackend=backend)
        return backend


if __name__ == "__main__":
    parser = AMBERParser()
    parser.parse()
