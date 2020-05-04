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

from builtins import next
from builtins import object
import numpy as np
import mdtraj
from mdtraj import FormatRegistry
from mdtraj.formats.xtc import load_xtc
from mdtraj.formats.trr import load_trr
from mdtraj.formats.hdf5 import load_hdf5
from mdtraj.formats.lh5 import load_lh5
from mdtraj.formats.netcdf import load_netcdf
from mdtraj.formats.mdcrd import load_mdcrd
from mdtraj.formats.dcd import load_dcd
from mdtraj.formats.binpos import load_binpos
from mdtraj.formats.pdb import load_pdb
from mdtraj.formats.arc import load_arc
from mdtraj.formats.openmmxml import load_xml
from mdtraj.formats.prmtop import load_prmtop
from mdtraj.formats.psf import load_psf
from mdtraj.formats.mol2 import load_mol2
from mdtraj.formats.amberrst import load_restrt, load_ncrestrt
from mdtraj.formats.lammpstrj import load_lammpstrj
from mdtraj.formats.dtr import load_dtr, load_stk
from mdtraj.formats.xyzfile import load_xyz
from mdtraj.formats.hoomdxml import load_hoomdxml
#from mdtraj.formats.tng import load_tng
from mdtraj.core.topology import Topology
from mdtraj.core.trajectory import _TOPOLOGY_EXTS
import mdtraj.formats
import ase.io
import ase.io.formats
from ase.io.formats import format2modulename, import_module
import logging
import os
import sys
logger = logging.getLogger(__name__)

class TrajectoryReader(object):
    """Used to parse various different atomic coordinate files.

    Reading is primarily done by MDTraj and if the file format
    is not supported by MDTraj, ASE is used instead. For file
    formats that are not supported by both codes, a user defined 
    function can be used. See user supported formats in the dictionary.

    Returns all coordinates as numpy arrays.
    """
    def __init__(self):
        self.trajfile = None # File name of the trajectory file with path if needed
        self.topofile = None # File name of the topology file with path if needed
        self.topo_n_atoms = None # Number of atoms at topology file
        self.trajformat = None # Format type of trajectory file
        self.topoformat = None # Format type of topology file
        self.trajhandler = None # The object parsing the trajectory file
        self.trajiter = None # The object parsing the trajectory file
        self.trajchunk = 100 # The object parsing the trajectory file
        self.topohandler = None # The object parsing the topology file
        self.trajcode = None # To explicitly define the parsing library (MDTraj/ASE) for trajectory files
        self.topocode = None # To explicitly define the parsing library (MDTraj/ASE) for topology files
        self.readfirst = False
        self.user_formats = {
            "cif":    {"Crystallographic Information File", self.custom_iread}
        }

    def load(self):
        """Loads the file handles for trajectory and/or topology
        """
        
        if self.topohandler is None:
            if self.topofile:
                self.check_topology_format_support()
        if self.topohandler is None:
            #logger.warning("The topology file format '{}' is not supported by TrajectoryReader.".format(self.topoformat))
            pass

        if self.trajhandler is None:
            self.check_trajectory_format_support()
            if self.trajhandler is None:
                #logger.warning("The trajectory file format '{}' is not supported by TrajectoryReader.".format(self.trajformat))
                #logger.warning("ASE could not read the file '{}' with format '{}'. The contents might be malformed or wrong format used.".format(filename, file_format))
                pass
        return self.trajhandler
    
    def load_topology(self):
        """Loads the file handles for topology only
        """
        
        if self.topohandler is None:
            if self.topofile:
                self.check_topology_format_support()
        if self.topohandler is None:
            logger.error("The topology file format '{}' is not supported by TrajectoryReader.".format(self.topoformat))

        return self.topohandler
    
    def get_file_format(self, filename, givenformat):
        """Returns extension of file
        """
        file_format = None
        if givenformat:
            fname = filename.split("/")[-1]
            # MDTraj expects that extensions start with dot
            if '.' not in fname:
                file_format = "." + fname.lower()
            else:
                file_format = fname.lower()
        else:
            extension = filename.split(".")[-1]
            file_format = "." + extension.lower()
        return file_format

    def check_topology_format_support(self):
        """Check if the given format is supported.
        """
        topofilename = os.path.basename(self.topofile)
        if self.topoformat is None:
            file_format = self.get_file_format(topofilename, self.topoformat)
            self.topoformat = file_format
        else:
            file_format = self.topoformat

        if (file_format in _TOPOLOGY_EXTS and 
            self.topocode is not "ase"):
            self.topohandler = self.load_mdtraj_topology(file_format)
            self.topocode = "mdtraj"

        # If MDTraj does not have support for the format 
        # or can not load the topology, use ASE instead
        ase_support = None
        if self.topohandler is None:
            logger.warning("MDTraj could not read the file '{}' with format '{}'. The contents might be malformed or wrong format used.".format(self.topofile, file_format))
            ase_support = self.get_ase_format_support(file_format)
            # May still have chance that ASE can recognize the 
            # format with its filetype checking function
            if ase_support is None:
                ase_support = ase.io.formats.filetype(topofilename)
                self.topohandler = self.load_ase_support(topofilename, file_format=ase_support)
                self.topocode = "ase"
            else:
                self.topohandler = self.load_ase_support(topofilename, file_format=file_format)
                self.topocode = "ase"

        # May add openbabel and/or Pybel support here.
        if self.topohandler is None:
            logger.error("The format '{}' is not supported by CoordinateReader.".format(file_format))
            return False
        else:
            return True

    def check_trajectory_format_support(self):
        """Check if the given format is supported.
        """
        mdtraj_failed = False
        trajfilename = os.path.basename(self.trajfile)
        if self.trajformat is None:
            file_format = self.get_file_format(trajfilename, self.trajformat)
            self.trajformat = file_format
        else:
            file_format = self.trajformat

        # First check whether MDtraj has support for the file type
        # trajhandler_check = FormatRegistry.loaders[file_format]
        trajhandler_check = None
        try:
            trajhandler_check = FormatRegistry.fileobjects[file_format]
        except KeyError:
            pass
        else:
            if trajhandler_check:
                self.trajhandler = self.mdtraj_iread(mdtraj_handler=trajhandler_check)
                self.trajcode = "mdtraj"

        # If MDTraj does not have support for the format 
        # or can not load the trajectory, use ASE instead
        ase_support = None
        if self.trajhandler is None:
            logger.warning("MDTraj could not read the file '{}' with format '{}'. The contents might be malformed or wrong format used.".format(self.trajfile, file_format))
            ase_support = self.get_ase_format_support(file_format)
            # May still have chance that ASE can recognize the 
            # format with its filetype checking function
            if ase_support is None:
                ase_support = ase.io.formats.filetype(trajfilename)
                self.trajhandler = self.load_ase_support(self.trajfile, file_format=ase_support)
                self.trajcode = "ase"
            else:
                self.trajhandler = self.load_ase_support(self.trajfile, file_format=file_format)
                self.trajcode = "ase"

        # May add openbabel and/or Pybel support here.

        if self.trajhandler is None:
            logger.error("The format '{}' is not supported by CoordinateReader.".format(file_format))
            return False
        else:
            return True

    def iread(self):
        """Returns an iterator that goes through the given trajectory file one
        configuration at a time.
        """
        if self.trajhandler is not None:
            iterator_object = iter(self.trajhandler)
            try:
                while True:
                    self.trajiter = next(iterator_object)
                    return self.trajiter
            except StopIteration:
                pass
            finally:
                del iterator_object
        else:
            return None

    def get_topology(self):
        """Returns an iterator that goes through the given trajectory file one
        configuration at a time.
        """
        return self.topohandler

    def get_unitcell(self):
        """Returns an iterator that goes through the given trajectory file one
        configuration at a time.
        """
        return self.trajhandler.cell_lengths

    def load_mdtraj_topology(self, file_format, **kwargs):
        """Borrowed from mdtraj.trajectory to remove 
           IOError and TypeError warnings
    
        Get the topology from a argument of indeterminate type
        If top is a string, we try loading a pdb, if its a trajectory
        we extract its topology.

        Returns
        -------
        topology : md.Topology
        """

        topology = None
        wrapkwargs = kwargs.copy()
        wrapkwargs.pop("top", None)
        wrapkwargs.pop("atom_indices", None)
        wrapkwargs.pop("frame", None)

        ext=file_format
        top=self.topofile

        if ext in ['.pdb', '.pdb.gz', '.h5','.lh5']:
            _traj = mdtraj.load_frame(top, 0, **wrapkwargs)
            topology = _traj.topology
        elif ext in ['.prmtop', '.parm7']:
            topology = load_prmtop(top, **wrapkwargs)
        elif ext in ['.psf']:
            topology = load_psf(top, **wrapkwargs)
        elif ext in ['.mol2']:
            topology = load_mol2(top, **wrapkwargs).topology
        elif ext in ['.gro']:
            topology = mdtraj.core.trajectory.load_gro(top, **wrapkwargs).topology
        elif ext in ['.arc']:
            topology = load_arc(top, **wrapkwargs).topology
        elif ext in ['.hoomdxml']:
            topology = load_hoomdxml(top, **wrapkwargs).topology
        elif isinstance(top, Trajectory):
            topology = top.topology
        elif isinstance(top, Topology):
            topology = top

        return topology

    def load_ase_support(self, file_name, file_format):
        """Get the topology from a argument of indeterminate type
        If top is a string, we try loading a pdb, if its a trajectory
        we extract its topology.

        Returns
        -------
        topology : generator
        """
        return ase.io.read(file_name, index=":", format=file_format)

    def get_ase_format_support(self, file_format):
        """ Inspired from ase.io.formats.get_ioformat

        Get support level of ASE for file formats

        Returns
        -------
        IOformat list
        """

        ioformat = None
        _format = file_format.replace('-', '_')
        module_name = format2modulename.get(file_format, _format)
        module_handler = None
        module_handler = import_module('ase.io.' + module_name)
        if module_handler:
            module_read = getattr(module_handler, 'read_' + _format, None)
            if module_read and not inspect.isgeneratorfunction(module_read):
                read = functools.partial(wrap_read_function, module_read)
                if module_read:
                    ioformat = all_formats[file_format]

        return ioformat

    def natoms(self):
        """Read the first configuration of the coordinate file to extract the
        number of atoms in it.
        """
        return self.trajiter.xyz.shape[1]

    def ase_iread(self):
        """Wrapper for ASE's iread function. Returns numpy arrays instead of
        Atoms objects.
        """
        # The newest ASE version found in Github has an iread function.
        # After reading the ASE source code, it seems that the ASE iread does
        # actually read the entire file into memory and the yields the
        # configurations from it. Should be checked at some point.
        def ase_generator(iterator):
            """Used to wrap an iterator returned by ase.io.iread so that it returns
            the positions instead of the ase.Atoms object.
            """
            for value in iterator:
                yield value.get_positions()

        iterator = self.trajhandler
        return ase_generator(iterator)

    def mdtraj_iread(self, mdtraj_handler=None, **kwargs):
        """Generator function that is used to read an atomic configuration file (MD
        trajectory, geometry optimization, static snapshot) from a file one frame
        at a time. Only the xyz positions are returned from the file, and no unit
        conversion is done, so you have to be careful with units.

        By using a generator pattern we can avoid loading the entire trajectory
        file into memory. This function will instead load a chunk of the file into
        memory (with MDTraj you can decide the chunk size, with ASE it seems to
        always be one frame), and serve individual files from that chunk. Once the
        frames in one chunk are iterated, the chunk will be garbage collected and
        memory is freed.

        Args:
            filename: String for the file path.
            file_format: String for the file format. If not given the format is
                automatically detected from the extension.

        Yields:
            numpy array containing the atomic positions in one frame.
        """

        # Try to open the file with MDTraj first. With a brief inspection it seems
        # that MDTraj is better performance wise, because it can iteratively load a
        # "chunk" of frames, and still serve the individual frames one by one. ASE
        # on the other hand will iteratively read frames one by one (unnecessary
        # IO).

        # Must use the low level MDTraj API to open files without topology.
        # mdtraj_supported_format = FormatRegistry.loaders[self.fileformat]
        if mdtraj_handler is not None:
            if self.trajchunk == 0:
                try:
                    loader = FormatRegistry.loaders[self.trajformat]
                except KeyError:
                    return
                # If chunk was 0 then we want to avoid filetype-specific code
                # in case of undefined behavior in various file parsers.
                # TODO: this will first apply stride, then skip!
                if self.trajformat not in _TOPOLOGY_EXTS:
                    topkwargs = kwargs.copy()
                    topkwargs['top'] = self.topohandler
                    # standard_names is a valid keyword argument only for files containing topologies
                    topkwargs.pop('standard_names', None)
                    positions = loader(self.trajfile, **topkwargs)
                else:
                    positions = loader(self.trajfile)
                for pos in positions:
                    yield pos
            elif self.trajformat in ('.pdb', '.pdb.gz'):
                # the PDBTrajectortFile class doesn't follow the standard API. Fixing it here
                try:
                    loader = FormatRegistry.loaders[self.trajformat]
                except KeyError:
                    return
                t = loader(filename)
                for i in range(0, len(t), self.trajchunk):
                    positions = t[i:i+self.trajchunk]
                    for pos in positions:
                        yield pos

            else:
                if self.topohandler is None:
                    n_atoms_set=None
                    if self.trajformat in ('.crd', '.mdcrd'):
                        return
                else:
                    n_atoms_set=self.topohandler.n_atoms
                try:
                    with (lambda x: mdtraj_handler(x, n_atoms=n_atoms_set)
                          if self.trajformat in ('.crd', '.mdcrd')
                          else mdtraj_handler(self.trajfile, mode="r"))(self.trajfile) as f:
                        empty = False
                        while not empty:
                            if self.trajformat not in _TOPOLOGY_EXTS:
                                data = f.read_as_traj(self.topohandler, n_frames=self.trajchunk)
                            else:
                                data = f.read_as_traj(n_frames=self.trajchunk)
                            if isinstance(data, tuple):
                                positions = data[0]
                            else:
                                positions = data
                            if len(positions) == 0:
                                empty = True
                            else:
                                for pos in positions:
                                    yield pos
                except IOError:
                    #logger.warning("MDTraj could not read the file '{}' with format '{}'. The contents might be malformed or wrong format used.".format(self.trajfile, self.trajformat))
                    return

    def custom_iread(self, file_handle, format):
        """
        """
        pass

if __name__ == "__main__":
    topo_file = sys.argv[1]
    traj_file = sys.argv[2]
    readdata = TrajectoryReader()
    readdata.trajfile = traj_file
    readdata.topofile = topo_file
#    readdata.trajformat = '.netcdf'
    traj_iterator = readdata.load()
    mytopology = readdata.get_topology()
    # print([ atom.name for atom in mytopology.atoms ])
    readdata.trajchunk=300
    readaccess = False
    # check whether mdcrd can be read with default format
    atom = readdata.iread()
    i=0
    if atom is None:
        # mdcrd trajectory file might be in NetCDF format
        # Trying to read with netcdf format
        readdata = TrajectoryReader()
        readdata.trajfile = traj_file
        readdata.topofile = topo_file
        readdata.trajformat = '.netcdf'
        readdata.trajchunk=300
        traj_iterator = readdata.load()
        mytopology = readdata.get_topology()
        atom = readdata.iread()
    
    if atom is not None:
        readaccess = True
    
#    if readaccess:
#        while True:
#            if atom is not None:
#                print(atom)
#                print(readdata.natoms())
#                i += 1
#            else:
#                break
#            atom = readdata.iread()
    print(i)
#    angles = mdtraj.compute_angles(mytopology)
#    dihedrals = mdtraj.compute_dihedrals(mytopology)
#    print(dihedrals)
    #print(' Chains: %s' % [res for res in mytopology.chains])
    #print(' Residues: %s' % [res for res in mytopology.residues])
    #print(' Atoms: %s' % [res for res in mytopology.atoms])
    #print(' Bonds: %s' % [res for res in mytopology.bonds])
    #print(mytopology.chain(0))
    table, bonds = mytopology.to_dataframe()
    print(bonds)
    print(len(bonds))
    array = table.get("resSeq")
    #print([res for res in array])
    mydict = table.to_dict(orient='list')
    mylist = mydict["resSeq"]
    print('mylist')
    print(mylist)
    print('mylist')
    atomindex = np.arange(len(mylist))
    passarray = np.zeros((len(mylist), 2), dtype=int)
    passarray[:,0] = atomindex
    passarray[:,1] = mylist
    print(table)
    print(passarray)
    #print(bonds)
#    for bond in bonds:
#        print(bond)
    mol_to_mol_bond = []
    mol_to_mol_bond_num = []
    atom_in_mol_bond = []
    atom_in_mol_bond_num = []
    index=0
    print(len([res for res in mytopology.bonds]))
    for res in mytopology.bonds:
        molname1, molatom1 = str(res[0]).split('-')
        molname2, molatom2 = str(res[1]).split('-')
        if molname1 in molname2: 
            atom_in_mol_bond.append(res)
            atom_in_mol_bond_num.append(bonds[index])
        else:
            mol_to_mol_bond.append(res)
            mol_to_mol_bond_num.append(bonds[index])
        index += 1
    print(mol_to_mol_bond)
    print(np.array(mol_to_mol_bond_num))

    atom_list = list(set(mydict["name"]))
    atom_type_dict = {}
    print(atom_list)
    print(len(atom_list))
    masses = {}
    for ielm in range(len(atom_list)):
        elm = atom_list[ielm]
        atom_type_dict.update({elm : ielm})
        for atom in mytopology.atoms:
            if elm == atom.name:
                masses.update({atom.name : atom.element})
    print([atom_list, masses])
    interDict = {}
    interTypeDict = {}
    for elm in atom_list:
        bondList = []
        typeList = []
        bondid = 0
        for bond in mytopology.bonds:
            molname1, molatom1 = str(bond[0]).split('-')
            molname2, molatom2 = str(bond[1]).split('-')
            if (elm == str(molatom1) or elm == str(molatom2)):
                bondList.append(list(bonds[bondid]))
                #typeList.append(list([molatom1,molatom2]))
                typeList.append(list([
                    atom_type_dict[str(molatom1)],
                    atom_type_dict[str(molatom2)]
                    ]))
                interDict.update({elm : bondList})
                interTypeDict.update({elm : typeList})
            bondid += 1

#    print(interDict)
    print(interDict["CH3"])
    print(interTypeDict["CH3"])
    print(len(interDict["CH3"][0]))

#    print([[atom.element.radius, atom.index, atom.name, atom.residue.name, atom.element] for atom in mytopology.atoms])

# Atom labels:
    print(table)
#    for atom in mytopology.atoms:
#        print(atom)
    labellist = mydict["element"]
    #print(labellist)
#    atomlabels = np.arange(len(labellist))
#    labelarray = np.zeros((len(labellist), 2), dtype=int)
#    labelarray[:,0] = atomlabels
#    labelarray[:,1] = labellist
#    print(labelarray)
#    print(readdata.unitcell_lengths)
    latticevectors = readdata.iread().unitcell_vectors
    box = readdata.iread().unitcell_lengths
    step = 0
    for latvec in latticevectors:
        print(step, latvec)
        step += 1
    for pbc in box:
        print(step, box)
        step += 1



            





