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

import numpy as np
import logging
import json
import os
import re
from collections import namedtuple

class MetaInfoMap(dict):
    """Map cache values to meta info
    """
    activeInfo=False
    infoPurpose=None
    defaultValue=None
    nameTranslate=None
    matchStr=None
    metaHeader=None
    metaName=None
    metaNameTag=None
    metaInfoType=None
    value=None
    valueSize=None
    sizeMetaName=None
    depends=None
    lookupdict=None
    subfunction=None
    activeSections=None
    autoSections=False

    def __init__(self, *args, **kwargs):
        super(MetaInfoMap, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    if k in self:
                        self[k] = v
        if kwargs:
            for k, v in kwargs.items():
                if k in self:
                    self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(MetaInfoMap, self).__setitem__(key, value)
        self.__dict__.update({key: value})

class FileInfoMap(dict):
    """Map cache values to meta info
    """
    activeInfo=False
    infoPurpose=None
    fileName=None
    fileFormat=None
    fileSupplied=False
    fileHolder=None
    nameTranslate=None
    matchStr=None
    metaHeader=None
    metaName=None
    metaNameTag=None
    metaInfoType=None
    value=None
    valueSize=None
    sizeMetaName=None
    depends=None
    lookupdict=None
    subfunction=None
    activeSections=None

    def __init__(self, *args, **kwargs):
        super(FileInfoMap, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    if k in self:
                        self[k] = v
        if kwargs:
            for k, v in kwargs.items():
                if k in self:
                    self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(FileInfoMap, self).__setitem__(key, value)
        self.__dict__.update({key: value})

class MapDictionary(dict):
    """
    Modified from the reference source below:
    https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
    Example:
    m = MapDictionary({'Name': 'mdtraj'}, format='.mdcrd', found=True, list=['Value'])
    """
    def __init__(self, *args, **kwargs):
        super(MapDictionary, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    if (isinstance(v, FileInfoMap) or
                        isinstance(v, MetaInfoMap)):
                        if v.nameTranslate:
                            v.metaName = v.nameTranslate(k)
                        else:
                            v.metaName = k
                    v.matchStr = k
                    metaStr = ''
                    if v.metaHeader:
                        metaStr = metaStr + v.metaHeader + '_'
                    if v.metaNameTag:
                        metaStr = metaStr + v.metaNameTag + '_'
                    metaStr = metaStr + v.metaName
                    self[metaStr] = v
                    if metaStr != k:
                        self.pop(k, None)

        if kwargs:
            for k, v in kwargs.items():
                if (isinstance(v, FileInfoMap) or
                    isinstance(v, MetaInfoMap)):
                    if v.metaTranslate:
                        v.metaName = v.nameTranslate(k)
                    else:
                        v.metaName = k
                v.matchStr = k
                metaStr = ''
                if v.metaHeader:
                    metaStr = metaStr + v.metaHeader + '_'
                if v.metaNameTag:
                    metaStr = metaStr + v.metaNameTag + '_'
                metaStr = metaStr + v.metaName
                self[metaStr] = v
                if metaStr != k:
                    self.pop(k, None)

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(MapDictionary, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(MapDictionary, self).__delitem__(key)
        del self.__dict__[key]

    def get_keys(self):
        return [val.metaName for val in self.__dict__.values()]

def get_unitDict(keyname):
    """ Unit dictionary for convertions

        Unit names will be converted to values.
        When defining units in translator dictionary,
        the unit names in the dictionary should be used.
        The unit convertion values are written for SI units.
        If you would like to change it, just add another key
        to the dictionary and change the key at parser base.
        Usage:
           Natively support python language in definitions.
           You can use any python math operator and function.
           Moreover, you can use space or - for multiplication
           and ^ for power of values.
        Example:
            kilogram/meter^2 can be written as
            kilo-gram/meter^2 or kilo-gram/meter**2
            and will be calculated as
            kilo*gram/meter**2

        For a quick check for AMBER units please see:
        http://ambermd.org/Questions/units.html
    """
    unitDict = {
        "si" : {
            "meter"          : "1.0",
            "kilo"           : "1.0e3",
            "gram"           : "1.0e-3",
            "second"         : "1.0",
            "joule"          : "1.0",
            "newton"         : "1.0",
            "kelvin"         : "1.0",
            "pascal"         : "1.0",
            "coulomb"        : "1.0",
            "volt"           : "1.0",
            "centi"          : "1.0e-2",
            "milli"          : "1.0e-3",
            "micro"          : "1.0e-6",
            "nano"           : "1.0e-9",
            "pico"           : "1.0e-12",
            "femto"          : "1.0e-15",
            "atto"           : "1.0e-18",
            "erg"            : "1.0e-7",
            "dyne"           : "1.0e-5",
            "bar"            : "1.0e-1",
            "angstrom"       : "1.0e-10",
            "kcal"           : "4184.096739614824",
            "mol"           : "0.602213737699784e24",
            "atmosphere"     : "1.01325e5",
            "electron"       : "1.602176565e-19",
            "atomicmassunit" : "1.66054e-27",
            "amu"            : "1.66054e-27",
            "bohr"           : "5.29177249e-11",
            "hartree"        : "4.35974e-18",
            "pascal"         : "1.0",
            }}
    if keyname == "amber":
        resDict = unitDict["si"]
    else:
        resDict = unitDict["si"]
    return resDict

def metaNameConverter(keyName):
    newName = keyName.lower().replace(" ", "").replace("-", "")
    newName = newName.replace("(", "").replace(")", "")
    newName = newName.replace("[", "").replace("]", "")
    newName = newName.replace(",", "").replace(".", "")
    newName = newName.replace("\\", "").replace("/", "")
    return newName

def get_fileListDict():
    """Loads dictionary for file namelist of AMBER.

    Returns:
        the list of defaults file namelists
    """
    # Default topology format of Amber is parm, prmtop file
    # As of Amber 9, default trajectory format of Amber is in binary NetCDF format.
    # The alternative is formatted ASCII (mdcrd) format and the format will be
    #     determined after parsing the input control parameters
    # Default input coordinate file format is auto-detected at run time by Amber.
    # The file format can be either formatted ASCII (inpcrd) or NetCDF.
    # The format will be determined by checking the file with load
    #     and iread functions that are supplied by TrajectoryReader.
    # Default restart file is also in NetCDF format.
    # Optionally, velocities and forces can be written to
    # mdvel (.mdvel) and mdfrc (.mdfrc) files, respectively
    startpage = {
        'nameTranslate'   :  metaNameConverter,
        'metaHeader'      : 'x_amber',
        'metaNameTag'     : 'mdin_file',
        'metaInfoType'    : 'C',
        'activeMetaNames' : [],
        'activeSections'  : ['x_amber_section_input_output_files']
        }
    namelist = {
        'MDIN'   : FileInfoMap(startpage),
        'MDOUT'  : FileInfoMap(startpage),
        'INPCRD' : FileInfoMap(startpage, activeInfo=True, fileFormat=['.inpcrd', '.ncrst'],
                               infoPurpose=['inputcoordinates', 'inputvelocities', 'inputunitcell']),
        'PARM'   : FileInfoMap(startpage, activeInfo=True, fileFormat=['.prmtop', '.parm7'],
                               infoPurpose=['topology', 'unitcell']),
        'RESTRT' : FileInfoMap(startpage),
        'REFC'   : FileInfoMap(startpage),
        'MDVEL'  : FileInfoMap(startpage, activeInfo=True, infoPurpose=['velocities']),
        'MDFRC'  : FileInfoMap(startpage, activeInfo=True, infoPurpose=['forces']),
        'MDEN'   : FileInfoMap(startpage),
        'MDCRD'  : FileInfoMap(startpage, activeInfo=True, fileFormat=['.netcdf', '.mdcrd'],
                               infoPurpose=['trajectory', 'unitcell']),
        'MDINFO' : FileInfoMap(startpage),
        'MTMD'   : FileInfoMap(startpage),
        'INPDIP' : FileInfoMap(startpage),
        'RSTDIP' : FileInfoMap(startpage),
        'INPTRA' : FileInfoMap(startpage)
        }
    return MapDictionary(namelist)

def get_nameListDict(deflist):
    """Loads namelist data of AMBER.

    Args:
        deflist: name list definition (cntrl/ewald/qmmm/wt).

    Returns:
        the list of namelists
    """
    startpage = {
        'nameTranslate'   :  metaNameConverter,
        'metaHeader'      : 'x_amber',
        'metaNameTag'     : 'mdin',
        'metaInfoType'    : 'C',
        'activeMetaNames' : [],
        'activeSections'  : ['x_amber_mdin_method']
        }
    cntrllist = {
        'imin' : MetaInfoMap(startpage, defaultValue=0),
        'nmropt' : MetaInfoMap(startpage, defaultValue=0),
        'ntx' : MetaInfoMap(startpage, defaultValue=1),
        'irest' : MetaInfoMap(startpage, defaultValue=0),
        'ntxo' : MetaInfoMap(startpage, defaultValue=2),
        'ntpr' : MetaInfoMap(startpage, defaultValue=50),
        'ntave' : MetaInfoMap(startpage, defaultValue=0),
        'ntwr' : MetaInfoMap(startpage, defaultValue='nstlim'),
        'iwrap' : MetaInfoMap(startpage, defaultValue=0),
        'ntwx' : MetaInfoMap(startpage, defaultValue=0),
        'ntwv' : MetaInfoMap(startpage, defaultValue=0),
        'ntwf' : MetaInfoMap(startpage, defaultValue=0),
        'ntwe' : MetaInfoMap(startpage, defaultValue=0),
        'ioutfm' : MetaInfoMap(startpage, defaultValue=1),
        'ntwprt' : MetaInfoMap(startpage, defaultValue=0),
        'idecomp' : MetaInfoMap(startpage, defaultValue=0),
        'ibelly' : MetaInfoMap(startpage, defaultValue=0),
        'ntr' : MetaInfoMap(startpage, defaultValue=0),
        'restraint_wt' : MetaInfoMap(startpage),
        'restraintmask' : MetaInfoMap(startpage),
        'bellymask' : MetaInfoMap(startpage),
        'maxcyc' : MetaInfoMap(startpage, defaultValue=1),
        'ncyc' : MetaInfoMap(startpage, defaultValue=10),
        'ntmin' : MetaInfoMap(startpage, defaultValue=1),
        'dx0' : MetaInfoMap(startpage, defaultValue=0.01),
        'drms' : MetaInfoMap(startpage, defaultValue=1.0E-4),
        'nstlim' : MetaInfoMap(startpage, defaultValue=1),
        'nscm' : MetaInfoMap(startpage, defaultValue=1000),
        't' : MetaInfoMap(startpage, defaultValue=0.0),
        'dt' : MetaInfoMap(startpage, defaultValue=0.001),
        'nrespa' : MetaInfoMap(startpage, defaultValue=0),
        'ntt' : MetaInfoMap(startpage, defaultValue=0),
        'temp0' : MetaInfoMap(startpage, defaultValue=300),
        'temp0les' : MetaInfoMap(startpage, defaultValue=-1),
        'tempi' : MetaInfoMap(startpage, defaultValue=0.0),
        'ig' : MetaInfoMap(startpage, defaultValue=-1),
        'tautp' : MetaInfoMap(startpage, defaultValue=1.0),
        'gamma_ln' : MetaInfoMap(startpage, defaultValue=0),
        'vrand' : MetaInfoMap(startpage, defaultValue=1000),
        'vlimit' : MetaInfoMap(startpage, defaultValue=0),
        'nkija' : MetaInfoMap(startpage, defaultValue=1),
        'idistr' : MetaInfoMap(startpage, defaultValue=0),
        'sinrtau' : MetaInfoMap(startpage, defaultValue=1.0),
        'ntp' : MetaInfoMap(startpage, defaultValue=0),
        'barostat' : MetaInfoMap(startpage, defaultValue=1),
        'mcbarint' : MetaInfoMap(startpage, defaultValue=100),
        'pres0' : MetaInfoMap(startpage, defaultValue=1.0),
        'comp' : MetaInfoMap(startpage, defaultValue=44.6),
        'taup' : MetaInfoMap(startpage, defaultValue=1.0),
        'csurften' : MetaInfoMap(startpage, defaultValue=0),
        'gamma_ten' : MetaInfoMap(startpage, defaultValue=0),
        'ninterface' : MetaInfoMap(startpage, defaultValue=2),
        'ntc' : MetaInfoMap(startpage, defaultValue=0),
        'tol' : MetaInfoMap(startpage, defaultValue=0.00001),
        'jfastw' : MetaInfoMap(startpage, defaultValue=0),
        'noshakemask' : MetaInfoMap(startpage, defaultValue=''),
        'ivcap' : MetaInfoMap(startpage, defaultValue=0),
        'fcap' : MetaInfoMap(startpage),
        'cutcap' : MetaInfoMap(startpage),
        'xcap' : MetaInfoMap(startpage),
        'ycap' : MetaInfoMap(startpage),
        'zcap' : MetaInfoMap(startpage),
        'iscale' : MetaInfoMap(startpage, defaultValue=0),
        'noeskp' : MetaInfoMap(startpage, defaultValue=1),
        'ipnlty' : MetaInfoMap(startpage, defaultValue=1),
        'mxsub' : MetaInfoMap(startpage, defaultValue=1),
        'scalm' : MetaInfoMap(startpage, defaultValue=100),
        'pencut' : MetaInfoMap(startpage, defaultValue=0.1),
        'tausw' : MetaInfoMap(startpage, defaultValue=0.1),
        'iemap' : MetaInfoMap(startpage, defaultValue=0),
        'gammamap' : MetaInfoMap(startpage, defaultValue=1),
        'ntf' : MetaInfoMap(startpage, defaultValue=1),
        'ntb' : MetaInfoMap(startpage, defaultValue=1),
        'dielc' : MetaInfoMap(startpage, defaultValue=1.0),
        'cut' : MetaInfoMap(startpage, defaultValue=8.0),
        'fswitch' : MetaInfoMap(startpage, defaultValue=-1),
        'nsnb' : MetaInfoMap(startpage, defaultValue=25),
        'ipol' : MetaInfoMap(startpage, defaultValue=0),
        'ifqnt' : MetaInfoMap(startpage, defaultValue=0),
        'igb' : MetaInfoMap(startpage, defaultValue=0),
        'irism' : MetaInfoMap(startpage, defaultValue=0),
        'ievb' : MetaInfoMap(startpage, defaultValue=0),
        'iamoeba' : MetaInfoMap(startpage, defaultValue=0),
        'lj1264' : MetaInfoMap(startpage),
        'efx' : MetaInfoMap(startpage, defaultValue=0),
        'efy' : MetaInfoMap(startpage, defaultValue=0),
        'efz' : MetaInfoMap(startpage, defaultValue=0),
        'efn' : MetaInfoMap(startpage, defaultValue=0),
        'efphase' : MetaInfoMap(startpage),
        'effreq' : MetaInfoMap(startpage)
        }

    wtlist = {
        'istep1' : MetaInfoMap(startpage, defaultValue=0),
        'istep2' : MetaInfoMap(startpage, defaultValue=0),
        'value1' : MetaInfoMap(startpage),
        'value2' : MetaInfoMap(startpage),
        'iinc' : MetaInfoMap(startpage),
        'imult' : MetaInfoMap(startpage)
        }

    ewaldlist = {
        'nfft1' : MetaInfoMap(startpage),
        'nfft2' : MetaInfoMap(startpage),
        'nfft3' : MetaInfoMap(startpage),
        'order' : MetaInfoMap(startpage, defaultValue=4),
        'verbose' : MetaInfoMap(startpage, defaultValue=0),
        'ew_type' : MetaInfoMap(startpage, defaultValue=0),
        'dsum_tol' : MetaInfoMap(startpage, defaultValue=1.0E-5),
        'rsum_tol' : MetaInfoMap(startpage, defaultValue=5.0E-5),
        'mlimit\(1\)' : MetaInfoMap(startpage),
        'mlimit\(2\)' : MetaInfoMap(startpage),
        'mlimit\(3\)' : MetaInfoMap(startpage),
        'ew_coeff' : MetaInfoMap(startpage),
        'nbflag' : MetaInfoMap(startpage, defaultValue=1),
        'skinnb' : MetaInfoMap(startpage, defaultValue=2.0),
        'nbtell' : MetaInfoMap(startpage, defaultValue=0),
        'netfrc' : MetaInfoMap(startpage, defaultValue=1),
        'vdwmeth' : MetaInfoMap(startpage, defaultValue=1),
        'eddmet' : MetaInfoMap(startpage, defaultValue=1),
        'eedtbdns' : MetaInfoMap(startpage, defaultValue=500),
        'column_fft' : MetaInfoMap(startpage, defaultValue=0),
        'ips' : MetaInfoMap(startpage, defaultValue=0),
        'raips' : MetaInfoMap(startpage, defaultValue=-1),
        'mipsx' : MetaInfoMap(startpage, defaultValue=-1),
        'mipsy' : MetaInfoMap(startpage, defaultValue=-1),
        'mipsz' : MetaInfoMap(startpage, defaultValue=-1),
        'mipso' : MetaInfoMap(startpage, defaultValue=4),
        'gridips' : MetaInfoMap(startpage, defaultValue=2),
        'dvbips' : MetaInfoMap(startpage, defaultValue=1.0E-8),
        'frameon' : MetaInfoMap(startpage, defaultValue=1),
        'chngmask' : MetaInfoMap(startpage, defaultValue=1),
        'indmeth' : MetaInfoMap(startpage),
        'diptol' : MetaInfoMap(startpage, defaultValue=0.0001),
        'maxiter' : MetaInfoMap(startpage, defaultValue=20),
        'dipmass' : MetaInfoMap(startpage, defaultValue=0.33),
        'diptau' : MetaInfoMap(startpage, defaultValue=11),
        'irstdip' : MetaInfoMap(startpage),
        'scaldip' : MetaInfoMap(startpage, defaultValue=1)
        }

    qmmmlist = {
        'qm_theory' : MetaInfoMap(startpage),
        'dftb_slko_path' : MetaInfoMap(startpage),
        'dftb_disper' : MetaInfoMap(startpage),
        'dftb_3rd_order' : MetaInfoMap(startpage),
        'dftb_chg' : MetaInfoMap(startpage),
        'dftb_telec' : MetaInfoMap(startpage),
        'dftb_maxiter' : MetaInfoMap(startpage),
        'qmcharge' : MetaInfoMap(startpage),
        'spin' : MetaInfoMap(startpage),
        'qmqmdx' : MetaInfoMap(startpage),
        'verbosity' : MetaInfoMap(startpage),
        'tight_p_conv' : MetaInfoMap(startpage),
        'scfconv' : MetaInfoMap(startpage),
        'pseudo_diag' : MetaInfoMap(startpage),
        'diag_routine' : MetaInfoMap(startpage),
        'printcharges' : MetaInfoMap(startpage),
        'print_eigenvalues' : MetaInfoMap(startpage),
        'qxd' : MetaInfoMap(startpage),
        'parameter_file' : MetaInfoMap(startpage),
        'peptide_corr' : MetaInfoMap(startpage),
        'itrmax' : MetaInfoMap(startpage),
        'ntpr' : MetaInfoMap(startpage),
        'grms_tol' : MetaInfoMap(startpage),
        'ndiis_attempts' : MetaInfoMap(startpage),
        'ndiis_matrices' : MetaInfoMap(startpage),
        'vshift' : MetaInfoMap(startpage),
        'errconv' : MetaInfoMap(startpage),
        'qmmm_int' : MetaInfoMap(startpage),
        'qmmask' : MetaInfoMap(startpage),
        'qmcut' : MetaInfoMap(startpage),
        'lnk_dis' : MetaInfoMap(startpage),
        'dftb_telec_step' : MetaInfoMap(startpage),
        'fockp_d1' : MetaInfoMap(startpage),
        'fockp_d2' : MetaInfoMap(startpage),
        'fockp_d3' : MetaInfoMap(startpage),
        'fockp_d4' : MetaInfoMap(startpage),
        'pseudo_diag_criteria' : MetaInfoMap(startpage),
        'damp' : MetaInfoMap(startpage),
        'kappa' : MetaInfoMap(startpage),
        'min_heavy_mass' : MetaInfoMap(startpage),
        'r_switch_hi' : MetaInfoMap(startpage),
        'r_switch_lo' : MetaInfoMap(startpage),
        'iqmatoms' : MetaInfoMap(startpage),
        'qmgb' : MetaInfoMap(startpage),
        'lnk_atomic_no' : MetaInfoMap(startpage),
        'lnk_method' : MetaInfoMap(startpage),
        'printbondorders' : MetaInfoMap(startpage),
        'buffercharge' : MetaInfoMap(startpage),
        'printdipole' : MetaInfoMap(startpage),
        'qmshake' : MetaInfoMap(startpage),
        'qmmmrij_incore' : MetaInfoMap(startpage),
        'qmqm_erep_incore' : MetaInfoMap(startpage),
        'qm_ewald' : MetaInfoMap(startpage),
        'qm_pme' : MetaInfoMap(startpage),
        'kmaxqx' : MetaInfoMap(startpage),
        'kmaxqy' : MetaInfoMap(startpage),
        'kmaxqz' : MetaInfoMap(startpage),
        'ksqmaxsq' : MetaInfoMap(startpage),
        'adjust_q' : MetaInfoMap(startpage),
        'density_predict' : MetaInfoMap(startpage),
        'fock_predict' : MetaInfoMap(startpage),
        'vsolv' : MetaInfoMap(startpage),
        'abfqmmm' : MetaInfoMap(startpage),
        'hot_spot' : MetaInfoMap(startpage),
        'qmmm_switch' : MetaInfoMap(startpage),
        'core_iqmatoms' : MetaInfoMap(startpage),
        'coremask' : MetaInfoMap(startpage),
        'buffermask' : MetaInfoMap(startpage),
        'centermask' : MetaInfoMap(startpage),
        'pot_ene' : MetaInfoMap(startpage),
        'tot' : MetaInfoMap(startpage),
        'vdw' : MetaInfoMap(startpage),
        'elec' : MetaInfoMap(startpage),
        'gb' : MetaInfoMap(startpage),
        'bond' : MetaInfoMap(startpage),
        'angle' : MetaInfoMap(startpage),
        'dihedral' : MetaInfoMap(startpage),
        'vdw_14' : MetaInfoMap(startpage),
        'elec_14' : MetaInfoMap(startpage),
        'constraint' : MetaInfoMap(startpage),
        'polar' : MetaInfoMap(startpage),
        'hbond' : MetaInfoMap(startpage),
        'surf' : MetaInfoMap(startpage),
        'scf' : MetaInfoMap(startpage),
        'disp' : MetaInfoMap(startpage),
        'dvdi' : MetaInfoMap(startpage),
        'angle_ub' : MetaInfoMap(startpage),
        'imp' : MetaInfoMap(startpage),
        'cmap' : MetaInfoMap(startpage),
        'emap' : MetaInfoMap(startpage),
        'les' : MetaInfoMap(startpage),
        'noe' : MetaInfoMap(startpage),
        'pb' : MetaInfoMap(startpage),
        'rism' : MetaInfoMap(startpage),
        'ct' : MetaInfoMap(startpage),
        'amd_boost' : MetaInfoMap(startpage)
        }

    startpage.update({
        'metaNameTag'     : 'parm'
        })
    parmlist = {
        'NATOM' : MetaInfoMap(startpage),
        'NTYPES' : MetaInfoMap(startpage),
        'NBONH' : MetaInfoMap(startpage),
        'MBONA' : MetaInfoMap(startpage),
        'NTHETH' : MetaInfoMap(startpage),
        'MTHETA' : MetaInfoMap(startpage),
        'NPHIH' : MetaInfoMap(startpage),
        'MPHIA' : MetaInfoMap(startpage),
        'NHPARM' : MetaInfoMap(startpage),
        'NPARM' : MetaInfoMap(startpage),
        'NNB' : MetaInfoMap(startpage),
        'NRES' : MetaInfoMap(startpage),
        'NBONA' : MetaInfoMap(startpage),
        'NTHETA' : MetaInfoMap(startpage),
        'NPHIA' : MetaInfoMap(startpage),
        'NUMBND' : MetaInfoMap(startpage),
        'NUMANG' : MetaInfoMap(startpage),
        'NPTRA' : MetaInfoMap(startpage),
        'NATYP' : MetaInfoMap(startpage),
        'NPHB' : MetaInfoMap(startpage),
        'IFPERT' : MetaInfoMap(startpage),
        'NBPER' : MetaInfoMap(startpage),
        'NGPER' : MetaInfoMap(startpage),
        'NDPER' : MetaInfoMap(startpage),
        'MBPER' : MetaInfoMap(startpage),
        'MGPER' : MetaInfoMap(startpage),
        'MDPER' : MetaInfoMap(startpage),
        'IFBOX' : MetaInfoMap(startpage),
        'NMXRS' : MetaInfoMap(startpage),
        'IFCAP' : MetaInfoMap(startpage),
        'NUMEXTRA' : MetaInfoMap(startpage),
        'NCOPY' : MetaInfoMap(startpage)
        }

    startpage.update({
        'metaNameTag'     : 'mdout',
        'activeSections'  : ['section_single_configuration_calculation']
        })
    mddatalist = {
        'NSTEP' : MetaInfoMap(startpage),
        'TIME\(PS\)' : MetaInfoMap(startpage),
        'TEMP\(K\)' : MetaInfoMap(startpage),
        'PRESS' : MetaInfoMap(startpage),
        'Etot' : MetaInfoMap(startpage),
        'EKtot' : MetaInfoMap(startpage),
        'EPtot' : MetaInfoMap(startpage),
        'BOND' : MetaInfoMap(startpage),
        'ANGLE' : MetaInfoMap(startpage),
        'DIHED' : MetaInfoMap(startpage),
        '1-4 NB' : MetaInfoMap(startpage),
        '1-4 EEL' : MetaInfoMap(startpage),
        'VDWAALS' : MetaInfoMap(startpage),
        'EELEC' : MetaInfoMap(startpage),
        'EHBOND' : MetaInfoMap(startpage),
        'RESTRAINT' : MetaInfoMap(startpage),
        'MNDOESCF' : MetaInfoMap(startpage),
        'AINT' : MetaInfoMap(startpage),
        'EGB' : MetaInfoMap(startpage),
        'VOLUME' : MetaInfoMap(startpage),
        'EKCMT' : MetaInfoMap(startpage),
        'VIRIAL' : MetaInfoMap(startpage),
        'Density' : MetaInfoMap(startpage)
        }

    startpage.update({
        'metaNameTag'     : 'parm',
        'activeSections'  : ['x_amber_mdin_method']
        })
    extralist = {
        'flags' : MetaInfoMap(startpage),
        'box_info' : MetaInfoMap(startpage),
        'unitcell_radius' : MetaInfoMap(startpage),
        'total_memory' : MetaInfoMap(startpage),
        'file_format' : MetaInfoMap(startpage),
        'file_version' : MetaInfoMap(startpage),
        'file_date' : MetaInfoMap(startpage),
        'file_time' : MetaInfoMap(startpage)
        }

    if deflist == 'mddata':
        namelist = mddatalist
    elif deflist == 'ewald':
        namelist = ewaldlist
    elif deflist == 'qmmm':
        namelist = qmmmlist
    elif deflist == 'wt':
        namelist = wtlist
    elif deflist == 'parm':
        namelist = parmlist
    elif deflist == 'extra':
        namelist = extralist
    else:
        namelist = cntrllist
    return MapDictionary(namelist)

def get_updateDictionary(self, defname):

    startpage = {
        'metaHeader'      : '',
        'metaNameTag'     : '',
        'activeMetaNames' : [],
        'activeSections'  : ['x_amber_mdin_method']
        }

    # ---------------------------------------------------------------
    #   Definitions of meta data values for section_sampling_method
    # ---------------------------------------------------------------
    sampling = {
        'ensemble_type' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 1']],
                 'assign' : 'minimization'},
                {'test' : [['imin', '== 0'],
                           ['ntt', '== 0']],
                 'assign' : 'NVE'},
                {'test' : [['imin', '== 0'],
                           ['ntt', '> 0'],
                           ['ntt', '< 3']],
                 'assign' : 'NVT'},
                {'test' : [['imin', '== 0'],
                           ['ntp', '> 0'],
                           ['ntt', '< 3']],
                 'assign' : 'NPT'},
                {'test' : [['imin', '== 0'],
                           ['ntt', '== 3']],
                 'assign' : 'Langevin'}
                ],
            lookupdict=self.cntrlDict
            ),
        'sampling_method' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 1']],
                 'assign' : 'geometry_optimization'},
                {'test' : [['imin', '== 0']],
                 'assign' : 'molecular_dynamics'}
                ],
            lookupdict=self.cntrlDict
            ),
#            'settings_geometry_optimization' : MetaInfoMap(startpage),
#            'settings_metadynamics' : MetaInfoMap(startpage),
#            'settings_molecular_dynamics' : MetaInfoMap(startpage),
#            'settings_Monte_Carlo' : MetaInfoMap(startpage),
#        'geometry_optimization_energy_change' : MetaInfoMap(startpage,
#            depends={
#                '' : {'imin' : '1'},
#                },
#            lookupdict=self.cntrlDict
#            ),
#       'geometry_optimization_geometry_change' : MetaInfoMap(startpage),
        'geometry_optimization_method' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 1'],
                           ['ntmin', '== 0']],
                 'assign' : 'CG'},
                {'test' : [['imin', '== 1'],
                           ['ntmin', '== 1']],
                 'assign' : 'SD + CG'},
                {'test' : [['imin', '== 1'],
                           ['ntmin', '== 2']],
                 'assign' : 'XMIN'},
                {'test' : [['imin', '== 1'],
                           ['ntmin', '== 3']],
                 'assign' : 'LMOD'},
                {'test' : [['imin', '== 1'],
                           ['ntmin', '== 4']],
                 'assign' : 'LMOD'},
                ],
            lookupdict=self.cntrlDict
            ),
#       'geometry_optimization_threshold_force' : MetaInfoMap(startpage),
        'x_amber_barostat_type' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 0'],
                           ['ntp', '!= 0'],
                           ['barostat', '== 1']],
                 'assign' : 'Berendsen'},
                {'test' : [['imin', '== 0'],
                           ['ntp', '!= 0'],
                           ['barostat', '== 2']],
                 'assign' : 'Monte Carlo barostat'}
                ],
            lookupdict=self.cntrlDict,
            #autoSections=True,
            activeSections=['settings_barostat']
            ),
        'x_amber_barostat_target_pressure' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 0'],
                           ['ntp', '!= 0']],
                 'value' : 'pres0'}
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            #autoSections=True,
            activeSections=['settings_barostat']
            ),
        'x_amber_barostat_tau' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 0'],
                           ['ntp', '!= 0']],
                 'value' : 'taup'}
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            #autoSections=True,
            activeSections=['settings_barostat']
            ),
        'x_amber_integrator_type' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 0']],
                 'assign' : 'verlet'},
                {'test' : [['imin', '== 1']],
                 'assign' : 'minimization'}
                ],
            lookupdict=self.cntrlDict,
            #autoSections=True,
            activeSections=['settings_integrator']
            ),
        'x_amber_integrator_dt' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 0']],
                 'value' : 'dt'}
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            #autoSections=True,
            activeSections=['settings_integrator']
            ),
        'x_amber_number_of_steps_requested' : MetaInfoMap(startpage,
            depends=[{'value' : 'nstlim'}],
            lookupdict=self.cntrlDict,
            valtype='int',
            #autoSections=True,
            activeSections=['settings_integrator']
            ),
        'x_amber_thermostat_type' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 0'], ['ntt', '== 1']],
                 'assign' : 'Constant Temperature Scaling with weak-coupling'},
                {'test' : [['imin', '== 0'], ['ntt', '== 2']],
                 'assign' : 'Andersen-like'},
                {'test' : [['imin', '== 0'], ['ntt', '== 3']],
                 'assign' : 'Langevin'},
                {'test' : [['imin', '== 0'], ['ntt', '== 9']],
                 'assign' : 'Optimized Isokinetic Nose-Hoover'},
                {'test' : [['imin', '== 0'], ['ntt', '== 10']],
                 'assign' : 'RESPA Stochastic Isokinetic Nose-Hoover'}
                ],
            lookupdict=self.cntrlDict,
            #autoSections=True,
            activeSections=['settings_thermostat']
            ),
        'x_amber_thermostat_target_temperature' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 0'], ['ntt', '> 0']],
                 'value' : 'temp0'}
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            #autoSections=True,
            activeSections=['settings_thermostat']
            ),
        'x_amber_thermostat_tau' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 0'], ['ntt', '== 1']],
                 'value' : 'tautp'},
                {'test' : [['imin', '== 0'], ['ntt', '== 2']],
                 'value' : 'tautp'},
                {'test' : [['imin', '== 0'], ['ntt', '== 9']],
                 'value' : 'gamma_ln'},
                {'test' : [['imin', '== 0'], ['ntt', '== 10']],
                 'value' : 'gamma_ln'},
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            #autoSections=True,
            activeSections=['settings_thermostat']
            ),
        'x_amber_langevin_gamma' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['imin', '== 0'], ['ntt', '== 3']],
                 'value' : 'gamma_ln'}
                ],
            lookupdict=self.cntrlDict,
            valtype='float',
            #autoSections=True,
            activeSections=['settings_thermostat']
            ),
        'x_amber_periodicity_type' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['ntb', '== 0']],
                 'assign' : 'no periodic boundaries'},
                {'test' : [['ntb', '> 0']],
                 'assign' : 'periodic boundaries'}
                ],
            lookupdict=self.cntrlDict,
            #autoSections=True,
            activeSections=['settings_thermostat']
            )
        }

    # ------------------------------------------------------------
    #   Definitions for section_single_configuration_calculation
    # ------------------------------------------------------------
    singleconfcalc = {
        #'atom_forces_type' : MetaInfoMap(startpage,
        #    depends=[{'assign' : 'Amber Force Field'}],
        #    lookupdict=self.mddataDict
        #    ),
        'energy_correction_entropy' : MetaInfoMap(startpage),
        #'energy_current' : MetaInfoMap(startpage,
        #    depends=[{'value' : 'EPtot'}],
        #    lookupdict=self.mddataDict
        #    ),
        'energy_current' : MetaInfoMap(startpage),
        'energy_electrostatic' : MetaInfoMap(startpage,
            depends=[{'value' : 'EELEC'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kcal/mol',
            lookupdict=self.mddataDict
            ),
        'energy_free_per_atom' : MetaInfoMap(startpage),
        'energy_free' : MetaInfoMap(startpage),
        #'energy_method_current' : MetaInfoMap(startpage,
        #    depends=[{'assign' : 'Force Field'}],
        #    lookupdict=self.mddataDict
        #    ),
        'energy_T0_per_atom' : MetaInfoMap(startpage),
        'energy_total_T0_per_atom' : MetaInfoMap(startpage),
        'energy_total_T0' : MetaInfoMap(startpage),
        'energy_total' : MetaInfoMap(startpage,
            depends=[{'value' : 'Etot'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kcal/mol',
            lookupdict=self.mddataDict
            ),
        'hessian_matrix' : MetaInfoMap(startpage),
        'single_configuration_calculation_converged' : MetaInfoMap(startpage),
        'single_configuration_calculation_to_system_ref' : MetaInfoMap(startpage),
        'single_configuration_to_calculation_method_ref' : MetaInfoMap(startpage),
        'time_calculation' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_cpu1_end' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_cpu1_start' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_date_end' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_date_start' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_wall_end' : MetaInfoMap(startpage),
        'time_single_configuration_calculation_wall_start' : MetaInfoMap(startpage),
        'stress_tensor_kind' : MetaInfoMap(startpage),
        'stress_tensor_value' : MetaInfoMap(startpage)
        }

    # section_single_energy_van_der_Waals of section_single_configuration_calculation
    singlevdw = {
        'energy_van_der_Waals_value' : MetaInfoMap(startpage,
            depends=[{'value' : 'VDWAALS'}],
            lookupdict=self.mddataDict,
            valtype='float',
            unitdict=self.unitDict,
            unit='kcal/mol',
            #autoSections=True,
            activeSections=['section_energy_van_der_Waals']
            ),
        }


    # Info for section_restricted_uri
    restrictions = {
        'number_of_restricted_uri_files' : MetaInfoMap(startpage,
            value=int(1),
            autoSections=True,
            activeSections=['section_restricted_uri']
            ),
        'restricted_uri_files' : MetaInfoMap(startpage,
            subfunction=self.parameter_file_name,
            autoSections=True,
            activeSections=['section_restricted_uri']
            ),
        'restricted_uri_license' : MetaInfoMap(startpage,
            value='AMBER License',
            autoSections=True,
            activeSections=['section_restricted_uri']
            ),
        'restricted_uri_restriction' : MetaInfoMap(startpage,
            value='any access',
            autoSections=True,
            activeSections=['section_restricted_uri']
            ),
        'restricted_uri_reason' : MetaInfoMap(startpage,
            value='propriety license',
            autoSections=True,
            activeSections=['section_restricted_uri']
            ),
        'restricted_uri_issue_authority' : MetaInfoMap(startpage,
            value='AMBER',
            autoSections=True,
            activeSections=['section_restricted_uri']
            ),
        }

    # ------------------------------------------
    #   Definitions for section_frame_sequence
    # ------------------------------------------
    frameseq = {
        'frame_sequence_conserved_quantity_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'NSTEP'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_conserved_quantity_stats' : MetaInfoMap(startpage),
        'frame_sequence_conserved_quantity' : MetaInfoMap(startpage,
            depends=[{'store' : 'RESTRAINT'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kcal/mol',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_continuation_kind' : MetaInfoMap(startpage),
        'frame_sequence_external_url' : MetaInfoMap(startpage),
        'frame_sequence_kinetic_energy_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'NSTEP'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_kinetic_energy_stats' : MetaInfoMap(startpage),
        'frame_sequence_kinetic_energy' : MetaInfoMap(startpage,
            depends=[{'store' : 'EKtot'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kcal/mol',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_local_frames_ref' : MetaInfoMap(startpage),
        'frame_sequence_potential_energy_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'NSTEP'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_potential_energy_stats' : MetaInfoMap(startpage),
        'frame_sequence_potential_energy' : MetaInfoMap(startpage,
            depends=[{'store' : 'EPtot'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='kcal/mol',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_pressure_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'NSTEP'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_pressure_stats' : MetaInfoMap(startpage),
        'frame_sequence_pressure' : MetaInfoMap(startpage,
            depends=[{'store' : 'PRESS'}],
            valtype='float',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_temperature_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'NSTEP'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_temperature_stats' : MetaInfoMap(startpage),
        'frame_sequence_temperature' : MetaInfoMap(startpage,
            depends=[{'store' : 'TEMP\(K\)'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='Kelvin',
            lookupdict=self.mddataDict
            ),
        'frame_sequence_time' : MetaInfoMap(startpage,
            depends=[{'store' : 'TIME\(PS\)'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='pico-second',
            lookupdict=self.mddataDict
            ),
        'x_amber_frame_sequence_volume_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'NSTEP'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_amber_frame_sequence_volume' : MetaInfoMap(startpage,
            depends=[{'store' : 'VOLUME'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='Angstrom**3',
            lookupdict=self.mddataDict
            ),
        'x_amber_frame_sequence_density_frames' : MetaInfoMap(startpage,
            depends=[{'store' : 'NSTEP'}],
            valtype='int',
            lookupdict=self.mddataDict
            ),
        'x_amber_frame_sequence_density' : MetaInfoMap(startpage,
            depends=[{'store' : 'Density'}],
            valtype='float',
            unitdict=self.unitDict,
            unit='1.0/Angstrom**3',
            lookupdict=self.mddataDict
            ),
        #'frame_sequence_to_sampling_ref' : MetaInfoMap(startpage),
        'geometry_optimization_converged' : MetaInfoMap(startpage,
            value=self.minConverged
            ),
        'number_of_conserved_quantity_evaluations_in_sequence' : MetaInfoMap(startpage,
            value=(lambda x: np.array(x).flatten().shape[0] if x is not None else None)(
                self.metaStorage.fetchAttr(
                    {'frame_sequence_conserved_quantity_frames' : None}
                    )['frame_sequence_conserved_quantity_frames']['val']
                )
            ),
        'number_of_frames_in_sequence' : MetaInfoMap(startpage,
            value=(lambda x: np.array(x).flatten().shape[0] if x is not None else None)(
                self.metaStorage.fetchAttr(
                    {'frame_sequence_potential_energy_frames' : None}
                    )['frame_sequence_potential_energy_frames']['val']
                )
            ),
        'number_of_kinetic_energies_in_sequence' : MetaInfoMap(startpage,
            value=(lambda x: np.array(x).flatten().shape[0] if x is not None else None)(
                self.metaStorage.fetchAttr(
                    {'frame_sequence_kinetic_energy_frames' : None}
                    )['frame_sequence_kinetic_energy_frames']['val']
                )
            ),
        'number_of_potential_energies_in_sequence' : MetaInfoMap(startpage,
            value=(lambda x: np.array(x).flatten().shape[0] if x is not None else None)(
                self.metaStorage.fetchAttr(
                    {'frame_sequence_potential_energy_frames' : None}
                    )['frame_sequence_potential_energy_frames']['val']
                )
            ),
        'number_of_pressure_evaluations_in_sequence' : MetaInfoMap(startpage,
            value=(lambda x: np.array(x).flatten().shape[0] if x is not None else None)(
                self.metaStorage.fetchAttr(
                {'frame_sequence_pressure_frames' : None}
                )['frame_sequence_pressure_frames']['val']
                )
            ),
        'number_of_temperatures_in_sequence' : MetaInfoMap(startpage,
            value=(lambda x: np.array(x).flatten().shape[0] if x is not None else None)(
                self.metaStorage.fetchAttr(
                    {'frame_sequence_temperature_frames' : None}
                    )['frame_sequence_temperature_frames']['val']
                )
            ),
        #'previous_sequence_ref' : MetaInfoMap(startpage)
        }

    # ----------------------------------
    #   Definitions for section_system
    # ----------------------------------
    # section_system
    sec_system = {
        #'topology_ref' : MetaInfoMap(startpage),
        'atom_velocities' : MetaInfoMap(startpage),
        'configuration_raw_gid' : MetaInfoMap(startpage),
        'local_rotations' : MetaInfoMap(startpage),
        'number_of_atoms' : MetaInfoMap(startpage,
            depends=[{'value' : 'NATOM'}],
            valtype='int',
            lookupdict=self.parmDict
            ),
        'number_of_sites' : MetaInfoMap(startpage),
        'number_of_symmetry_operations' : MetaInfoMap(startpage),
        'reduced_symmetry_matrices' : MetaInfoMap(startpage),
        'reduced_symmetry_translations' : MetaInfoMap(startpage),
        'SC_matrix' : MetaInfoMap(startpage),
        'spacegroup_3D_choice' : MetaInfoMap(startpage),
        'spacegroup_3D_hall' : MetaInfoMap(startpage),
        'spacegroup_3D_international' : MetaInfoMap(startpage),
        'spacegroup_3D_number' : MetaInfoMap(startpage),
        'spacegroup_3D_origin_shift' : MetaInfoMap(startpage),
        'spacegroup_3D_pointgroup' : MetaInfoMap(startpage),
        'spacegroup_3D_std_lattice' : MetaInfoMap(startpage),
        'spacegroup_3D_std_positions' : MetaInfoMap(startpage),
        'spacegroup_3D_std_types' : MetaInfoMap(startpage),
        'spacegroup_3D_trasformation_matrix' : MetaInfoMap(startpage),
        'spacegroup_3D_wyckoff' : MetaInfoMap(startpage),
        'symmorphic' : MetaInfoMap(startpage),
        'system_name' : MetaInfoMap(startpage,
            subfunction=self.topology_system_name,
            #functionbase=self
            ),
        'time_reversal_symmetry' : MetaInfoMap(startpage)
        }

    # section_configuration_core of section_system
    configuration_core = {
        'number_of_electrons' : MetaInfoMap(startpage,
            value=0,
            valtype='float',
            ),
        'atom_labels' : MetaInfoMap(startpage,
            #subfunction=self.system_atom_labels()
            ),
        'atom_positions' : MetaInfoMap(startpage,
            #subfunction=self.system_atom_positions()
            ),
        'configuration_periodic_dimensions' : MetaInfoMap(startpage,
            depends=[
                {'test' : [['ntb', '== 0']],
                 'assign' : np.asarray([False, False, False])},
                {'test' : [['ntb', '> 0']],
                 'assign' : np.asarray([True, True, True])}
                ],
            lookupdict=self.cntrlDict,
            ),
        'embedded_system' : MetaInfoMap(startpage),
        'lattice_vectors' : MetaInfoMap(startpage,
            #subfunction=self.system_lattice_vectors()
            ),
        'simulation_cell' : MetaInfoMap(startpage,
            #subfunction=self.system_simulation_cell()
            )
        }

    # section_spacegroup_3D_operation of section_system
    spacegroup_op = {
            'spacegroup_3D_rotation' : MetaInfoMap(startpage),
            'spacegroup_3D_translation' : MetaInfoMap(startpage)
        }

    # section_system_to_system_refs of section_system
    sys_to_sys = {
            'system_to_system_kind' : MetaInfoMap(startpage),
            'system_to_system_ref' : MetaInfoMap(startpage)
        }

    # --------------------------------------------------------
    #   Definitions of meta data values for section_topology
    # --------------------------------------------------------
    # section_topology of section_run
    topology = {
        'atom_to_molecule' : MetaInfoMap(startpage,
            subfunction=self.topology_atom_to_mol,
            #functionbase=self
            ),
        'molecule_to_molecule_type_map' : MetaInfoMap(startpage),
        'number_of_topology_atoms' : MetaInfoMap(startpage,
            depends=[{'value' : 'NATOM'}],
            valtype='int',
            lookupdict=self.parmDict
            ),
        'number_of_topology_molecules' : MetaInfoMap(startpage,
            subfunction=self.topology_num_topo_mol,
            valtype='int',
            ),
        'topology_force_field_name' : MetaInfoMap(startpage,
            value='AMBER Force Field',
            )
        }

    # section_atom_type of section_topology
    atom_type = {
        'atom_type_charge' : MetaInfoMap(startpage),
        'atom_type_mass' : MetaInfoMap(startpage),
        'atom_type_name' : MetaInfoMap(startpage)
        }

    # section_constraint of section_topology
    constraint = {
        'constraint_atoms' : MetaInfoMap(startpage),
        'constraint_kind' : MetaInfoMap(startpage),
        'constraint_parameters' : MetaInfoMap(startpage),
        'number_of_atoms_per_constraint' : MetaInfoMap(startpage),
        'number_of_constraints' : MetaInfoMap(startpage)
        }

    # section_interaction of section_topology
    interaction = {
        'interaction_atoms' : MetaInfoMap(startpage),
        'interaction_kind' : MetaInfoMap(startpage),
        'interaction_parameters' : MetaInfoMap(startpage),
        'number_of_atoms_per_interaction' : MetaInfoMap(startpage),
        'number_of_interactions' : MetaInfoMap(startpage)
        }

    # -------------------------------------------------------------
    #   Definitions of meta data values for section_molecule_type
    # -------------------------------------------------------------
    # section_molecule_type of section_topology
    mol_type = {
        'molecule_type_name' : MetaInfoMap(startpage),
        'number_of_atoms_in_molecule' : MetaInfoMap(startpage),
        'settings_atom_in_molecule' : MetaInfoMap(startpage)
        }

    # section_molecule_constraint of section_molecule_type
    mol_constraint = {
        'molecule_constraint_atoms' : MetaInfoMap(startpage),
        'molecule_constraint_kind' : MetaInfoMap(startpage),
        'molecule_constraint_parameters' : MetaInfoMap(startpage),
        'number_of_atoms_per_molecule_constraint' : MetaInfoMap(startpage),
        'number_of_molecule_constraints' : MetaInfoMap(startpage)
        }

    # section_molecule_interaction of section_molecule_type
    mol_interaction = {
        'molecule_interaction_atoms' : MetaInfoMap(startpage),
        'molecule_interaction_kind' : MetaInfoMap(startpage),
        'molecule_interaction_parameters' : MetaInfoMap(startpage),
        'number_of_atoms_per_molecule_interaction' : MetaInfoMap(startpage),
        'number_of_molecule_interactions' : MetaInfoMap(startpage)
        }

    # section_atom_in_molecule of section_molecule_type
    atom_in_mol = {
        'atom_in_molecule_charge' : MetaInfoMap(startpage),
        'atom_in_molecule_name' : MetaInfoMap(startpage),
        'atom_in_molecule_to_atom_type_ref' : MetaInfoMap(startpage)
        }

    if defname == 'system':
        dictionary = systemDict
    elif defname == 'topology':
        dictionary = topology
    elif defname == 'singleconfcalc':
        dictionary = singleconfcalc
    elif defname == 'restrictions':
        dictionary = restrictions
    elif defname == 'frameseq':
        dictionary = frameseq
    elif defname == 'atom_type':
        dictionary = atom_type
    elif defname == 'molecule_type':
        dictionary = moltypeDict
    elif defname == 'interaction':
        dictionary = interaction
    elif defname == 'sampling':
        dictionary = sampling
    elif defname == 'singlevdw':
        dictionary = singlevdw
    else:
        dictionary = singleconfcalclist
    return MapDictionary(dictionary)

def set_excludeList(self):
    """Sets the exclude list for x_amber

    Returns:
        the list of names
    """
    excludelist = [
        'x_amber_mdin_verbatim_writeout',
        'x_amber_dumm_text',
        'x_amber_dummy',
        'x_amber_mdin_wt'
        ]
    excludelist.extend(['x_amber_mdin_file_%s' % fileNL.lower() for fileNL in self.fileDict.keys()])
    excludelist.extend(['x_amber_mdin_%s' % cntrlNL.lower() for cntrlNL in self.cntrlDict.keys()])
    excludelist.extend(['x_amber_mdin_%s' % ewaldNL.lower() for ewaldNL in self.ewaldDict.keys()])
    excludelist.extend(['x_amber_mdin_%s' % qmmmNL.lower() for qmmmNL in self.qmmmDict.keys()])
    return excludelist

def set_includeList():
    """Sets the include list for x_amber

    Returns:
        the list of names
    """
    includelist = []
    return includelist

def getList_MetaStrInDict(sourceDict):
    """Returns a list that includes all meta name
       strings for the given dictionary.
       Meta name strings are not actual meta names but
       used as the keywords in the parsing.
    """
    return [sourceDict[item].matchStr for item in sourceDict]

def getDict_MetaStrInDict(sourceDict):
    """Returns a dict that includes all meta name
       strings and corresponding values for the given dictionary.
       Meta name strings are not actual meta names but
       used as the keywords in the parsing.
    """
    newDict = {}
    for key, value in sourceDict.items():
        newDict.update({sourceDict[key].matchStr : key})
    return newDict


