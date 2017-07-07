import setup_paths
import numpy as np
import logging
import json
import os
import re
from collections import namedtuple

class MetaInfoMap(object):
    """Map cache values to meta info
    """
    def __init__(self, activeInfo=False, defaultValue=None, metaName=None, 
            metaNameTag=None, metaInfoType=None, metaInfoValues=None, 
            metaInfoDepends=None, activeSections=None, subFunction=None):
        self.activeInfo = activeInfo
        self.defaultValue = defaultValue
        self.metaName = metaName
        self.metaNameTag = metaNameTag
        self.metaInfoType = metaInfoType
        self.metaInfoValues = metaInfoValues
        self.metaInfoDepends = metaInfoDepends
        self.activeSections = activeSections
        self.subFunction = subFunction

class FileInfoMap(dict):
    """Map cache values to meta info
    """
    activeInfo=False
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
                    self[v.metaHeader + '_' + v.metaNameTag + '_' + v.metaName] = v

        if kwargs:
            for k, v in kwargs.items():
                if (isinstance(v, FileInfoMap) or
                    isinstance(v, MetaInfoMap)):
                    if v.metaTranslate:
                        v.metaName = v.nameTranslate(k)
                    else:
                        v.metaName = k
                self[v.metaHeader + '_' + v.metaNameTag + '_' + v.metaName] = v

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

def metaNameConverter(keyName):
    newName = keyName.lower().replace(" ", "").replace("-", "")
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
        'activeMetaNames' : []
        'activeSections'  : ['x_amber_section_input_output_files']
        }
    namelist = {
            'MDIN'   : FileInfoMap(startpage),
            'MDOUT'  : FileInfoMap(startpage), 
            'INPCRD' : FileInfoMap(startpage, activeInfo=True, fileFormat=['.inpcrd', '.ncrst']), 
            'PARM'   : FileInfoMap(startpage, activeInfo=True, fileFormat=['.prmtop']),
            'RESTRT' : FileInfoMap(startpage),
            'REFC'   : FileInfoMap(startpage),
            'MDVEL'  : FileInfoMap(startpage, activeInfo=True),
            'MDFRC'  : FileInfoMap(startpage, activeInfo=True),
            'MDEN'   : FileInfoMap(startpage),
            'MDCRD'  : FileInfoMap(startpage, activeInfo=True, fileFormat=['.netcdf', '.mdcrd']),
            'MDINFO' : FileInfoMap(startpage),
            'MTMD'   : FileInfoMap(startpage),
            'INPDIP' : FileInfoMap(startpage),
            'RSTDIP' : FileInfoMap(startpage), 
            'INPTRA' : FileInfoMap(startpage)
            }
    return MapDictionary(namelist)

def get_nameList(deflist):
    """Loads namelist data of AMBER.

    Args:
        deflist: name list definition (cntrl/ewald/qmmm/wt).

    Returns:
        the list of namelists
    """
    cntrllist = [
        'imin',         'nmropt',        'ntx',       'irest',     'ntxo',    'ntpr',          
        'ntave',        'ntwr',          'iwrap',     'ntwx',      'ntwv',    'ntwf',
        'ntwe',         'ioutfm',        'ntwprt',    'idecomp',   'ibelly',  'ntr',   
        'restraint_wt', 'restraintmask', 'bellymask', 'maxcyc',    'ncyc',    'ntmin',
        'dx0',          'drms',          'nstlim',    'nscm',      't',       'dt',
        'nrespa',       'ntt',           'temp0',     'temp0les',  'tempi',   'ig',
        'tautp',        'gamma_ln',      'vrand',     'vlimit',    'nkija',   'idistr',
        'sinrtau',      'ntp',           'barostat',  'mcbarint',  'pres0',   'comp',
        'taup',         'csurften',      'gamma_ten', 'ntc',       'tol',     'jfastw',
        'noshakemask',  'ivcap',         'fcap',      'outcap',    'xcap',    'ycap',
        'zcap',         'iscale',        'noeskp',    'ipnlty',    'mxsub',   'scalm',
        'pencut',       'tausw',         'iemap',     'gammamap',  'ntf',     'ntb',
        'dielc',        'cut',           'fswitch',   'nsnb',      'ipol',    'ifqnt',
        'igb',          'irism',         'ievb',      'iamoeba',   'lj1264',  'efx',
        'efy',          'efz',           'efn',       'efphase',   'effreq'
        ]

    wtlist = [
        'istep1',       'istep2',      'value1',       'value2',       'iinc',     'imult'
        ]
    
    ewaldlist = [
        'nfft1',       'nfft2',      'nfft3',       'order',       'verbose',     'ew_type',
        'dsum_tol',    'rsum_tol',   'mlimit(1)',   'mlimit(2)',   'mlimit(3)',   'ew_coeff',
        'nbflag',      'skinnb',     'nbtell',      'netfrc',      'vdwmeth',     'eddmet',
        'eedtbdns',    'column_fft', 'ips',         'raips',       'mipsx',       'mipsy',
        'mipsz',       'mipso',      'gridips',     'dvbips',      'frameon',     'chngmask',
        'indmeth',     'diptol',     'maxiter',     'dipmass',     'diptau',      'irstdip',
        'scaldip'
        ]

    qmmmlist = [
        'qm_theory',      'dftb_slko_path',       'dftb_disper',  'dftb_3rd_order', 'dftb_chg',
        'dftb_telec',     'dftb_maxiter',         'qmcharge',     'spin',           'qmqmdx',
        'verbosity',      'tight_p_conv',         'scfconv',      'pseudo_diag',    'diag_routine',
        'printcharges',   'print_eigenvalues',    'qxd',          'parameter_file', 'peptide_corr',
        'itrmax',         'ntpr',                 'grms_tol',     'ndiis_attempts', 'ndiis_matrices', 
        'vshift',         'errconv',              'qmmm_int',     'qmmask',         'qmcut', 
        'lnk_dis',        'dftb_telec_step',      'fockp_d1',     'fockp_d2',       'fockp_d3', 
        'fockp_d4',       'pseudo_diag_criteria', 'damp',         'kappa',          'min_heavy_mass', 
        'r_switch_hi',    'r_switch_lo',          'iqmatoms',     'qmgb',           'lnk_atomic_no',  
        'lnk_method',     'printbondorders',      'buffercharge', 'printdipole',    'qmshake', 
        'qmmmrij_incore', 'qmqm_erep_incore',     'qm_ewald',     'qm_pme',         'kmaxqx', 
        'kmaxqy',         'kmaxqz',               'ksqmaxsq',     'adjust_q',       'density_predict', 
        'fock_predict',   'vsolv',                'abfqmmm',      'hot_spot',       'qmmm_switch', 
        'core_iqmatoms',  'coremask',             'buffermask',   'centermask',     'pot_ene', 
        'tot',            'vdw',                  'elec',         'gb',             'bond', 
        'angle',          'dihedral',             'vdw_14',       'elec_14',        'constraint', 
        'polar',          'hbond',                'surf',         'scf',            'disp', 
        'dvdi',           'angle_ub',             'imp',          'cmap',           'emap', 
        'les',            'noe',                  'pb',           'rism',           'ct', 
        'amd_boost'
        ]
    
    parmlist = [
        'NATOM',    'NTYPES', 'NBONH',  'MBONA',  'NTHETH', 'MTHETA',
        'NPHIH',    'MPHIA',  'NHPARM', 'NPARM',  'NNB',    'NRES',
        'NBONA',    'NTHETA', 'NPHIA',  'NUMBND', 'NUMANG', 'NPTRA',
        'NATYP',    'NPHB',   'IFPERT', 'NBPER',  'NGPER',  'NDPER',
        'MBPER',    'MGPER',  'MDPER',  'IFBOX',  'NMXRS',  'IFCAP',
        'NUMEXTRA', 'NCOPY'
        ]

    mddatalist = [
        'NSTEP',    'TIME', 'TEMP',  'PRESS',  'Etot', 'EKtot',
        'EPtot',    'BOND',  'ANGLE', 'DIHED',  '1-4 NB',    '1-4 EEL',
        'VDWAALS',    'EELEC', 'EHBOND',  'RESTRAINT', 'MNDOESCF', 'AINT',
        'EGB',    'VOLUME',   'EKCMT', 'VIRIAL',  'Density'
        ]

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
    else:
        namelist = cntrllist
    return namelist

def set_excludeList(self):
    """Sets the exclude list for x_amber

    Returns:
        the list of names
    """
    excludelist = [
        'x_amber_mdin_verbatim_writeout',
        'x_amber_dumm_text',
        'x_amber_dummy',
        ]
    excludelist.extend(['x_amber_mdin_file_%s' % fileNL.lower() for fileNL in self.fileDict.keys()])
    excludelist.extend(['x_amber_mdin_%s' % cntrlNL.lower() for cntrlNL in get_nameList('cntrl')])
    excludelist.extend(['x_amber_mdin_%s' % ewaldNL.lower() for ewaldNL in get_nameList('ewald')])
    excludelist.extend(['x_amber_mdin_%s' % qmmmNL.lower() for qmmmNL in get_nameList('qmmm')])
    return excludelist

def set_includeList():
    """Sets the include list for x_amber

    Returns:
        the list of names
    """
    includelist = [
        'x_amber_mdin_wt'
        ]
    return includelist


