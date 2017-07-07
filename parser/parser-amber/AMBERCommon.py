import setup_paths
import numpy as np
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion import convert_unit
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from AMBERDictionary import get_nameList, get_fileListDict, set_excludeList, set_includeList
import logging
import json
import os
import re

PARSER_INFO_DEFAULT = {
        'name':'amber-parser', 
        'version': '0.0.1'
}

META_INFO_PATH = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 
    "../../../../nomad-meta-info/meta_info/nomad_meta_info/amber.nomadmetainfo.json"))

LOGGER = logging.getLogger("nomad.AMBERParser")

def get_metaInfo(self):
    metaInfoEnv, warnings = loadJsonFile(filePath=META_INFO_PATH, 
                                         dependencyLoader=None, 
                                         extraArgsHandling=InfoKindEl.ADD_EXTRA_ARGS, 
                                         uri=None)
    metaInfoEnv = set_metaInfoEnv(metaInfoEnv, 'mdin', get_nameList('cntrl'), 'C', ["x_amber_mdin_method"])
    metaInfoEnv = set_metaInfoEnv(metaInfoEnv, 'mdin', get_nameList('ewald'), 'C', ["x_amber_mdin_method"])
    metaInfoEnv = set_metaInfoEnv(metaInfoEnv, 'mdin', get_nameList('qmmm'), 'C', ["x_amber_mdin_method"])
    metaInfoEnv = set_metaInfoEnv(metaInfoEnv, 'parm', get_nameList('parm'), 'C', ["x_amber_mdin_method"])
    metaInfoEnv = set_section_metaInfoEnv(metaInfoEnv, 'section', ['input_output_files'], 'type_section', False, ["section_run"])
#    metaInfoEnv = set_section_metaInfoEnv(metaInfoEnv, 'mdin_file', 
#                                          ['input_output_files'], 'type_abstract_document_content', 
#                                          False, ["x_amber_section_input_output_files"])
    metaInfoEnv = setDict_metaInfoEnv(metaInfoEnv, self.fileDict)
    metaInfoEnv = set_metaInfoEnv(metaInfoEnv, 'mdout', get_nameList('mddata'), 'C', ["section_single_configuration_calculation"])
#    metaInfoEnv = set_metaInfoEnv(metaInfoEnv, 'mdout', get_nameList('mddata'), 'C', ["x_amber_mdout_single_configuration_calculation"])
    metaInfoEnv = set_metaInfoEnv(metaInfoEnv, 'parm', 
                                  [ 'flags',           'box_info',
                                    'unitcell_radius', 'total_memory',
                                    'file_format',     'file_version', 
                                    'file_date',       'file_time'], 
                                  'C', ["x_amber_mdin_method"])
    return metaInfoEnv

def set_section_metaInfoEnv(infoEnv, metaNameTag, newList, listTypStr, repeatingSection, supraNames):
    """Modifies meta info data.

    Args:
        metaInfoEnv: meta info environment json type data.

    Returns:
        metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
    """
    for newName in newList:
        newName = newName.lower().replace(" ", "").replace("-", "")
        if 'x_amber_%s_%s' % (metaNameTag, newName) not in infoEnv.infoKinds:
            infoEnv.addInfoKindEl(InfoKindEl(
                description='auto generated section meta info data',
                name='x_amber_%s_%s' % (metaNameTag, newName),
                kindStr=listTypStr,
                repeats=repeatingSection,
                superNames=supraNames)) 

    return infoEnv

def setDict_metaInfoEnv(infoEnv, nameDict):
    """Modifies meta info data.

    Args:
        metaInfoEnv: meta info environment json type data.
        nameDict: dictionary for name info and data.

    Returns:
        metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
    """
    for keyName in nameDict.keys():
        if '%s' % (keyName) not in infoEnv.infoKinds:
            infoEnv.addInfoKindEl(InfoKindEl(
                name='%s' % (keyName),
                description='auto generated meta info data',
                dtypeStr=nameDict[keyName].metaInfoType,
                shape=[],
                superNames=nameDict[keyName].activeSections)) 

    return infoEnv

def set_metaInfoEnv(infoEnv, metaNameTag, newList, listTypStr, supraNames):
    """Modifies meta info data.

    Args:
        metaInfoEnv: meta info environment json type data.

    Returns:
        metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
    """
    for newName in newList:
        newName = newName.lower().replace(" ", "").replace("-", "")
        if 'x_amber_%s_%s' % (metaNameTag, newName) not in infoEnv.infoKinds:
            infoEnv.addInfoKindEl(InfoKindEl(
                name='x_amber_%s_%s' % (metaNameTag, newName),
                description='auto generated meta info data',
                dtypeStr=listTypStr,
                shape=[],
                superNames=supraNames)) 

    return infoEnv

def write_mdin(self, backend, metaInfoEnv, metaNameStart, valuesDict, writeCheck, location, logger):
    """Writes the values of the mdin file.

    Write the last occurrence of a keyword, i.e. [-1], since aims uses the last occurrence of a keyword.

    ATTENTION
    backend.superBackend is used here instead of only the backend to write the JSON values,
    since this allows to bybass the caching setting which was used to collect the values for processing.
    However, this also bypasses the checking of validity of the metadata name by the backend.
    The scala part will check the validity nevertheless.

    Args:
        backend: Class that takes care of writing and caching of metadata.
        metainfoenv: Loaded metadata.
        valuesdict: Dictionary that contains the cached values of a section.
        writecheck: Boolean that determines whether the keywords related to the ?? should be written.
        location: A string that is used to specify where more than one setting was found.
        logger: Logging object where messages should be written to.
    """
    # list of excluded metadata for writeout
    # namelist values at metdata are only needed to detect 
    # the type of simulation from cntrl not to writeout to backend.
    excludelist = set_excludeList(self)
    includelist = set_includeList()
    # write settings
    for k,v in valuesDict.items():
        if (k.startswith('x_amber_mdin_') or 
            k.startswith('x_amber_mdout_') or
            k.startswith('x_amber_parm_')):
            if (k     in excludelist and 
                k not in includelist):
                continue
            # default writeout
            else:
                # convert keyword values of mdin which are strings to lowercase for consistency
                if isinstance(v[-1], str):
                    value = v[-1].lower()
                else:
                    value = v[-1]
                backend.superBackend.addValue(k, value)
    #if writecheck:
    #    return 1

    # distinguish between cntrl namelist 
    if metaNameStart == 'x_amber_mdin':
        Write = False
    elif metaNameStart == 'x_amber_mdinout':
        Write = True
    else:
        logger.error("Unknown metaNameStart %s in function in %s. Please correct." % (metaNameStart, os.path.basename(__file__)))
        return

class AMBERParserBase(object):
    """Base class for Amber parsers"""
    def __init__(self,cachingLevelForMetaName=None, coverageIgnoreList=None,
                 re_program_name=None):
        self.re_program_name = re_program_name
        self.fileDict = get_fileListDict()
        self.parserInfo = PARSER_INFO_DEFAULT.copy()
        self.metaInfoEnv = get_metaInfo(self)
#        self.cachingLevelForMetaName = {
#                               'x_amber_trajectory_file_detect': CachingLevel.Cache,
#                               'x_amber_geometry_optimization_cdetect': CachingLevel.Cache,
#                               'x_amber_section_MD_detect': CachingLevel.Ignore,
#                               'x_amber_single_configuration_calculation_detect': CachingLevel.Cache,
#                              }
        self.coverageIgnoreList = [
            # ignore empty lines
            r"\s*",
            # table separators
            #r"^\s*[=%-]+\s*$",
            #r"^\s*%\s*%\s*$",
        ]
        self.coverageIgnore = None

    def parse(self):
        self.coverageIgnore = re.compile(r"^(?:" + r"|".join(self.coverageIgnoreList) + r")$")
        mainFunction(mainFileDescription=self.mainFileDescription(), 
                     metaInfoEnv=self.metaInfoEnv, 
                     parserInfo=self.parserInfo,
                     cachingLevelForMetaName=self.cachingLevelForMetaName,
                     superContext=self)

    def adHoc_amber_program_name(self, parser):
        if self.re_program_name is not None:
            if not self.re_program_name.match(
                    parser.lastMatch['x_amber_program_name']):
                raise Exception(
                    "mainFile program name was: %s, unsuited for %s" % (
                        parser.lastMatch['x_amber_program_name'],
                        type(self).__name__))

    def mainFileDescription(self):
        # assemble matchers and submatchers
        return SM(name='Root',
            startReStr="",
            forwardMatch=True,
            weak=True,
            subMatchers=self.build_subMatchers()
            ) # END Root

    def build_subMatchers(self):
        return []

#class MetaInfoStorageBase(object):
#    """Base class for meta info storage"""
#    def __init__(self,cachingLevelForMetaName=None, coverageIgnoreList=None):
#        self.re_program_name = re_program_name


