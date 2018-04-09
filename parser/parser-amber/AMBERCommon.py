import setup_paths
import numpy as np
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.unit_conversion.unit_conversion import convert_unit
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from AMBERDictionary import get_unitDict, get_nameListDict, get_fileListDict, set_excludeList, set_includeList
from MetaInfoStorage import COMMON_META_INFO_PATH, PUBLIC_META_INFO_PATH
import MetaInfoStorage as mStore
import logging
import json
import os
import re

PARSERNAME = "AMBER"
PROGRAMNAME = "amber"
PARSERVERSION = "0.0.3"
PARSERMETANAME = PARSERNAME.lower()
PARSERTAG = 'x_' + PARSERMETANAME

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
    metaInfoEnv = set_section_metaInfoEnv(metaInfoEnv, 'section', ['input_output_files'], 'type_section', False, ["section_run"])
    metaInfoEnv = setDict_metaInfoEnv(metaInfoEnv, self.fileDict)
    metaInfoEnv = setDict_metaInfoEnv(metaInfoEnv, self.cntrlDict)
    metaInfoEnv = setDict_metaInfoEnv(metaInfoEnv, self.ewaldDict)
    metaInfoEnv = setDict_metaInfoEnv(metaInfoEnv, self.qmmmDict)
    metaInfoEnv = setDict_metaInfoEnv(metaInfoEnv, self.parmDict)
    metaInfoEnv = setDict_metaInfoEnv(metaInfoEnv, self.mddataDict)
    metaInfoEnv = setDict_metaInfoEnv(metaInfoEnv, self.extraDict) 
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

class AMBERParserBase(object):
    """Base class for Amber parsers"""
    def __init__(self,cachingLevelForMetaName=None, coverageIgnoreList=None,
                 re_program_name=None):
        self.metaStorage = mStore.Container('section_run')
        self.metaStorageRestrict = mStore.Container('section_restricted_uri')
        exclude_dict = { 
            'section_run' : [
            'section_processor_info', 
            'section_processor_log', 
            'section_springer_material',
            'section_repository_info'
            ]}
        jsonmetadata = mStore.JsonMetaInfo(
                COMMON_META_INFO_PATH, 
                PUBLIC_META_INFO_PATH,
                META_INFO_PATH
                )
        self.metaStorage.build(jsonmetadata, 'section_run', exclude_dict)
        self.metaStorageRestrict.build(jsonmetadata, 'section_restricted_uri', exclude_dict)
        self.re_program_name = re_program_name
        self.unitDict = get_unitDict('si')
        self.fileDict = get_fileListDict()
        self.cntrlDict = get_nameListDict('cntrl')
        self.ewaldDict = get_nameListDict('ewald')
        self.qmmmDict = get_nameListDict('qmmm')
        self.wtDict = get_nameListDict('wt')
        self.parmDict = get_nameListDict('parm')
        self.mddataDict = get_nameListDict('mddata')
        self.extraDict = get_nameListDict('extra')
        self.parserInfo = PARSER_INFO_DEFAULT.copy()
        self.metaInfoEnv = get_metaInfo(self)
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


