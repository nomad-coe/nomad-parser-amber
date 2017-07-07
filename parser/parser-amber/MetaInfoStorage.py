import setup_paths
import numpy as np
from nomadcore.unit_conversion.unit_conversion import convert_unit
import logging
import json
import os
import re
from collections import namedtuple

COMMON_META_INFO_PATH = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 
    "../../../../nomad-meta-info/meta_info/nomad_meta_info/common.nomadmetainfo.json"))

PUBLIC_META_INFO_PATH = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 
    "../../../../nomad-meta-info/meta_info/nomad_meta_info/public.nomadmetainfo.json"))

class Container(object):
    """The container class for nested data storage
    """
    def __init__(self, *args):
        self.Name = []
        for arg in args:
            if isinstance(arg, str):
                self.Name.append(arg)
        self.Active = False
        self.OpenBackend = False
        self.CloseBackend = False
        self.gid = None
        self.Containers = []
        self.References = []
        self.ReferencedFrom = []
        self.Storage = None
        self.Indent = '   '
        self.Color = None
        self.PrintOnlyActive = None
        self.PrintOnlyNames = None

    def add(self, *args):
        for arg in args:
            if isinstance(arg, Container):
                self.Containers.append(arg) 
            if isinstance(arg, Storage):
                if self.Storage:
                    self.Storage(arg.__dict__)
                else:
                    self.Storage = arg
            if isinstance(arg, dict):
                if self.Storage:
                    if isinstance(self.Storage, dict):
                        self.Storage.update(arg)
                    elif isinstance(self.Storage, Storage):
                        self.Storage(arg)
                else:
                    self.Storage = Storage(arg)
            if isinstance(arg, JsonMetaInfo):
                for item in arg.jsonList:
                    self.add(item) 

    def build(self, metadata, startsection, umaskdict):
        if isinstance(metadata, JsonMetaInfo):
            attrdict = metadata.attributes(startsection)
            if attrdict:
                self.add(attrdict)
            childs = metadata.siblings(startsection)
            if startsection in umaskdict.keys():
                excludes = umaskdict[startsection]
            else:
                excludes = []
            for section in childs:
                if section not in excludes:
                    newContainer = Container(section)
                    newContainer.build(metadata, section, umaskdict)
                    if newContainer.Storage is None:
                        newContainer.add()
                    self.add(newContainer)
            return True
        else:
            return False

    def populate(self, metadata, startsection, umaskdict, updatedict):
        if isinstance(metadata, JsonMetaInfo):
            self.build(metadata, startsection, umaskdict)
            self.update(updatedict)

    def updaterefs(self, *args):
        for arg in args:
            if isinstance(arg, Container):
                self.References.extend(arg.Name) 

    def update(self, *args):
        for arg in args:
            if isinstance(arg, dict):
                if arg["startSection"]:
                    if self.Name in arg["startSection"]:
                        self.accumulateValues(self, arg)
                    else:
                        if self.Containers:
                            for module in self.Containers:
                                module.update(arg)
                else:
                    self.accumulateValues(self, arg)

    def updateBackend(self, backend, startsection=None, autoopenclose=False):
        if startsection:
            if (self.Name in startsection or 
                self.Name == startsection):
                if autoopenclose:
                    with self.autosection(self.Name):
                        self.updateBackendValues(self, backend)
                else:
                    self.updateBackendValues(self, backend)
            else:
                if self.Containers:
                    for module in self.Containers:
                        module.updateBackend(backend, startsection, autoopenclose)
        else:
            if autoopenclose:
                with self.autosection(self.Name):
                    self.updateBackendValues(self, backend)
            else:
                self.updateBackendValues(self, backend)

    def autosection(self, sectionname):
        self.opensection(sectionname)
        yield self.gid
        self.closesection(sectionname)

    def opensection(self, backend, sectionname):
        self.gid = backend.openSection(sectionname)
        yield self.gid

    def closesection(self, backend, sectionname):
        backend.closeSection(sectionname, self.gid)

    def fetchAttr(self, resdict):
        for item in resdict:
            if self.Storage:
                if item in self.Storage.__dict__:
                    resdict.update({item: self.Storage.__dict__[item]})
                else:
                    if self.Containers:
                        for module in self.Containers:
                            resdict.update(module.fetchAttr(resdict))
            else:
                if self.Containers:
                    for module in self.Containers:
                        resdict.update(module.fetchAttr(resdict))
        return resdict

    def updateBackendValues(self, backend):
        if self.Storage:
            self.updateBackendStorage()
            self.Active = False
        if self.Containers:
            for module in self.Containers:
                module.updateBackendValues(backend)

    def accumulateValues(self, *args):
        for arg in args:
            if isinstance(arg, dict):
                if self.Storage:
                    self.accumulateDict(arg["dictionary"])
                if self.Containers:
                    for module in self.Containers:
                        module.accumulateValues(arg)
                if "activeSections" in arg:
                    if self.Name in arg["activeSections"]:
                        self.Active = True
                if "muteSections" in arg:
                    if self.Name in arg["muteSections"]:
                        self.Active = False

    def checkUpdateValue(self, item):
        updateValue = False
        if item.lookupdict:
            for depdict in item["depends"]:
                depmeet = 0
                for depk, depv in depdict:
                    if depk in item.lookupdict:
                        if item.lookupdict[depk].value==depv:
                            depmeet += 1
                if depmeet == len(depdict.keys()):
                    updateValue = True
        else:
            for depdict in item["depends"]:
                depmeet = 0
                depcopydict = depdict.copy()
                attrdict = self.fetchAttr(depcopydict)
                for depk, depv in depdict:
                    if depk in attrdict:
                        if attrdict[depk]==depv:
                            depmeet += 1
                if depmeet == len(depdict.keys()):
                    updateValue = True
        return updateValue

    def updateBackendStorage(self, backend):
        for itemk in self.Storage.__dict__:
            if self.Storage.__dict__[itemk]["act"]:
                self.Storage.__dict__[itemk]["act"] = False
                value = self.Storage.__dict__[itemk]["val"]
                if isinstance(value, np.ndarray):
                    backend.addArrayValues(itemk, value)
                elif isinstance(value, (list, tuple)):
                    backend.addArrayValues(itemk, np.asarray(value))
                else:
                    backend.addValue(itemk, value)

    def accumulateDict(self, checkDict):
        for itemk in checkDict:
            if itemk in self.Storage.__dict__:
                itemv = checkDict[itemk]
                updateValue = False
                if itemv["depends"]:
                    updateValue = self.checkUpdateValue(itemv)
                else:
                    updateValue = True
                if updateValue:
                    self.Storage.__dict__[itemk]["val"]=itemv["value"]
                    self.Storage.__dict__[itemk]["act"] = True
                    self.Active = True
                    if "valueSize" in itemv:
                        if "sizeMetaName" in itemv:
                            self.Storage.__dict__[itemv["sizeMetaName"]]=itemv["valueSize"]
                    if "subfunction" in itemv:
                        newValue = itemv["subfunction"](self, itemk, itemv)
                        self.Storage.__dict__[itemk["val"]]=newvalue

    def __str__(self, caller=None, decorate='', color=None, printactive=None, onlynames=None):
        string = ''
        if onlynames is None:
            if self.PrintOnlyNames:
                onlynames = self.PrintOnlyNames
        if printactive is None:
            if self.PrintOnlyActive:
                printactive = self.PrintOnlyActive
        if color:
            color = int(color) + 1
        else:
            if self.Color:
                color = int(self.Color)
        printok = False
        if printactive:
            if self.Active:
                printok = True
        else:
            printok = True
        if caller:
            if printok:
                if color:
                    string = '%s\033[9%sm`-->[' % (decorate + self.Indent, str(color%6 + 1))
                    string = string + ','.join(['%s' % (name) for name in self.Name]) + ']\n\033[0m'
                else:
                    string = '%s`-->[' % (decorate + self.Indent)
                    string = string + ','.join(['%s' % (name) for name in self.Name]) + ']\n'
        else:
            if printok:
                if color:
                    string = '%s\033[9%sm-->[' % (decorate + self.Indent, str(color%6 + 1)) 
                    string = string + ','.join(['%s' % (name) for name in self.Name]) + ']\n\033[0m'
                else:
                    string = '%s-->[' % (decorate + self.Indent) 
                    string = string + ','.join(['%s' % (name) for name in self.Name]) + ']\n'
        if printok:
            if color:
                string = string + '%s\033[9%sm|\033[0m   `.\n' % (decorate + self.Indent, str(color%6 + 1)) 
            else:
                string = string + '%s|   `.\n' % (decorate + self.Indent) 
        if self.Storage:
            for key in self.Storage.__dict__:
                printattrok = False
                if printactive:
                    if self.Storage.__dict__[key]["act"]:
                        printattrok = True
                else:
                    printattrok = True
                if printattrok and printok:
                    if color:
                        if onlynames:
                            string = string + '%s\033[9%sm|\033[0m    |__.%s\n' % (decorate 
                                    + self.Indent, str(color%6 + 1), key)
                        else:
                            string = string + '%s\033[9%sm|\033[0m    |__.%s : Active=%s Value=%s\n' % (decorate 
                                    + self.Indent, str(color%6 + 1), key, self.Storage.__dict__[key]["act"], 
                                    self.Storage.__dict__[key]["val"])
                    else:
                        if onlynames:
                            string = string + '%s|    |__.%s : Active=%s Value=%s\n' % (decorate + 
                                    self.Indent, key)
                        else:
                            string = string + '%s|    |__.%s : Active=%s Value=%s\n' % (decorate + 
                                    self.Indent, key, self.Storage.__dict__[key]["act"], 
                                    self.Storage.__dict__[key]["val"])
            if color:
                string = string + '%s\033[9%sm|\033[0m\n' % (decorate + self.Indent, str(color%6 + 1))
            else:
                string = string + '%s|\n' % (decorate + self.Indent)
        if self.Containers:
            for module in self.Containers:
                if color:
                    string = string + '%s\033[9%sm|\033[0m\n' % (decorate + self.Indent, 
                            str(color%6 + 1)) + module.__str__(self.Name, 
                            '%s\033[9%sm|\033[0m' % (decorate + self.Indent, str(color%6 + 1)), 
                            color, printactive, onlynames)
                else:
                    string = string + '%s|\n' % (decorate + self.Indent) + module.__str__(self.Name, 
                            '%s|' % (decorate + self.Indent), 
                            color, printactive, onlynames)
        return string

class Storage(dict):
    """ Storage for meta info document types.
        Sections are build by Container class
    """
    def __init__(self, *args, **kwargs):
        super(Storage, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Storage, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Storage, self).__delitem__(key)
        del self.__dict__[key]


class JsonMetaInfo(object):
    """ Json file loader for meta info data of NOMAD.
        Loads data and extracts values of items 
        with specified superNames
    """
    def __init__(self, *args):
        self.jsonList = None
        for filepath in args:
            try:
                with open(filepath, encoding="utf-8") as f:
                    jsonDict = json.load(f)
            except:
                logging.exception("Error while loading file %s" % filePath)
                raise
            typeStr = jsonDict.get("type","nomad_meta_info_1_0")
            typeRe = re.compile(r"nomad_meta_info_(?P<major>[0-9]+)_(?P<minor>[0-9]+)$")
            m = typeRe.match(typeStr)
            if not m:
                raise Exception("unexpected type '%s', expected nomad_meta_info_1_0" % typeStr)
            newJsonList = jsonDict.get("metaInfos",[])
            if self.jsonList:
                self.jsonList = self.jsonList + newJsonList
            else:
                self.jsonList = newJsonList

    def attributes(self, sectionname):
        attributes = {}
        for item in self.jsonList:
            superlist = item['superNames']
            itemname = item['name']
            try:
                kindname = item['kindStr']
            except:
                kindname = []
            try:
                dtyp = item['dtypeStr']
            except:
                dtyp = []
            try:
                size = item['shape']
            except:
                size = []
            try:
                unit = item['units']
            except:
                unit = []
            try:
                refs = item['referencedSections']
            except:
                refs = []
            if ('type_section' in kindname or 
                sectionname not in superlist):
                continue
            attrvalues = {
                    'act' : False, 
                    'val' : None, 
                    'kind': kindname, 
                    'dtyp': dtyp,
                    'unit': unit,
                    'size': size,
                    'refs': refs
                    }
            attributes.update({itemname: attrvalues})
        return attributes

    def siblings(self, sectionname):
        siblings = []
        searchList = []
        nameList = []
        for item in self.jsonList:
            superlist = item['superNames']
            itemname = item['name']
            try:
                kindname = item['kindStr']
            except:
                kindname = []
            if ('type_section' in kindname or 
                'type_abstract_document_content' in kindname): 
                if sectionname in superlist:
                    searchList.append(itemname)
                if itemname not in nameList:
                    nameList.append(itemname)
        for name in searchList:
            if (('section' in name or 'settings' in name ) and 
                set([name]) not in set(nameList) and 
                self.isparent(name)):
                siblings.append(name)
        return siblings

    def isparent(self, itemname):
        haschild = False
        for item in self.jsonList:
            if itemname in item['superNames']:
                haschild = True
        return haschild

    def rootsections(self):
        rootname = []
        searchList = []
        nameList = []
        for item in self.jsonList:
            superlist = item['superNames']
            itemname = item['name']
            try:
                kindname = item['kindStr']
            except:
                kindname = []
            if 'type_section' not in kindname:
                continue
            if not superlist:
                searchList.append(itemname)
            if itemname not in nameList:
               nameList.append(itemname)
        for name in searchList:
            if ('section' in name and set([name]) not in set(nameList)):
                rootname.append(name)
        return rootname

    def fetchdict(self, itemname, pattern):
        resDict = {}
        for item in self.jsonList:
            val = dict(item)
            itemProperty = item[itemname]
            if pattern in itemProperty:
                resDict.update({item["name"]: val})
        return resDict

if __name__ == "__main__":
    run = Container('section_run')
    exclude_dict = { 
            'section_run' : [
            'section_processor_info', 
            'section_processor_log', 
            'section_springer_material',
            'section_repository_info'
            ]}

    jsonmetadata = JsonMetaInfo(COMMON_META_INFO_PATH, PUBLIC_META_INFO_PATH)

    updateDict = {
            'startSection' : [['section_topology']],
            'muteSections' : [['section_interaction']],
            'dictionary' : {
                'molecule_constraint_atoms' : {'value' : 100, 'depends' : {}},
                'interaction_atoms' : {'value' : 10, 'depends' : {}},
                'topology_force_field_name' : {'value' : "ReaxFF", 'depends' : {}}
                }
            }
    
    run.populate(jsonmetadata, 'section_run', exclude_dict, updateDict)
    run.Color = 4
    for container in run.Containers:
        if 'section_topology' in container.Name:
            select = container
    #select.Color = 4
    #select.PrintOnlyActive = 1
    run.PrintOnlyNames = 1
    print(run)


