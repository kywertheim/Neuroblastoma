# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 12:41:46 2019

@author: Kenneth Y. Wertheim
"""

import json
from json import JSONEncoder
from PyNB.Environment import Environment
from PyNB.Neuroblastoma import Neuroblastoma
from PyNB.ValidationData import ValidationData

'''
Basic class to allow desired objects 'Neuroblastoma' and 'Environment' to be serialised
'''
class NBEncoder(JSONEncoder):
    def default(self, object):
        if isinstance(object, Neuroblastoma):
            return object.__dict__
        elif isinstance(object, Environment):
            return object.__dict__
        elif isinstance(object, ValidationData):
            return object.__dict__
        else:
            # call base class implementation which takes care of
            # raising exceptions for unsupported types
            return json.JSONEncoder.default(self, object)

'''
Decodes JSON dictionary into 'Neuroblastoma' and 'Environment' objects
This should scale dynamically as member variables are added/removed
'''
def as_complex(dict):
    if 'environment' in dict:
        env = Environment()
        env_dict = env.__dict__
        for key, val in dict['environment'].items():
            if key in env_dict:
                setattr(env, key, val)
        dict['environment'] = env
    if 'Neuroblastoma' in dict:
        nb_dict = Neuroblastoma().__dict__
        nbList = []
        for _nb in dict['Neuroblastoma']:
            nb = Neuroblastoma()
            for key, val in _nb.items():
                if key in nb_dict:
                    setattr(nb, key, val)
            nbList.append(nb)
        dict['Neuroblastoma'] = nbList
    if 'validation' in dict:
        vd = ValidationData()
        vd_dict = vd.__dict__
        for key, val in dict['validation'].items():
            if key in vd_dict:
                setattr(vd, key, val)
        dict['validation'] = vd
    return dict

'''
Returns config loaded from JSON at specified path
@param filePath to JSON input
@return Dictionary containing config.Environment, Neuroblastoma and run configurations
'''
def load(filePath, suppressVersionWarnings=False):
    try:
        inFile = open(filePath,"r") 
        ret = json.load(inFile, object_hook=as_complex)
        inFile.close()
        #Versioning warnings
        if 'version' in ret:
            if ret['version'][0] != Environment.VERSION()[0] or ret['version'][1] != Environment.VERSION()[1] or ret['version'][2] != Environment.VERSION()[2]:
                Environment.VERSION_WARNING = True
                if not suppressVersionWarnings:
                    print('Warning: Input file version (v%s.%s.%s) does not match runtime (v%s.%s.%s)' %(
                    ret['version'][0], ret['version'][1], ret['version'][2],
                Environment.VERSION()[0], Environment.VERSION()[1], Environment.VERSION()[2]))
            elif ret['version'][3]:
                Environment.VERSION_WARNING = True                
                if not suppressVersionWarnings:
                    print("Warning: Version from input file (%s) contains warn flag." % (filePath))
        else:
            Environment.VERSION_WARNING = True             
            if not suppressVersionWarnings:
                print("Warning: Version missing from input file: %s" % (filePath))
        return ret;
    except IOError:
        raise IOError("Could not open file '%s' for JSON import." % (filePath))


def save(filePath, environment, steps, nbList=None, validation=None , pretty=True, seed=None, distributed=None):
    #Construct the export structure
    exportState = {}
    exportState['version'] = Environment.VERSION()
    exportState['config'] = {}
    
    if isinstance(environment, Environment):
        exportState['config']['environment'] = environment
    else:
        raise ValueError("environment must be type PyNB.Environment")
        
    if seed!=None and isinstance(seed, int):
        exportState['config']['seed'] = seed
    elif seed!=None:
        raise ValueError("seed must be type int")
        
    if distributed!=None and isinstance(distributed, bool):
        exportState['config']['distributed'] = distributed
    elif seed!=None:
        raise ValueError("distributed must be type bool")
        
    if isinstance(steps, int):
        exportState['config']['steps'] = steps
    else:
        raise ValueError("steps must be type int")
    
    if nbList!=None:    
        if isinstance(nbList, list):
            for nb in nbList:
                if not isinstance(nb, Neuroblastoma):
                    raise ValueError("steps must be type [PyNB.Neuroblastoma]")
            exportState['Neuroblastoma'] = nbList
        else:
            raise ValueError("steps must be type [PyNB.Neuroblastoma]")
       
    if validation!=None:
        if not isinstance(validation, ValidationData):
            raise ValueError("validation must be type PyNB.ValidationData")
        exportState['validation'] = validation
       
    #Perform export
    try:
        outFile = open(filePath,"w+") 
        if pretty==True:
            json.dump(exportState, outFile, cls=NBEncoder, indent=4)
        else:
            json.dump(exportState, outFile, cls=NBEncoder)   
        outFile.close()
    except IOError:
        raise IOError("Could not open file '%s' for JSON export."%(filePath))