"""
this is a wrapper class to easily calculate fingerpritns or descriptors by RDKit
20200905 modif lib to show errors with broken molecules
20210504 add jrgui module, change interface
"""

#https://github.com/KanHatakeyama/polySMILES/blob/main/polysmiles/AutoDescriptor.py

# import rdkit library
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
import pandas as pd
import numpy as np
from numpy import inf
from rdkit.Avalon.pyAvalonTools import GetAvalonFP
#from .JRWrapper import JRWrapper
import os
import zlib
import base64
import joblib
from joblib import Memory


memory = Memory("cache", verbose=0)


#make mol object from smiles
def mol_from_smiles(smiles,assert_smiles=True):
    m=Chem.MolFromSmiles(smiles)
    
    if assert_smiles:
        assert m, ("invalid smiles!", smiles)    
    else:
        if m is None:
            print("failed to purse: ", smiles )
    return m

#draw chemical structure from smiles
def draw_SMILES(smiles):
    """
    smiles: SMILES (string)
    return: image
    """
    m=mol_from_smiles(smiles)
    return Draw.MolToImage(m) 


def auto_correct_descs(descs):
    descs=np.nan_to_num(descs)
    descs[descs > 10**5] = 0
    descs[descs < -10**5] = 0
    return descs

def check_smiles_list(sm_list):
    for smiles in sm_list:
        mol_from_smiles(smiles,assert_smiles=False)

#descriptor class
class RDKitDescriptors:
    def __init__(self,auto_correct=True,dict_mode=True):
        self.desc_list = [desc_name[0] for desc_name in Descriptors.descList]
        self.calculator = MoleculeDescriptors.MolecularDescriptorCalculator(self.desc_list)
        self.desc_list = ["RDKit_desc_"+desc_name for desc_name in self.desc_list]
        self.auto_correct=auto_correct
        self.dict_mode=dict_mode
        
    def _desc_calc(self,smiles):
        m=mol_from_smiles(smiles)
        return self.calculator.CalcDescriptors(m)
    
    #calc descriptors    
    def calc(self,smiles):
        """
        smiles: smiles
        dict_mode: if true, return type will be a dict, otherwise, list
        """

            #calculate newly
        descs=self._desc_calc(smiles)

        if self.auto_correct:
            descs=auto_correct_descs(descs)      

        desc_dict={k:v for k,v in zip(self.desc_list,descs)}

        if self.dict_mode:
            return desc_dict
        else:
            return descs

    #calc descriptors from smiles list
    def calc_list(self,ls,pandas_mode=True):
        """
        ls: list of smiles
        pandas_mode: if true, return type will be dataframe, otherwise list
        """
        temp_mode=self.dict_mode
        self.dict_mode=False
        res_list=[self.calc(i) for i in ls]
        self.dict_mode=temp_mode
        
        if pandas_mode:
            df=pd.DataFrame(res_list)
            df.columns=self.desc_list
            return df
        
        return res_list
    

    def auto_calc(self,arg,pandas_mode=True):
        if type(arg) is type(""):
            return self.calc(arg)
        elif type(arg) is type([]):
            return self.calc_list(arg,pandas_mode=pandas_mode)
        else:
            assert False,("unexpected type: ",type(arg))
            
            
#calculate fingerprnt
#default FP function
def default_FP_Func(smiles):
    m=mol_from_smiles(smiles)
    fp=GetAvalonFP(m) 
    fp=[int(i) for i in fp]
    return fp

class Fingerprint(RDKitDescriptors):
    def __init__(self,fp_func=default_FP_Func,dict_mode=True):
        super(RDKitDescriptors, self).__init__()
        self.calculator=fp_func
        self.fp_len=len(fp_func("C"))
        self.desc_list=["FP_"+str(i) for i in range(self.fp_len)]
        self.auto_correct=False
        self.dict_mode=dict_mode
        
    def _desc_calc(self,m):
        return self.calculator(m)
        
"""        
class GroupContMethod(RDKitDescriptors):
    def __init__(self,dict_mode=True,auto_correct=True):
        super(RDKitDescriptors, self).__init__()
        self.dict_mode=dict_mode
        self.calculator=JRWrapper()
        self.auto_correct=auto_correct
        self.desc_list=list(self._desc_calc("C",dict_mode=True).keys())
        
    def _desc_calc(self,smiles,dict_mode=False):
        res=self.calculator.compute_phys_properties(smiles)
        
        if dict_mode:
            return res
        else:
            return list(res.values())
     """   
        
class AutoDescriptor():
    def __init__(self,calculators=[
        #GroupContMethod(),
        RDKitDescriptors()]):
        self.calculators=calculators
        self.descriptor_names=list(self.__call__("C").columns)
        self.descriptor_names.remove("SMILES")
        
    def __call__(self,smiles_list):
        if type(smiles_list) is type(""):
            smiles_list=[smiles_list]
        elif type(smiles_list) is type([]):
            pass
        else:
            assert False,("unexpected type: ",type(smiles_list))
        
        integ_pd=pd.DataFrame()
        integ_pd["SMILES"]=smiles_list
        for num,calculator in enumerate(self.calculators):
            df=(calculator.auto_calc(smiles_list))
            integ_pd=pd.concat([integ_pd,df],axis=1)
        
        return integ_pd
    
class AutoFingerprint:
    def __init__(self):
        self.ad=AutoDescriptor([Fingerprint()])
        self.descriptor_names=list(self.__call__("C").columns)
        self.descriptor_names.remove("SMILES")
    def __call__(self,smiles_list):
        return self.ad(smiles_list)


class MordredDescriptor(RDKitDescriptors):
    def __init__(self,dict_mode=True,auto_correct=True,ignore_3D=True):
        from mordred import Calculator, descriptors
        super(RDKitDescriptors, self).__init__()
        self.dict_mode=dict_mode
        self.calculator=Calculator(descriptors,ignore_3D=ignore_3D)
        self.auto_correct=auto_correct
        self.desc_list=list(self.calculator.pandas([mol_from_smiles("C")]).columns)
        self.desc_list = ["Mordred_desc_"+desc_name for desc_name in self.desc_list]
        
    def _desc_calc(self,smiles,dict_mode=False):
        res=mord_calc(smiles,self.calculator)
        if dict_mode:
            return res
        else:
            return list(res)

#use cache for mordred
@memory.cache
def mord_calc(smiles,calculator):
    m=mol_from_smiles(smiles)
    res=list(calculator(m))
    res=np.array(res).astype(float)
    return res

