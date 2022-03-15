from .data_converter import mol2dgl_single
from .auto_trainer import load_trained_model
from rdkit import Chem
import dgl
import os
import joblib

from joblib import Memory
memory = Memory("cache", verbose=0,mmap_mode='r')



class NeuralDescriptor:
    def __init__(self,settings):
        self.model=load_trained_model(settings)
        #self.desc_path=settings["neural_desc_path"]

        self.desc_dict={}
        #if not os.path.exists(self.desc_path) or settings["init_neural_descriptors"]:
        #    self.desc_dict={}
        #else:
        #    self.desc_dict=joblib.load(self.desc_path)
        self.calc_fingerprint=self.calc_descriptors

    def __call__(self,smiles):
        if smiles=="Blank":
            smiles="C"

        if smiles in self.desc_dict:
            return self.desc_dict[smiles]
        
        desc=self.predict(smiles)
        self.desc_dict[smiles]=desc
        #joblib.dump(self.desc_dict,self.desc_path)

        return desc

    def predict(self,smiles):
        return cache_predict(smiles,self.model)


    def calc_descriptors(self, SMILES_list):
        """
        calcualte descrtiptor list for list of smiles

        Parameters
        ---------------
        SEMILES_list: list of string
            smiles list
        
        Returns
        -------------
        desc_list: list of descs
        available_index: list of int
            list of index whose smiles can be calculated appropriately

        """
        desc_list=[]
        for smiles in SMILES_list:
            try:
                desc=self.__call__(smiles)
            except:
                desc=None

            desc_list.append(desc)

        available_index = [
            True if fp is not None else False for fp in desc_list]

        desc_list = [list(fp) for fp in desc_list if fp is not None]

        return desc_list, available_index

@memory.cache
def cache_predict(smiles,model):
    mol=Chem.MolFromSmiles(smiles)
    graphs = mol2dgl_single([mol])
    test_bg = dgl.batch(graphs)
    test_bg.set_e_initializer(dgl.init.zero_initializer)
    test_bg.set_n_initializer(dgl.init.zero_initializer)
    hidden = model(test_bg).detach().numpy()
    return hidden[0]

