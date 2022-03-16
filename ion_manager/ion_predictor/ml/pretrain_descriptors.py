
import pandas as pd
import numpy as np
import joblib
import random
from ..chem_utils.Autodescriptor import MordredDescriptor
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
import yaml

random.seed(0)


def dump(
        settings
):
    """
    automatically calculate descriptors for training
    from smiles data (path) to descriptor data (descriptor path)


    """

    path = settings["pretrain_smiles_path"]
    num_learning_molecules = settings["num_learning_molecules"]
    descriptor_path = settings["descriptor_path"]
    smiles_path = settings["smiles_path"]

    # read smiles list
    smiles_list = list(pd.read_csv(path)["SMILES"].values)

    if num_learning_molecules == -1:
        num_learning_molecules = len(smiles_list)

    # select random ones
    smiles_list = random.sample(smiles_list, num_learning_molecules)

    #########
    # calc descriptors

    # init descriptor module
    desc_calcualtor = MordredDescriptor()
    desc_calcualtor.dict_mode = False

    # calcilaion. it takes time
    descriptor_list = [desc_calcualtor.calc(sm) for sm in tqdm(smiles_list)]

    # scaling
    scaler = StandardScaler()
    descriptor_array = scaler.fit_transform(np.array(descriptor_list))

    # save
    joblib.dump(descriptor_array, descriptor_path, compress=9)
    joblib.dump(smiles_list, smiles_path, compress=9)
