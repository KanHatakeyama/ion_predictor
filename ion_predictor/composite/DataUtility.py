"""
utility functions to process composite database

"""

from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd
import numpy as np
import re


def flatten_list(nested_list):
    return [e for inner_list in nested_list for e in inner_list]


def change_dict_key(d, old_key, new_key, default_value=None):
    d[new_key] = d.pop(old_key, default_value)


def is_nan(v):
    return v != v


def list_to_slash(ratio_list, punct="/"):
    ratio_list = [str(i) for i in ratio_list]
    return punct.join(ratio_list)


def slash_to_list(ratio_list, mode="float"):
    ratio_list = str(ratio_list).split("/")

    if mode == "float":
        return [float(i) for i in ratio_list]
    else:
        return [i for i in ratio_list]


def normalize_ratio(s):

    ratio_list = slash_to_list(s)
    tot = sum(ratio_list)
    ratio_list = [str(i/tot) for i in ratio_list]

    return list_to_slash(ratio_list)


def calc_wt_ratio_from_mol_wt(single_record, mol_wt_list, key="SMILES", punct="."):

    # if molar ratio is recorded, calculate weight ratio
    if not is_nan(single_record["mol_ratio"]):

        mol_list = str(single_record["mol_ratio"]).split("/")
        mol_list = [float(i) for i in mol_list]
        wt_list = [str(i*j) for i, j in zip(mol_list, mol_wt_list)]
        single_record["wt_ratio"] = "/".join(wt_list)

    # single component
    if single_record[key].find(punct) < 0:
        single_record["wt_ratio"] = "1"

    # multiple component
    else:
        # if no information is given about watio, just put 1/1/1..
        if is_nan(single_record["wt_ratio"]):
            semiles_count = single_record[key].count(punct)+1
            array_of_one = ["1"]*semiles_count
            single_record["wt_ratio"] = "/".join(array_of_one)

        # normalize weight ratio
        single_record["wt_ratio"] = normalize_ratio(single_record["wt_ratio"])

    return single_record


# calculate molar weight of compounds excluding Mg and Ca
def calc_mol_weight(smiles):
    mol = Chem.MolFromSmiles(smiles)
    try:
        wt = Descriptors.MolWt(mol)
        f_charge = Chem.rdmolops.GetFormalCharge(mol)
    except:
        wt = 100
        f_charge = 0
    # lithium is often omitted in the compound database. add the weight of lithium
    if f_charge < 0:
        wt = wt-6.941*f_charge

    # Ca and Mg is often used just to reperesent the repeating units
    wt = wt-40.078*smiles.count("Ca")
    wt = wt-24.305*smiles.count("Mg")

    return wt


def sort_compounds_by_weight_ratio(comp_list, wt_list):

    ordering_list = [[k, v] for k, v in zip(wt_list, comp_list)]
    ordering_list.sort(reverse=True)
    a, b = list(zip(*ordering_list))
    return b, a


def unnest_dataframe(df, explode, axis):
    if axis == 1:
        idx = df.index.repeat(df[explode[0]].str.len())
        df1 = pd.concat([
            pd.DataFrame({x: np.concatenate(df[x].values)}) for x in explode], axis=1)
        df1.index = idx

        return df1.join(df.drop(explode, 1), how='left')
    else:
        df1 = pd.concat([
            pd.DataFrame(df[x].tolist(), index=df.index).add_prefix(x) for x in explode], axis=1)
        return df1.join(df.drop(explode, 1), how='left')


def get_column_names(df, name):
    return list(df.columns[[True if re.match(name, i) else False for i in df.columns]])
