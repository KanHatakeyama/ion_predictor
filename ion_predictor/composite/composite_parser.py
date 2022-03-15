import pandas as pd
from ion_predictor.composite.DataUtility import calc_wt_ratio_from_mol_wt, slash_to_list, flatten_list, unnest_dataframe,sort_compounds_by_weight_ratio
import yaml
import re
from tqdm import tqdm

from joblib import Memory
memory = Memory("cache", verbose=0)


setting_path="settings.yaml"
with open(setting_path) as file:
    settings= yaml.safe_load(file)
MAX_SMILES=settings["max_smiles"]

def process_composite_dataframe(data_df,compound_parser,neural_calculator):
    """
    automatically process the database
    """
    composite_dict = data_df.to_dict(orient="index")

    # calculate weight ratio
    composite_dict = {k: process_composite(
        composite_dict[k],compound_parser,neural_calculator) for k in (composite_dict)}

    # export as dataframe
    converted_df = pd.DataFrame(composite_dict).T

    unnest_list = ["SMILES_wt_list", "structureList",
                    "wt_ratio", "descriptor_list", "MWList"]
    converted_df = unnest_dataframe(
        converted_df, unnest_list, axis=0)
    
    #return converted_df


    # unnest FP
    unNest_FP_list = list(converted_df.columns[[True if re.match(
        "descriptor_list", i) else False for i in converted_df.columns]])
    rename_dict = {k: k+"_" for k in unNest_FP_list}
    converted_df = converted_df.rename(columns=rename_dict)

    converted_df = unnest_dataframe(
        converted_df, rename_dict.values(), axis=0)

    return converted_df

#@memory.cache
def process_composite(single_record,compound_parser,neural_calculator):
    comp_list = single_record["composition"].split("/")
    wt_list = [compound_parser(i)["averageWeight"]
                for i in comp_list]

    single_record = calc_wt_ratio_from_mol_wt(
        single_record, wt_list, key="composition", punct="/")

    single_record["wt_ratio"] = slash_to_list(single_record["wt_ratio"])[:MAX_SMILES]

    # load smiles information and split smiles
    smiles_list = [compound_parser(i)["SMILES"]
                    for i in comp_list]
    single_record["smiles_list"] = flatten_list(
        [s.split(".") for s in smiles_list])[:MAX_SMILES]
    
    # load smiles wt information
    SMILES_wt_list = [compound_parser(i)
                    ["wt_ratio"] for i in comp_list][:MAX_SMILES]
    
    # calculate weight ratio of each smiles in a composite
    wt_ratio= single_record["wt_ratio"]
    actual_smiles_weight_ratio=[]
    smiles_wt_list_float=[[float(j) for j in i.split("/")] for i in SMILES_wt_list]

    for wt, sm_wt_list in zip(wt_ratio,smiles_wt_list_float):
        for sm_wt in sm_wt_list:
            actual_smiles_weight_ratio.append(wt*sm_wt)

        single_record["SMILES_wt_list"]=actual_smiles_weight_ratio
    
    #sort by weight ratio
    single_record["smiles_list"],actual_smiles_weight_ratio=sort_compounds_by_weight_ratio(
        single_record["smiles_list"],actual_smiles_weight_ratio)
    
    #normalize ratio
    single_record["SMILES_wt_list"]=actual_smiles_weight_ratio

    # Fill blank SMILEs as "Blank"
    smiles_list = single_record["smiles_list"]
    single_record["smiles_list"] = [smiles_list[i] if i < len(
        smiles_list) else "Blank" for i in range(MAX_SMILES)]

    # load structure information
    single_record["structureList"] = [
        compound_parser(i)["Structure"] for i in comp_list][:MAX_SMILES]

    # MWList
    single_record["MWList"] = [
        compound_parser(i)["Mw"] for i in comp_list][:MAX_SMILES]

    # load fingerprint
    try:
        single_record["descriptor_list"] = [list(neural_calculator(i))
                                    for i in single_record["smiles_list"]][:MAX_SMILES]
    except:
        print("error", single_record["smiles_list"])

    single_record.pop("mol_ratio")

    return single_record