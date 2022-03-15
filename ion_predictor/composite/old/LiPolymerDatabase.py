import pandas as pd
from .DataUtility import calc_wt_ratio_from_mol_wt, slash_to_list, flatten_list, unnest_dataframe,sort_compounds_by_weight_ratio
import re
import numpy as np

# maximum number of smiles to be processed in one composite
MAX_SMILES = 6


class LiPolymerDatabase:
    """
    a special class to process the lithium conducting polymer database:
    """

    def __init__(self, path_to_excel, compound_database,
                 ignore_temperature_dependency=True,
                 temperature_range=(-190, 340)):
        """
        Parameters
        ---------------
        path_to_excel: string
            path to the excel database
        compound_database: CompoundDatabase class
            please initialize compound database class before this class
        ignore_temperature_dependency: bool
            if true, temperature dependency of a single conductor is ignored (treated as duplication)
        temperature_range: (float,float)
            range of temperature to be processed
        """
        self.load_list = ["ID", "composition", "mol_ratio", "wt_ratio", "inorg_name",
                          "inorg_contain_ratio(wt)", "Temperature", "Conductivity"]

        self.data_df = pd.read_excel(path_to_excel, sheet_name="Database")[
            self.load_list]
        self.compound_database = compound_database

        self.temperature_range = temperature_range
        self.data_df = self.data_df[self.data_df["Temperature"]
                                    < self.temperature_range[1]]
        self.data_df = self.data_df[self.data_df["Temperature"]
                                    > self.temperature_range[0]]

        if ignore_temperature_dependency:
            self.data_df = self.data_df.drop_duplicates(
                subset=["composition", "mol_ratio", "wt_ratio", "inorg_name", "inorg_contain_ratio(wt)"])

    def calc_compound_weight_ratio(self, single_record):
        """
        calcualte the weight ratio of a compound in a composite

        Parameters
        --------------
        single_record: data series
            a single line of the database (dataframe) 

        Returns
        --------------
        single_record: data series
            a processed data
        """
        comp_list = single_record["composition"].split("/")
        wt_list = [self.compound_database.comp_dict[i]["averageWeight"]
                   for i in comp_list]

        single_record = calc_wt_ratio_from_mol_wt(
            single_record, wt_list, key="composition", punct="/")

        single_record["wt_ratio"] = slash_to_list(single_record["wt_ratio"])[:MAX_SMILES]

        # load smiles information and split smiles
        smiles_list = [self.compound_database.comp_dict[i]["SMILES"]
                       for i in comp_list]
        single_record["smiles_list"] = flatten_list(
            [s.split(".") for s in smiles_list])[:MAX_SMILES]
        
        # load smiles wt information
        SMILES_wt_list = [self.compound_database.comp_dict[i]
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
        #norm_ratio=[actual_smiles_weight_ratio[i]/actual_smiles_weight_ratio[i+1] for i in range(len(actual_smiles_weight_ratio)-1)]
        #temp_list=[0]
        #temp_list.extend(list(np.log10(norm_ratio)))
        #single_record["SMILES_wt_list"]=temp_list
        single_record["SMILES_wt_list"]=actual_smiles_weight_ratio
    
        #print (single_record["SMILES_wt_list"])
        
        # Fill blank SMILEs as "Blank"
        smiles_list = single_record["smiles_list"]
        single_record["smiles_list"] = [smiles_list[i] if i < len(
            smiles_list) else "Blank" for i in range(MAX_SMILES)]

        # load structure information
        single_record["structureList"] = [
            self.compound_database.comp_dict[i]["Structure"] for i in comp_list][:MAX_SMILES]

        # MWList
        single_record["MWList"] = [
            self.compound_database.comp_dict[i]["Mw"] for i in comp_list][:MAX_SMILES]

        # load fingerprint
        try:
            single_record["fp_list"] = [self.compound_database.fp_dict[i]
                                        for i in single_record["smiles_list"]][:MAX_SMILES]
        except:
            print("error", single_record["smiles_list"])

        single_record.pop("mol_ratio")

        return single_record

    def process_database(self):
        """
        automatically process the database
        """
        self.DBDict = self.data_df.to_dict(orient="index")

        # calculate weight ratio
        self.DBDict = {k: self.calc_compound_weight_ratio(
            self.DBDict[k]) for k in self.DBDict}

        # export as dataframe
        self.converted_df = pd.DataFrame(self.DBDict).T

        unnest_list = ["SMILES_wt_list", "structureList",
                      "wt_ratio", "fp_list", "MWList"]
        self.converted_df = unnest_dataframe(
            self.converted_df, unnest_list, axis=0)

        # unnest FP
        unNest_FP_list = list(self.converted_df.columns[[True if re.match(
            "fp_list", i) else False for i in self.converted_df.columns]])
        rename_dict = {k: k+"_" for k in unNest_FP_list}
        self.converted_df = self.converted_df.rename(columns=rename_dict)

        self.converted_df = unnest_dataframe(
            self.converted_df, rename_dict.values(), axis=0)
