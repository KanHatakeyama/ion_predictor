from .DataUtility import calc_mol_weight, calc_wt_ratio_from_mol_wt
from .DataUtility import list_to_slash, slash_to_list, is_nan, sort_compounds_by_weight_ratio


class CompoundParser:
    def __init__(self,data_df):
        #self.data_df=data_df
        temp_list=data_df.to_dict(orient="records")
        self.all_comp_dict={}
        for i in temp_list:
            self.all_comp_dict[str(i["ID"])]=i

        self.cache_dict={}
    def __call__(self,compound_id):
 
        if compound_id in self.cache_dict:
            return self.cache_dict[compound_id]

        try:
            comp_data=self.all_comp_dict[compound_id]
        except:
            raise ValueError(f"{compound_id} not found!")
        process_comp_data(comp_data)

        self.cache_dict[compound_id]=comp_data
        return comp_data
        #return cache_call(compound_id,self.all_comp_dict)

"""
@memory.cache
def cache_call(compound_id,all_comp_dict):
        try:
            compound_id=int(compound_id)
        except:
            pass
        try:
            #comp_data=data_df[data_df["ID"]==compound_id].to_dict(orient="records")[0]
            comp_data=all_comp_dict[compound_id]
        except:
            raise ValueError(f"{compound_id} not found!")
        process_comp_data(comp_data)
        return comp_data

"""

def process_comp_data(comp_data):
       # calculate weight ratio
        smiles_list = comp_data["SMILES"].split(".")
        mol_wt_list = [calc_mol_weight(sm) for sm in smiles_list]

        comp_data = calc_wt_ratio_from_mol_wt(comp_data, mol_wt_list)

        # calculate average molecular weight of the compound
        comp_data["averageWeight"] = sum(
            [i*j for i, j in zip(slash_to_list(comp_data["wt_ratio"]), mol_wt_list)])

        # calculate Mw
        if not is_nan(comp_data["Mw"]):
            pass
        else:

            # calculate from Mn and MwMn or polymn_deg
            v = comp_data["Mn"]*comp_data["MwMn"]
            comp_data["Mw"] = v if not is_nan(
                v) else comp_data["Polymn_deg"]*comp_data["averageWeight"]

        comp_data.pop("mol_ratio")
        comp_data.pop("Mn")
        comp_data.pop("MwMn")
        comp_data.pop("Polymn_deg")

        if len(smiles_list) > 1:
            smiles_list, wt_list = sort_compounds_by_weight_ratio(
                smiles_list, slash_to_list(comp_data["wt_ratio"]))
            comp_data["SMILES"] = list_to_slash(smiles_list, punct=".")
            comp_data["wt_ratio"] = list_to_slash(wt_list)

