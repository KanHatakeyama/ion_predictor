import pandas as pd
from .CompoundParser import CompoundParser 
from .composite_parser import process_composite_dataframe
import copy
import numpy as np
from ..ml.NeuralDescriptor import NeuralDescriptor

def clean_df(original_df):
    parsed_df=copy.copy(original_df)
    for c in ["smiles_list","composition"]:
        if c in list(parsed_df.columns):
            parsed_df=parsed_df.drop(c,axis=1)
        

    parsed_df["Conductivity"]=np.log10(parsed_df["Conductivity"].astype(float))
    #parsed_df=parsed_df.dropna(how='all')
    #parsed_df=parsed_df.loc[:,~(parsed_df.nunique()==1)]


    parsed_df["Temperature"]=parsed_df["Temperature"].astype(float)
    parsed_df["inorg_contain_ratio(wt)"]=parsed_df["inorg_contain_ratio(wt)"].astype(float)

    return parsed_df


def load_ion_excel(settings,
                    compound_sheet_name="CompoundDatabase",
                    composite_sheet_name="Database",
                    compound_df=None,
                    composite_df=None):

    compound_load_list = ['ID',  'SMILES', 'mol_ratio', 'wt_ratio', 'Mw',
                        'Mn', 'MwMn', 'Polymn_deg', 'Structure']

    composite_load_list = ["ID", "composition", "mol_ratio", "wt_ratio", "inorg_name",
                        "inorg_contain_ratio(wt)", "Temperature", "Conductivity"]


    excel_path=settings["excel_path"]

    if compound_df is None:
        compound_df = pd.read_excel(excel_path, sheet_name=compound_sheet_name)
            
    if composite_df is None:
        composite_df = pd.read_excel(excel_path, sheet_name=composite_sheet_name)

    compound_df=compound_df[compound_load_list]
    composite_df=composite_df[composite_load_list]

    return auto_df_load(compound_df,composite_df,settings)


#load compound data
def auto_df_load(compound_df,composite_df,settings):
    #init calculators
    compound_parser=CompoundParser(compound_df)
    neural_calculator=NeuralDescriptor(settings)

    #load main database
    raw_df=process_composite_dataframe(composite_df,compound_parser,neural_calculator)

    parsed_df=clean_df(raw_df)

    return parsed_df