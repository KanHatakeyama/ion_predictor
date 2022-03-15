import pandas as pd
import numpy as np
from .DataUtility import flatten_list, calc_mol_weight, calc_wt_ratio_from_mol_wt
from .DataUtility import list_to_slash, change_dict_key, slash_to_list, is_nan, sort_compounds_by_weight_ratio


class CompoundDatabase:
    """
    a special class to analyze the compound database recorded in the excel sheet
    """

    def __init__(self, path_to_excel, fp, sheet_name="CompoundDatabase"):
        """
        Parameters
        -----------------
        path_to_excel: string
            path the the excel sheet
        fp: Fingerprint class
            a class to convert smiles to fingerprint
        sheet_name: string
            sheet name of the compound database
        """

        load_list = ['ID',  'SMILES', 'mol_ratio', 'wt_ratio', 'Mw',
                     'Mn', 'MwMn', 'Polymn_deg', 'Structure']

        self.data_df = pd.read_excel(
            path_to_excel, sheet_name=sheet_name)[load_list]
        self.fp = fp
        self.n_bits = len(self.fp.calc_fingerprint("C")[0][0])

    def calc_fp(self):
        """
        calc fingerptint of the compounds recorded in the "SMILES" column
        """
        smiles_list = list(set(self.data_df["SMILES"]))

        # split SMILES with dot  e.g., [Mg]OCC[Mg].[Mg]C(C)(C(O)=O)C[Mg] to ["[Mg]OCC[Mg]"" , "[Mg]C(C)(C(O)=O)C[Mg]"]
        split_smiles_list = [i.split(".") for i in smiles_list]
        split_smiles_list = list(set(flatten_list(split_smiles_list)))

        _, available_index = self.fp.calc_fingerprint(split_smiles_list)
        split_smiles_list = np.array(split_smiles_list)[available_index]
        fp_list, _ = self.fp.calc_fingerprint(split_smiles_list)

        self.fp_dict = {k: v for k, v in zip(split_smiles_list, fp_list)}
        self.fp_dict["Blank"] = [0]*self.n_bits

        mol_weight_list = [calc_mol_weight(i) for i in split_smiles_list]
        self.mol_weight_dict = {k: v for k,
                                v in zip(smiles_list, mol_weight_list)}

    def calc_weight_info(self, single_record):
        """
        calculate molecular weight of the compounds.

        Parameters
        ------------------
        single_record: data series
            one line of the database

        Returns
        ------------------
        single_record: data series
            one processed line of the database

        """

        # calculate weight ratio
        smiles_list = single_record["SMILES"].split(".")
        mol_wt_list = [calc_mol_weight(sm) for sm in smiles_list]

        single_record = calc_wt_ratio_from_mol_wt(single_record, mol_wt_list)

        # calculate average molecular weight of the compound
        single_record["averageWeight"] = sum(
            [i*j for i, j in zip(slash_to_list(single_record["wt_ratio"]), mol_wt_list)])

        # calculate Mw
        if not is_nan(single_record["Mw"]):
            pass
        else:

            # calculate from Mn and MwMn or polymn_deg
            v = single_record["Mn"]*single_record["MwMn"]
            single_record["Mw"] = v if not is_nan(
                v) else single_record["Polymn_deg"]*single_record["averageWeight"]

        single_record.pop("mol_ratio")
        single_record.pop("Mn")
        single_record.pop("MwMn")
        single_record.pop("Polymn_deg")

        if len(smiles_list) > 1:
            smiles_list, wt_list = sort_compounds_by_weight_ratio(
                smiles_list, slash_to_list(single_record["wt_ratio"]))
            single_record["SMILES"] = list_to_slash(smiles_list, punct=".")
            single_record["wt_ratio"] = list_to_slash(wt_list)

        return single_record

    def process_compound_database(self):
        """
        automatically process the compound database
        """
        self.comp_dict = self.data_df.to_dict(orient="index")

        # calculate weight ratio, etc
        self.comp_dict = {k: self.calc_weight_info(
            self.comp_dict[k]) for k in self.comp_dict}

        # change dict ID to compound ID
        for i in list(self.comp_dict.keys()):
            change_dict_key(self.comp_dict, i, str(self.comp_dict[i]["ID"]))
