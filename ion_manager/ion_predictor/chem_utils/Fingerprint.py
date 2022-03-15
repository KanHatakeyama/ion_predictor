"""
wrapper class of rdkit to calculate fingerprint

"""
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Avalon.pyAvalonTools import GetAvalonFP


class Fingerprint:
    # calculate fingerprint from SMILES

    def __init__(self):
        self.fpFunc = self.fingerprint

    def fingerprint(self, mol):
        return GetAvalonFP(mol, nBits=512)
        #return AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)

    def SMILES_to_bit(self, SMILES):
        """
        calculate bit for smiles
        
        Parameters
        ---------------
        SEMILES: string
            smiles
        
        Returns
        -------------
        bit: list of bit
            e.g, (10101011..)
        """
        try:
            mol = Chem.MolFromSmiles(SMILES)
            fp = self.fpFunc(mol)
            bit = fp.ToBitString()
            bit = [1 if i == "1" else 0 for i in bit]
        except:
            bit = None
            print("error ", SMILES)

        return bit

    def calc_fingerprint(self, SMILES_list):
        """
        calcualte bit list for list of smiles

        Parameters
        ---------------
        SEMILES_list: list of string
            smiles list
        
        Returns
        -------------
        fingerprint_list: list of (list of bit)
            e.g, [(10101011..) , (111...),...]     
        available_index: list of int
            list of index whose smiles can be calculated appropriately

        """
        fingerprint_list = [self.SMILES_to_bit(sm) for sm in SMILES_list]
        available_index = [
            True if fp is not None else False for fp in fingerprint_list]

        fingerprint_list = [fp for fp in fingerprint_list if fp is not None]

        return fingerprint_list, available_index

    def convert_bit_list_to_string(self, fp):
        return "".join([str(i) for i in fp])


def get_Tanimoto(fp1, fp2):
    """
    Calculate Tanimoto score

    Parameters
    ----------------
    fp1: bit array
        fingerprint 1
    fp2: bit array
        fingerprint 2
    
    Returns
    ---------------
    return: float
        Tanimoto score of fp1 and fp2

    """

    one_count = 0
    zero_count = 0
    different_count = 0
    if len(fp1) != len(fp2):
        print("caution! different FP length")
        return 0

    for i, j in zip(fp1, fp2):
        if i == j:
            if i == 1:
                one_count += 1
            else:
                zero_count += 1
        else:
            pass
            # different_count+=1
    if (len(fp1)-zero_count) == 0:
        return 0
    else:
        return one_count / (len(fp1)-zero_count)


#not used
def get_normal_similarity(fp1, fp2):
    one_count = 0
    zero_count = 0
    different_count = 0
    if len(fp1) != len(fp2):
        print("caution! different FP length")
        return 0

    for i, j in zip(fp1, fp2):
        if i == j:
            if i == 1:
                one_count += 1
            else:
                zero_count += 1
        else:
            pass
            different_count += 1

    return (one_count+zero_count) / (len(fp1))
