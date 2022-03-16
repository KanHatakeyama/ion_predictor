import base64
import io
from rdkit import Chem
from rdkit.Chem import Draw
import json
from django.utils.html import mark_safe
import networkx as nx

try:
    from PolyMolParser.dict_parse import parse_mol_text
    mol_parse = parse_mol_text
except:
    def mol_parse(k):
        return "error loading PolyMolParser"

# generate molecule image (as string) from smiles string


def smiles_to_buffer_img(sm: str, size: int = None) -> str:
    # fragmentate polymers
    if sm.find("M  END") >= 0:
        d = mol_parse(sm)
        if type(d) is str:
            return ""

        smiles_list = extract_smiles(d)
        sm = ".".join(smiles_list)

    mol = Chem.MolFromSmiles(sm)

    if size is None:
        img = Draw.MolToImage(mol)
    else:
        img = Draw.MolToImage(mol, size=(size, size))
    buffer = io.BytesIO()
    img.save(buffer, format="PNG")
    base64Img = base64.b64encode(buffer.getvalue()).decode().replace("'", "")
    return base64Img


def extract_smiles(obj):
    if type(obj) is str:
        return obj

    global smiles_list
    smiles_list = []

    def inner_loop(obj):
        global smiles_list

        if type(obj) is dict:
            for k, v in obj.items():
                if type(v) is dict:
                    inner_loop(v)

                if k == "SMILES":
                    smiles_list.append(v)
    inner_loop(obj)
    return smiles_list


def smiles_info(smiles):
    # in the case of normal smiles
    if smiles.find("M  END") == -1:
        return mark_safe(smiles)

    # in the case of polymer-type MOL
    parsed_dict = mol_parse(smiles)

    # in the case of module error
    if type(parsed_dict) is str:
        return parsed_dict

    tag = json.dumps(parsed_dict, indent=2, default=nx_encode)
    tag = tag.replace("  ", "___")
    tag = tag[1:-1]
    return (tag)


def nx_encode(object):
    if isinstance(object, nx.Graph):
        return ""
