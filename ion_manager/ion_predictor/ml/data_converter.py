
from dgl import DGLGraph
import torch
import os
from rdkit import Chem


ELEM_LIST = ['C', 'N', 'O', 'S', 'F', 'Si', 'P', 'Cl', 'Br','Na', 'Ca', 'I', 'B', 'K', 'H',  'unknown']
ATOM_FDIM = len(ELEM_LIST) + 6 + 5 + 1
MAX_ATOMNUM =60
BOND_FDIM = 5 
MAX_NB = 10
PAPER = os.getenv('PAPER', False)

def onek_encoding_unk(x, allowable_set):
    if x not in allowable_set:
        x = allowable_set[-1]
    return [x == s for s in allowable_set]


def atom_features(atom):
    return (onek_encoding_unk(atom.GetSymbol(), ELEM_LIST)
            + onek_encoding_unk(atom.GetDegree(), [0,1,2,3,4,5])
            + onek_encoding_unk(atom.GetFormalCharge(), [-1,-2,1,2,0])
            + [atom.GetIsAromatic()])

def bond_features(bond):
    bt = bond.GetBondType()
    return (torch.Tensor([bt == Chem.rdchem.BondType.SINGLE, bt == Chem.rdchem.BondType.DOUBLE, bt == Chem.rdchem.BondType.TRIPLE, bt == Chem.rdchem.BondType.AROMATIC, bond.IsInRing()]))

def mol2dgl_single(mols):
    cand_graphs = []
    n_nodes = 0
    n_edges = 0
    bond_x = []

    for mol in mols:
        n_atoms = mol.GetNumAtoms()
        n_bonds = mol.GetNumBonds()
        g = DGLGraph()       
        #g=dgl.graph()
        nodeF = []
        for i, atom in enumerate(mol.GetAtoms()):
            assert i == atom.GetIdx()
            nodeF.append(atom_features(atom))
        g.add_nodes(n_atoms)

        bond_src = []
        bond_dst = []
        for i, bond in enumerate(mol.GetBonds()):
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            begin_idx = a1.GetIdx()
            end_idx = a2.GetIdx()
            features = bond_features(bond)

            bond_src.append(begin_idx)
            bond_dst.append(end_idx)
            bond_x.append(features)
            bond_src.append(end_idx)
            bond_dst.append(begin_idx)
            bond_x.append(features)
        g.add_edges(bond_src, bond_dst)
        g.ndata['h'] = torch.Tensor(nodeF)
        cand_graphs.append(g)
    return cand_graphs