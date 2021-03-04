from rdkit.Chem import Lipinski
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import MolWt
from gendock.functional_groups import load_functional_groups

"""

Rules:
Max H-bond donors = 5
Max H-bond acceptors = 10
Max mass = 500
Max logP = 5

"""


def check_lipinski(mol):
    fgs = load_functional_groups()
    h_donors = Lipinski.NumHDonors(mol.rdmol)
    h_acceptors = Lipinski.NumHAcceptors(mol.rdmol)
    log_p = MolLogP(mol.rdmol)
    wt = MolWt(mol.rdmol)
    if h_donors <= 5 and h_acceptors <= 5 and log_p < 5:
        if wt >= 450:
            mol.join(fgs['terminal_fg'].get_random())
            return True, False
        else:
            return True, False
    else:
        return False, False
