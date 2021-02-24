from rdkit.Chem import Lipinski
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Descriptors import MolWt
from gendock.functional_groups import terminal_fg

"""

Rules:
Max H-bond donors = 5
Max H-bond acceptors = 10
Max mass = 500
Max logP = 5

"""


def check_lipinski(mol):
    h_donors = Lipinski.NumHDonors(mol)
    h_acceptors = Lipinski.NumHAcceptors(mol)
    log_p = MolLogP(mol)
    wt = MolWt(mol)
    if h_donors <= 5 and h_acceptors <= 5 and log_p < 5:
        if wt >= 450:
            mol.join(terminal_fg.get_random())
            return True, False
        else:
            return True, False
    else:
        return False, False
