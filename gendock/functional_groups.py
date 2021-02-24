from rdkit.Chem import MolFromSmiles
from gendock.docking import DockMol
from random import choice


class FGLibrary:
    def __init__(self, fgs, position):
        self.position = position
        self.fgs = fgs
        return

    def get_random(self):
        return choice(self.fgs)

    def return_all(self):
        return self.fgs


s_list = [
    "[*]c1ccccc1",
    "C1=CN=C([*])C=C1",
    "C1=CN=CN=C1[*]",
    "c1([*])cccc2cccnc21",
    "c1([*])c2c(cccc2)nc3ccccc31",
    "C1([*])CCc2c1cccc2",
    "C([*])(C1)Cc2c1cccc2",
    "C(CC1[*])c2c1cccc2",
    "c1([*])c2c(cccc2)ccc1",
    "c1c2c(cccc2)ccc1[*]",
    "c1([*])ccc2c(c1)cccc2",
    "c1ccc2c(c1)cc([*])cc2",
    "c1([*])cc2c(cccc2)cc1",
    "NC1=C(C2=C(C=C1)N[*])C(C3=C(O)C=CC(O)=C3C2=O)=O",
    "NC1=C(C(C2=C(C=C3)O)=O)C(C(C2=C3O)=O)=C(N[*])C=C1",
    "O=C1C2=C(C=NN2)N=C([*])N1",
    "O=C1C2=C(C=NN2[*])N=CN1",
    "C1([*])=NC(C=NN2)=C2C=N1",
    "C1=NC(C([*])=NN2)=C2C=C1",
    "C1=CC(C=NN2[*])=C2C=C1",
    "C1([*])=NC(C=CN2)=C2C=C1",
]


start_fg = FGLibrary([DockMol(MolFromSmiles(i)) for i in s_list], 'start')
# non-terminal groups

nt_list = [
    "c1c([*])cccc1",
    "c1cc([*])ccc1",
    "c1ccc([*])cc1",
    "OC[*]",
    "C(=O)OC[*]",
    "CC(=O)C[*]",
    "C[*]",
    "CS(=O)(=O)O[*]",
    "C/C=C/[*]",
    "CC#C[*]",
    "CS[*]",
    "CC(=O)N[*]",
    "CN[*]",
    "C1OC1[*]",
    "C1SC1[*]",
    "C(=O)OC(=O)[*]",
    "CC[*]",
    "N1CCN([*])CC1",
]

mid_fg = FGLibrary([DockMol(MolFromSmiles(i)) for i in nt_list], 'mid')

# terminal groups

t_list = [
    "Br",
    "Cl",
    "F",
    "C(=O)",
    "O",
    "CC#C",
    "C#N",
    "C(O)(=O)",
    "[N+]([O-])(=O)",
    "c1ccccc1",
    "N1C=CC=C1",
    "c1ccccn1",
    "C1CCCCC1",
    "C(N=[N+]=[N-])",
    "N",
    "C(Cl)(=O)",
    "C(=O)C",
    "S",
    "N(C)C",
    "N1C=NC=N1",
    "c1ccccc1(Cl)",
    "N1CCN(C)CC1",
    "N1CCOCC1",
    "c1cc(F)c(F)cc1(F)",
    "CO",
    "COC",
    "C1CC(O)C(CO)C1",
    "C1CC1"
]

terminal_fg = FGLibrary([DockMol(MolFromSmiles(i)) for i in t_list], 'terminal')
