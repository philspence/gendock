import gc

#starting functional groups

class S_FuncGroup(object):
    def __init__ (self, smiles, mass):
        self.smiles = smiles
        self.mass = mass

phenyl12 = S_FuncGroup(smiles = "c1([*])c([**])cccc1", mass = 77)
phenyl13 = S_FuncGroup(smiles = "c1([*])cc([**])ccc1", mass = 77)
phenyl14 = S_FuncGroup(smiles = "c1([*])ccc([**])cc1", mass = 77)
phenyl134 = S_FuncGroup(smiles = "C1([*])=CC=C([**])C=C1([***])", mass = 76)
pyridine245 = S_FuncGroup(smiles = "C1([*])=CN=C([**])C=C1([***])", mass = 77)
pyrazine245 = S_FuncGroup(smiles = "C1([*])=CN=C([**])N=C1([***])", mass = 78)
#ethene = S_FuncGroup(smiles = "C([*])=C[**]", mass = 28)
#ethene_four = S_FuncGroup(smiles = "C([*])([**])=C([***])[****]", mass = 28)
quinoline28 = S_FuncGroup(smiles = "c1([*])cccc2c([**])ccnc21", mass = 127)
acridine9 = S_FuncGroup(smiles = "c1([*])c2c(cccc2)nc3ccccc31", mass = 128)
indan1 = S_FuncGroup(smiles = "C1([*])CCc2c1cccc2", mass = 117)
indan2 = S_FuncGroup(smiles = "C([*])(C1)Cc2c1cccc2", mass = 117)
indan13 = S_FuncGroup(smiles = "C([*])(CC1[**])c2c1cccc2", mass = 117)
naphpth1 = S_FuncGroup(smiles = "c1([*])c2c(cccc2)ccc1", mass = 127)
naphth12 = S_FuncGroup(smiles = "c1([*])c2c(cccc2)ccc1[**]", mass = 127)
naphth27 = S_FuncGroup(smiles = "c1([*])ccc2c(c1)cc([**])cc2", mass = 127)
naphth2 = S_FuncGroup(smiles = "c1([*])cc2c(cccc2)cc1", mass = 127)
mitox = S_FuncGroup(smiles = "N([*])C1=C(C(C2=C(C=C3)O)=O)C(C(C2=C3O)=O)=C(N[**])C=C1", mass = 268)
sildenafil = S_FuncGroup(smiles = "O=C1C2=C(C([*])=NN2[**])N=C([***])N1", mass = 136)
sild2 = S_FuncGroup(smiles = "C1([*])=NC(C([**])=NN2[***])=C2C=N1", mass = 120)
sild3 = S_FuncGroup(smiles = "C1([*])=NC(C([**])=NN2[***])=C2C=C1", mass = 119)
sild4 = S_FuncGroup(smiles = "C1([*])=CC(C([**])=NN2[***])=C2C=C1", mass = 118)
sild5 = S_FuncGroup(smiles = "C1([*])=NC(C([**])=CN2[***])=C2C=C1", mass = 118)
indol = S_FuncGroup(smiles = "C1([*])=NC(C([**])=CN2[***])=C2C=C1", mass = 120)

S_FuncGroup_list = []
for i in gc.get_objects():
    if isinstance(i, S_FuncGroup):
        S_FuncGroup_list.append(i)

#nonterminal functional groups

class NT_FuncGroup(object):
    def __init__ (self, smiles, mass):
        self.smiles = smiles
        self.mass = mass
        
phenyl12 = NT_FuncGroup(smiles = "c1c([*])cccc1", mass = 77)
phenyl13 = NT_FuncGroup(smiles = "c1cc([*])ccc1", mass = 77)
phenyl14 = NT_FuncGroup(smiles = "c1ccc([*])cc1", mass = 77)
alcohol = NT_FuncGroup(smiles = "OC([*])", mass = 30)
ester = NT_FuncGroup(smiles = "C(=O)OC([*])", mass = 44)
ether = NT_FuncGroup(smiles = "CC(=O)C([*])", mass = 57)
methylene = NT_FuncGroup(smiles = "C([*])", mass = 14)
sulphonic_ester = NT_FuncGroup(smiles = "CS(=O)(=O)O([*])", mass = 95)
ethene = NT_FuncGroup(smiles = "C/C=C/[*]", mass = 41)
ethyne = NT_FuncGroup(smiles = "CC#C[*]", mass = 39)
thioether = NT_FuncGroup(smiles = "CS([*])", mass = 47)
ketone = NT_FuncGroup(smiles = "CC(=O)C([*])", mass = 57)
amide = NT_FuncGroup(smiles = "CC(=O)N([*])", mass = 58)
sec_amine = NT_FuncGroup(smiles = "CN([*])", mass = 30)
sec_aldimine = NT_FuncGroup(smiles = "C=N[*]", mass = 42)
epoxide = NT_FuncGroup(smiles = "C1OC1[*]", mass = 44)
ethylene_sulphide = NT_FuncGroup(smiles = "C1SC1[*]", mass = 56)
azo = NT_FuncGroup(smiles = "N=N[*]", mass = 28)
imine = NT_FuncGroup(smiles = "C=N[*]", mass = 26)
acid_anhydride = NT_FuncGroup(smiles = "C(=O)OC(=O)[*]", mass = 72)
alkane_twoC = NT_FuncGroup(smiles = "CC[*]", mass = 28)
dinitroxane = NT_FuncGroup(smiles = "N1CCN([*])CC1", mass = 86)

NT_FuncGroup_list = []
for i in gc.get_objects():
    if isinstance(i, NT_FuncGroup):
        NT_FuncGroup_list.append(i)

#terminal functional groups
                
class T_FuncGroup(object):
    def __init__ (self, smiles, mass):
        self.smiles = smiles
        self.mass = mass

bromide = T_FuncGroup(smiles = "Br", mass = 80)
chloride = T_FuncGroup(smiles = "Cl", mass = 35.5)
fluoride = T_FuncGroup(smiles = "F", mass = 19)
methanal = T_FuncGroup(smiles = "C(=O)", mass = 29)
term_methanol = T_FuncGroup(smiles = "O", mass = 17)
ethyne = T_FuncGroup(smiles = "CC#C", mass = 40)
nitrile = T_FuncGroup(smiles = "C#N", mass = 26)
methanoic_acid = T_FuncGroup(smiles = "C(O)(=O)", mass = 45)
nitro = T_FuncGroup(smiles = "[N+]([O-])(=O)", mass = 46)
term_phenyl = T_FuncGroup(smiles = "c1ccccc1", mass = 77)
term_benzyl = T_FuncGroup(smiles = "c1ccccc1", mass = 91)
term_pyrrole = T_FuncGroup(smiles = "N1C=CC=C1", mass = 66)
term_pyridine2 = T_FuncGroup(smiles = "c1ccccn1", mass = 78)
term_cyclohex = T_FuncGroup(smiles = "C1CCCCC1", mass = 83)
term_azidide = T_FuncGroup(smiles = "C(N=[N+]=[N-])", mass = 42)
prim_amine = T_FuncGroup(smiles = "[N]", mass = 16)
acid_chloride = T_FuncGroup(smiles = "C(Cl)(=O)", mass = 62)
acetyl = T_FuncGroup(smiles = "C(=O)C", mass = 43)
thiol = T_FuncGroup(smiles = "S", mass = 32)
isocyanate = T_FuncGroup(smiles = "N=C=O", mass = 42)
thiazole = T_FuncGroup(smiles = "C1=NC=CS1", mass = 85)
triazole = T_FuncGroup(smiles = "N1C=NC=N1", mass = 69)
chlorobenz2 = T_FuncGroup(smiles = "c1ccccc1(Cl)", mass = 112)
methyl_dinitroxane14 = T_FuncGroup(smiles = "N1CCN(C)CC1", mass = 102)
morpholine = T_FuncGroup(smiles = "N1CCOCC1", mass = 88)
trifluorobenz = T_FuncGroup(smiles = "c1cc(F)c(F)cc1(F)", mass = 132)
methylhyroxy = T_FuncGroup(smiles = "CO", mass = 31)
methylether = T_FuncGroup(smiles = "COC", mass = 45)
ribose = T_FuncGroup(smiles = "C1CC(O)C(CO)C1", mass = 116)
cyclopropane = T_FuncGroup(smiles = "C1CC1", mass = 40)

T_FuncGroup_list = []
for i in gc.get_objects():
    if isinstance(i, T_FuncGroup):
        T_FuncGroup_list.append(i)
