from rdkit import Chem
from rdkit.Chem import ChemicalFeatures, AllChem
from rdkit import RDConfig
import os
import csv

#set definitions of features
fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
factory = ChemicalFeatures.BuildFeatureFactory(fdefName)
csv_file = 'features.csv'

#import sdf file and get features for each mol
mols = []
suppl = Chem.SDMolSupplier('NCIDiv6-all.sdf')
mols = []
ligand_cnt = 1


def calcDistance(m1, m2):
    dist = ((m1[0] - m2[0])**2 + (m1[1] - m2[1])**2 + (m1[2] - m2[2])**2) ** 0.5
    return dist

for mol in suppl:
    eachMol = []
    eachMol.append(Chem.MolToSmiles(mol)) #get smiles and NSC
    if mol.HasProp('E_NSC'):
        eachMol.append(mol.GetProp('E_NSC'))
    elif mol.HasProp('NSC'):
        eachMol.append(mol.GetProp('NSC'))
    else:
        eachMol.append('N/A')
    feats = factory.GetFeaturesForMol(mol)
    nfeats = len(feats)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    Hd_cnt = 0
    Ha_cnt = 0
    feat_cnt = 0
    d_pos = []
    a_pos = []
    while feat_cnt < nfeats:
        if feats[feat_cnt].GetFamily() == 'Donor':
            Hd_cnt += 1
            d_pos.append(feat_cnt)
        if feats[feat_cnt].GetFamily() == 'Acceptor':
            Ha_cnt += 1
            a_pos.append(feat_cnt)
        feat_cnt += 1
    eachMol.append(Hd_cnt)
    eachMol.append(Ha_cnt)
    if len(d_pos) > 2:
        dist = calcDistance(list(feats[d_pos[0]].GetPos()), list(feats[d_pos[1]].GetPos()))
        print(dist)
    mols.append(eachMol)
    print(ligand_cnt)
    ligand_cnt += 1

with open(csv_file, "w") as f:
    writer = csv.writer(f, delimiter=",")
    for row in mols:
        writer.writerow(row)