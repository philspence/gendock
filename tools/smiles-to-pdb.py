import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

file = "smiles-strings.txt"

mols = []

with open(file) as f:
     for line in f:
        mols.append(line)#

nmol = 1
for m in mols:
    lname = "PDBmols/ligand_"+str(nmol)+".pdb"
    mol = AllChem.MolFromSmiles(m)
    mol = AllChem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    AllChem.MolToPDBFile(mol, lname)
    opt = 1
    while opt > 0:
       opt = AllChem.MMFFOptimizeMolecule(mol)
    print("Completed ligand_"+str(nmol)+".pdb")
    nmol += 1




    
