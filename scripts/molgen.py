import random
import os
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from scripts.functional_groups import *

def findLargestRingNum(s):
    max_num = 0
    for x in range(len(s)):
        if s[x].isdigit() and int(s[x]) > max_num:
            max_num = int(s[x])
    return max_num

def choose_FuncGroup(fg_list):
    fg_nxt = random.choice(fg_list)
    return fg_nxt
    
def sort_RingNum(ng, mol_smiles):
    maxring_current = int(findLargestRingNum(mol_smiles))#finds largest num in current mol
    maxring_new = int(findLargestRingNum(ng))#finds largest num in new group, should always be 1
    if maxring_new >= 1: #basically, if it has a ring then...
        maxring_nxt = int(maxring_current)+1 #add 1 to current largest num in current mol
        if maxring_nxt >= 10: #if the new number is greater than 10...
            maxring_nxt = str("%"+str(maxring_nxt)) #change it from 10 to %10
        ng = ng.replace(str(maxring_new), str(maxring_nxt))#replace 1 in new group to current largest+1
    return ng

def replace_R(mol, fg, fg_mass, mol_m, cap):
    global mol_mass
    mol = mol.replace("[*]", fg)
    mol_mass = mol_m + fg_mass
    if "[**]" in mol:
        if cap == 0:
            fg = choose_FuncGroup(NT_FuncGroup_list)
        else:
            fg = choose_FuncGroup(T_FuncGroup_list)
        mol_smiles = fg.smiles
        mol_smiles = sort_RingNum(fg.smiles, mol)
        mol_smiles = mol_smiles.replace("[*]", "[**]")
        mol = mol.replace("[**]", mol_smiles)#replace incoming group [*] with [**]
        mol_mass += fg.mass
    if "[***]" in mol:
        if cap == 0:
            fg = choose_FuncGroup(NT_FuncGroup_list)
        else:
            fg = choose_FuncGroup(T_FuncGroup_list)
        mol_smiles = fg.smiles
        mol_smiles = sort_RingNum(fg.smiles, mol)
        mol_smiles = mol_smiles.replace("[*]", "[***]")
        mol = mol.replace("[***]", mol_smiles)#replace incoming group [*] with [***]
        mol_mass += fg.mass
    if "[****]" in mol:
        if cap == 0:
            fg = choose_FuncGroup(NT_FuncGroup_list)
        else:
            fg = choose_FuncGroup(T_FuncGroup_list)
        mol_smiles = fg.smiles
        mol_smiles = sort_RingNum(fg.smiles, mol)
        mol_smiles = mol_smiles.replace("[*]", "[****]")
        mol = mol.replace("[****]", mol_smiles) #replace incoming group [*] with [****]
        mol_mass += fg.mass
    return (mol, mol_mass)

#MOLECULE GENERATION
def molgenerate(xp_num, num_mols, target_mass, input_smiles, ff):
        mass_gap = 75
        preT_target_mass = target_mass - mass_gap
        ligand_num = 0
        to_gen = 1 #number of molecule that needs to be generated (obv start is 1)
        targets_gen = 0 #how many mols have been gen'd (obv 0 to start)
        while (targets_gen < num_mols):#keeps going till target is met
            print("Generating molecule "+str(to_gen)+"...")
            to_cap = 0
            if input_smiles == 0: #generate from scratch, default is 0
                fg_nxt = choose_FuncGroup(S_FuncGroup_list)
                mol_smiles = fg_nxt.smiles
                mol_mass = fg_nxt.mass
            else: #use the input smiles to start, taken from args entered into moldock.py
                mol_smiles = input_smiles
                mol_mass = Descriptors.MolWt(Chem.MolFromSmiles(mol_smiles))
            while mol_mass <= preT_target_mass:
                fg_nxt = choose_FuncGroup(NT_FuncGroup_list)
                fg_nxt.smiles = sort_RingNum(fg_nxt.smiles, mol_smiles) # fg_nxt.smiles = ... gets around local and global scopes
                (mol_smiles, mol_mass) = replace_R(mol_smiles, fg_nxt.smiles, fg_nxt.mass, mol_mass, to_cap)
            to_cap = 1
            fg_nxt = choose_FuncGroup(T_FuncGroup_list)
            fg_nxt.smiles = sort_RingNum(fg_nxt.smiles, mol_smiles)
            (mol_smiles, mol_mass) = replace_R(mol_smiles, fg_nxt.smiles, fg_nxt.mass, mol_mass, to_cap)          
            #increase the nums
            targets_gen += 1
            to_gen +=1
            #make filenames
            ligand_name = "ligand_"
            ligand_num += 1
            ligand_name_pdb = os.path.join('data', xp_num, 'output_mols', str(ligand_name+str(ligand_num)+".pdb")) #e.g. 'output_mols/ligand_1.pdb'
            ligand_name_sdf = ligand_name_pdb.replace('pdb', 'sdf')
            temp_pdb = ligand_name_pdb.replace(str(ligand_name+str(ligand_num)), "temp")
            
            #get mol from SMILES
            new_mol = Chem.MolFromSmiles(mol_smiles)
            #print(mol_smiles)
            #add Hs and embed before PDB
            print(mol_smiles)
            print("Adding hydrogens...")
            output_mol = AllChem.AddHs(new_mol)
            AllChem.EmbedMolecule(output_mol)
            print("Optimizing Geometry")
            opt = 1
            while opt > 0:
                if ff == "MMFF94":
                    opt = AllChem.MMFFOptimizeMolecule(output_mol)
                else:
                    opt = AllChem.UFFOptimizeMolecule(output_mol)
            print("Done.")
            #make dir if doesn't exist
            data_path = os.path.join('data', xp_num, 'output_mols')
            if not os.path.isdir(data_path):
                os.makedirs(data_path)
            print("Saving ligand_"+str(ligand_num))
            Chem.MolToPDBFile(output_mol, ligand_name_pdb)
            print("Done.")
                

	
