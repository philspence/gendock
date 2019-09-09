import random
import os
import pybel as pb

from scripts.functional_groups import *

#todo : change over from pybel to rdkit

def findLargestRingNum(s):
    max_num = 0
    for x in range(len(s)):
        if s[x].isdigit() and int(s[x]) > max_num:
            max_num = int(s[x])
    return max_num

def choose_FuncGroup(dict):
    rand = random.choice(list(dict.items()))
    return rand #rand[0] = smiles, rand[1] = mass
    
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
            fg = choose_FuncGroup(nt_dict)
        else:
            fg = choose_FuncGroup(t_dict)
        mol_smiles = fg[0]
        mol_smiles = sort_RingNum(fg[0], mol)
        mol_smiles = mol_smiles.replace("[*]", "[**]")
        mol = mol.replace("[**]", mol_smiles)#replace incoming group [*] with [**]
        mol_mass += int(fg[1])
    if "[***]" in mol:
        if cap == 0:
            fg = choose_FuncGroup(nt_dict)
        else:
            fg = choose_FuncGroup(t_dict)
        mol_smiles = fg[0]
        mol_smiles = sort_RingNum(fg[0], mol)
        mol_smiles = mol_smiles.replace("[*]", "[***]")
        mol = mol.replace("[***]", mol_smiles)#replace incoming group [*] with [***]
        mol_mass += int(fg[1])
    if "[****]" in mol:
        if cap == 0:
            fg = choose_FuncGroup(nt_dict)
        else:
            fg = choose_FuncGroup(t_dict)
        mol_smiles = fg[0]
        mol_smiles = sort_RingNum(fg[0], mol)
        mol_smiles = mol_smiles.replace("[*]", "[****]")
        mol = mol.replace("[****]", mol_smiles) #replace incoming group [*] with [****]
        mol_mass += int(fg[1])
    return (mol, mol_mass)

#MOLECULE GENERATION
def molgenerate(xp_num, num_mols, target_mass, input_smiles):
        mols = []
        mass_gap = 75
        preT_target_mass = target_mass - mass_gap
        ligand_num = 0
        to_gen = 1 #number of molecule that needs to be generated (obv start is 1)
        targets_gen = 0 #how many mols have been gen'd (obv 0 to start)
        while (targets_gen < num_mols):#keeps going till target is met
            print("Generating molecule "+str(to_gen)+"...")
            to_cap = 0
            if input_smiles == 0: #generate from scratch, default is 0
                fg_nxt = choose_FuncGroup(s_dict)
                mol_smiles = fg_nxt[0]
                mol_mass = int(fg_nxt[1])
            else: #use the input smiles to start, taken from args entered into moldock.py
                mol_smiles = input_smiles
                tempmol = pb.readstring("smi", mol_smiles)
                mol_mass = tempmol.molwt
            while mol_mass <= preT_target_mass:
                fg_nxt = choose_FuncGroup(nt_dict)
                fg_nxt[0] = sort_RingNum(fg_nxt[0], mol_smiles) # fg_nxt[0] = ... gets around local and global scopes
                (mol_smiles, mol_mass) = replace_R(mol_smiles, fg_nxt[0], int(fg_nxt[1]), mol_mass, to_cap)
            to_cap = 1
            fg_nxt = choose_FuncGroup(t_dict)
            fg_nxt[0] = sort_RingNum(fg_nxt[0], mol_smiles)
            (mol_smiles, mol_mass) = replace_R(mol_smiles, fg_nxt[0], int(fg_nxt[1]), mol_mass, to_cap)          
            #increase the nums
            targets_gen += 1
            to_gen +=1
            #get mol from SMILES
            new_mol = pb.readstring("smi", mol_smiles)
            print("Optimizing Geometry")
            newmol.make3D()
            mols.append(newmol)
        sdf_path = os.path.join('data', xp_num, str(xp_num)+'.sdf')
        sdf_file = pb.Outputfile("sdf", sdf_path)
        if not os.path.exists(os.path.join('data', xp_num)):
            os.makedirs(os.path.join('data', xp_num))
        for m in mols:
            sdf_file.write(m)
	
