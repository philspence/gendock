import os
from rdkit import Chem
from rdkit.Chem import AllChem

from scripts.functional_groups import *

def findLargestRingNum(s):
    max_num = 0
    for x in range(len(s)):
        if s[x].isdigit() and int(s[x]) > max_num:
            max_num = int(s[x])
    return max_num

def sort_RingNum(ng, mol_smiles):
    maxring_current = int(findLargestRingNum(mol_smiles)) # finds largest num in current mol
    maxring_new = int(findLargestRingNum(ng)) # finds largest num in new group, should always be 1
    if maxring_new >= 1: # basically, if it has a ring then...
        maxring_nxt = int(maxring_current)+1 # a dd 1 to current largest num in current mol
        if maxring_nxt >= 10: # if the new number is greater than 10...
            maxring_nxt = str("%"+str(maxring_nxt)) # change it from 10 to %10
        ng = ng.replace(str(maxring_new), str(maxring_nxt)) # replace 1 in new group to current largest+1
    return ng

def replace_R(mol, fg):
    mol = mol.replace("[*]", fg)
    if "[**]" in mol:
        mol_smiles = fg
        mol_smiles = sort_RingNum(fg, mol)
        mol_smiles = mol_smiles.replace("[*]", "[**]")
        mol = mol.replace("[**]", mol_smiles) # replace incoming group [*] with [**]
    if "[***]" in mol:
        mol_smiles = fg
        mol_smiles = sort_RingNum(fg, mol)
        mol_smiles = mol_smiles.replace("[*]", "[***]")
        mol = mol.replace("[***]", mol_smiles) # replace incoming group [*] with [***]
    if "[****]" in mol:
        mol_smiles = fg
        mol_smiles = sort_RingNum(fg, mol)
        mol_smiles = mol_smiles.replace("[*]", "[****]")
        mol = mol.replace("[****]", mol_smiles) # replace incoming group [*] with [****]
    return mol

def getAverages(dict):
    tot = 0
    num = 0
    for k, v in dict.items():
        tot += v
        num += 1
    avg = tot / num
    return avg

def createList(dict):
    list = []
    for i in dict:
        list.append(i)
    return list

def writeMol(mol, sdf_path):
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)

#MOLECULE GENERATION
def molgenerate(xp_num, to_gen, target_mass, input_smiles):
    avgSGp = getAverages(s_dict)  # gets avg of the masses of each dictionary
    avgNTGp = getAverages(nt_dict)
    avgTGp = getAverages(t_dict)
    app_mass = target_mass - avgSGp - avgTGp  # subs avg mass of s and t groups from target
    NTGps_to_use = int(app_mass // avgNTGp) # finds out how many nt groups, '//' = rounded down to nearest integer
    # if this rounds to 0  then no nt groups will be used just s and t groups!

    perm_max = [len(s_dict)] + [len(nt_dict)] * NTGps_to_use + [len(t_dict)] #this is the max num of each group that will be gen'd
    total_cmpds = 1
    for i in perm_max:
        total_cmpds = total_cmpds * i #calculates the total permutations i.e. total cmpds to gen
    print("This will generate " + str(total_cmpds) + " compounds...")
    print("Generating combinations...")
    perm_length = len(perm_max)
    pos_cntr = 0
    init_list = [0] * perm_length #make the first item i.e. [0,0,0,0] if using 4 groups
    perm_array = [init_list] #add it as the first list in the array
    while pos_cntr < perm_length: #pos_cntr will be the place edit i.e. [x,0,0,0] and length if the length of the array, i.e. 4
        cnt = 1 #has to be 1 not 0 for < perm_max[x] to work
        temp_array = perm_array.copy() #copy the array as it stands to a new one for it to be edited before being added back
        while cnt < perm_max[pos_cntr]: #counter has to be less than the max. available groups (i.e max. permutations)
            for list in temp_array: #for each list in the array change a value and add it to the permutations array
                temp_list = list.copy()
                temp_list[pos_cntr] = cnt
                perm_array.append(temp_list)
            cnt += 1
        pos_cntr += 1
    print('Done.')
    #convert the dicts to lists:
    s_list = createList(s_dict.keys())
    nt_list = createList(nt_dict.keys())
    t_list = createList(t_dict.keys())
    print("Generating compounds. This can take a while...")
    #NEED TO ADD ABILITY TO PARSE MORE THAN ONE WILDCARD!
    for perm in perm_array: #for each permutation in the permutations array
        mol_smiles = s_list[perm[0]] #mol_smiles = the nth (the value of [x,0,0,0] in the array (perm[0]) item in the s_list
        nt_grp = 1 #set for the first nt_grp
        while nt_grp <= NTGps_to_use: #if ntgrp to use = 1 then will only pass once, if 0 it won't pass at all
            fg = nt_list[perm[nt_grp]] #get the correct group from the list
            fg = sort_RingNum(fg, mol_smiles) #sort the ring numbers
            mol_smiles = replace_R(mol_smiles, fg) #replace the position of the wildcard with the new group
            nt_grp += 1 #repeat
        fg = t_list[perm[nt_grp]] #same as above but the the t group i.e. 1 more than the last nt group
        fg = sort_RingNum(fg, mol_smiles)
        mol_smiles = replace_R(mol_smiles, fg)
        mol = Chem.MolFromSmiles(mol_smiles)
        finalMol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(finalMol)
        AllChem.UFFOptimizeMolecule(finalMol, maxIters=500)
        sdf_path = open(os.path.join('data', xp_num, str(xp_num) + '.sdf'), 'a') #open in append mode, saves holding 1000s of compounds in memory
        writeMol(finalMol, sdf_path) #write them to the sdf
