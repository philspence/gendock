import os
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import random
import scripts.functional_groups as fgs

def findLargestRingNum(s):
    max_num = 0
    for x in range(len(s)):
        if s[x].isdigit() and int(s[x]) > max_num:
            max_num = int(s[x])
    return max_num

def replace_R(mol, fg):
    wc = Chem.MolFromSmiles("[*]")
    newmol = Chem.ReplaceSubstructs(mol, wc, fg, replaceAll=True)
    print(Chem.MolToSmiles(mol))
    Chem.SanitizeMol(newmol[0])
    return newmol[0]

def getAverages(list):
    tot = 0
    num = 0
    for i in list:
        tempmol = Chem.MolFromSmiles(i)
        tot += Descriptors.MolWt(tempmol)
        num += 1
    avg = tot / num
    return avg

def writeMol(mol, sdf_path):
    writer = Chem.SDWriter(sdf_path)
    writer.write(mol)

#MOLECULE GENERATION
def molgenerate(name, target_mass, nligands, mol):
    if mol is not None:
        avgSGp = Descriptors.MolWt(mol)
        smi = Chem.MolToSmiles(mol)
        fgs.s_list = [smi]
    else:
        avgSGp = getAverages(fgs.s_list)  # gets avg of the masses of each dictionary
    avgNTGp = getAverages(fgs.nt_list)
    avgTGp = getAverages(fgs.t_list)
    app_mass = target_mass - avgSGp - avgTGp  # subs avg mass of s and t groups from target
    NTGps_to_use = int(app_mass // avgNTGp)  # finds out how many nt groups, '//' = rounded down to nearest integer
    # if this rounds to 0  then no nt groups will be used just s and t groups!
    perm_max = [len(fgs.s_list)] + [len(fgs.nt_list)] * NTGps_to_use + [len(fgs.t_list)] #this is the max num of each group that will be gen'd
    total_cmpds = 1
    for i in perm_max:
        total_cmpds = total_cmpds * i #calculates the total permutations i.e. total cmpds to gen
    print("This will generate " + str(total_cmpds) + " compounds...")
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
    if nligands > 0:
        print('Generating ' + str(nligands) + ' compounds. This can take a while...')
    else:
        print("Generating all compounds. This can take a while...")
    if nligands > 0:
        temp_array = random.sample(perm_array, nligands)
        perm_array = temp_array
    #NEED TO ADD ABILITY TO PARSE MORE THAN ONE WILDCARD!
    ligand_num = 1
    for perm in perm_array: #for each permutation in the permutations array
        mol = Chem.MolFromSmiles(fgs.s_list[perm[0]]) #mol = the nth (the value of [x,0,0,0] in the array (perm[0]) item in the s_list
        nt_grp = 1 #set for the first nt_grp
        while nt_grp <= NTGps_to_use: #if ntgrp to use = 1 then will only pass once, if 0 it won't pass at all
            fg = Chem.MolFromSmiles(fgs.nt_list[perm[nt_grp]]) #get the correct group from the list
            mol = replace_R(mol, fg) #replace the position of the wildcard with the new group
            nt_grp += 1 #repeat
        fg = Chem.MolFromSmiles(fgs.t_list[perm[nt_grp]]) #same as above but the the t group i.e. 1 more than the last nt group
        mol = replace_R(mol, fg)
        # mol.SetProp("Ligand_Number", ligand_num)
        Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        sdf_file = os.path.join('data', name, name + '.sdf')
        sdf_open = open(sdf_file, 'a') #open in append mode, saves holding 1000s of compounds in memory
        writeMol(mol, sdf_open) #write them to the sdf
        print('Generated ligand ' + str(ligand_num) + '.')
        ligand_num += 1
    print("Generation complete. Molecules saved in " + sdf_file)
