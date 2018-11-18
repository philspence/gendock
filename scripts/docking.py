import os
import numpy as np
import csv
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def get_energy(file_name):  
    file = open(file_name) 
    lines = file.readlines() 
    file.close() 
    line = lines[1] 
    result = float(line.split(':')[1].split()[0])  
    return result

def get_Mw(file_name):
    mol = Chem.MolFromPDBFile(file_name)
#   mol2 = AllChem.AddHs(mol)
    mw = Descriptors.MolWt(mol)    
    return mw

def get_LogP(file_name):
    mol = Chem.MolFromPDBFile(file_name)
    logp = Descriptors.MolLogP(mol)
    return logp

def get_SMILES(file_name):
    mol = Chem.MolFromPDBFile(file_name)
    smiles = Chem.MolToSmiles(mol)    
    return smiles

def process(molpdb, r_num, l_num, xp_num):
    print("Trying to acquire ligand attributes and binding energy...")
    try:
        temp_array = []
        recept_num = 1
        temp_array.append('ligand_'+str(l_num))
        temp_array.append(get_Mw(molpdb))
        temp_array.append(get_LogP(molpdb))
        temp_array.append(get_SMILES(molpdb))
        while recept_num <= r_num:
            ligand_file = os.path.join('data', str(xp_num), 'vina_files', 'ligand_'+str(l_num)+'-r_'+str(recept_num)+'.pdbqt')
            temp_array.append(str(get_energy(ligand_file)))
            recept_num += 1
        csv_fname = str(xp_num)+"_results.csv"
        csv_file = os.path.join('data', xp_num, csv_fname)
        with open(csv_file, "a") as f:
            writer = csv.writer(f, delimiter=",")
            writer.writerow(temp_array)
        print("Ligand attributes and binding energy acquired and appended to csv file")
    except:
        print("Failed to get ligand attributes")
        temp_array = []
        recept_num = 1
        temp_array.append('ligand_'+str(l_num))
        temp_array.append('N/A')
        temp_array.append('N/A')
        temp_array.append('N/A')
        while recept_num <= r_num:
            temp_array.append('N/A')
            recept_num += 1
        csv_fname = str(xp_num)+"_results.csv"
        csv_file = os.path.join('data', xp_num, csv_fname)
        with open(csv_file, "a") as f:
            writer = csv.writer(f, delimiter=",")
            writer.writerow(temp_array)

def moldocking(xp_num, num_mols, num_recept, path_to_py_scripts, start_ligand):
    data_path = os.path.join('data', xp_num)
    ligand_num = start_ligand
    while ligand_num <= num_mols:
        #set filenames
        print("Preparing to dock ligand_"+str(ligand_num)+"...")
        ligand_pdb = os.path.join(data_path, 'output_mols', "ligand_"+str(ligand_num)+".pdb")
        ligand_pdbqt = ligand_pdb.replace("pdb", "pdbqt")
        ligand_pdbqt = ligand_pdbqt.replace("output_mols", "vina_ligands")
        #make dir
        vina_path = os.path.join(data_path, 'vina_ligands')
        if not os.path.isdir(vina_path):
            os.makedirs(vina_path)
        #prepare ligand with openbabel
        print("Converting to pdbqt...")
        try:
            prep_ligand = os.path.join(path_to_py_scripts, 'prepare_ligand4.py')
            os.system("pythonsh "+prep_ligand+" -l "+str(ligand_pdb)+" -o "+str(ligand_pdbqt))
            print("Done.")
        except:
            print("Failed. Skipping ligand...")
            ligand_num += 1
        #make dirs
        vina_files_dir = os.path.join(data_path, 'vina_files')
        if not os.path.isdir(vina_files_dir):
            os.makedirs(vina_files_dir)
        recept_num = 1
        while recept_num <= num_recept:
            #run vina
            print("Docking ligand "+str(ligand_num)+" with receptor "+str(recept_num)+"...")
            f_out_pdbqt = os.path.join(str(vina_files_dir),'ligand_'+str(ligand_num)+"-r_"+str(recept_num)+".pdbqt")
            f_out_log = os.path.join(str(vina_files_dir), "ligand_"+str(ligand_num)+"-r_"+str(recept_num)+".txt")
            try:
                vina_command = "vina --config receptor/r"+str(recept_num)+"-vina-config.txt --ligand "+str(ligand_pdbqt)+" --out "+str(f_out_pdbqt)+" --log "+str(f_out_log)
            except:
                recept_num += 1
            os.system(vina_command)
            recept_num += 1
        process(ligand_pdb, num_recept, ligand_num, xp_num)
        ligand_num += 1
        
        
        
