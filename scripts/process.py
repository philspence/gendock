import os
import numpy as np
import csv
import openbabel as ob
import pybel as pb
from rdkit import Chem
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

def process_mols(xp_num, l_target, r_target):
    data_path = os.path.join('data', xp_num)
    ligand_num = 1
    results_array = []
    while ligand_num <= l_target:
        try:
            print("Processing Ligand "+str(ligand_num))
            recept_num = 1
            temp_array = []
            temp_array.append('ligand_'+str(ligand_num))
            pdb_file = os.path.join(str(data_path), 'output_mols', 'ligand_'+str(ligand_num)+'.pdb')
            temp_array.append(get_Mw(pdb_file))
            temp_array.append(get_LogP(pdb_file))
            temp_array.append(get_SMILES(pdb_file))
            while recept_num <= r_target:
                ligand_file = os.path.join('data', str(xp_num), 'vina_files', 'ligand_'+str(ligand_num)+'-r_'+str(recept_num)+'.pdbqt')
                temp_array.append(str(get_energy(ligand_file)))
                recept_num += 1
            results_array.append(temp_array)
            print(results_array)
            ligand_num += 1
        except:
            temp_array = []
            temp_array.append('ligand_'+str(ligand_num))
            temp_array.append('N/A')
            temp_array.append('N/A')
            temp_array.append('N/A')
            while recept_num <= r_target:
                temp_array.append('N/A')
                recept_num += 1
            results_array.append(temp_array)
            print(results_array)
            ligand_num += 1
            pass
    headers = ["Ligand", "Mw", "LogP", "SMILES"]
    r_count = 1
    while r_count <= r_target:
        r_header = str("Receptor "+str(r_count)+" Binding Affinity (kcal/mol)")
        headers.append(r_header)
        r_count += 1
    csv_fname = str(xp_num)+"_results.csv"
    csv_file = os.path.join('data', xp_num, csv_fname)
    with open(csv_file, "w+") as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(headers)
        writer.writerows(results_array)
    print("File saved as "+str(csv_file))


    
        
        
    
    
