import os
import numpy as np
import csv
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate

def get_energy(file_name):
    file = open(file_name) 
    lines = file.readlines() 
    file.close() 
    line = lines[1] 
    result = float(line.split(':')[1].split()[0])  
    return result

def save_learn_data(mol, ligand_file, xp_num):
    AllChem.EmbedMolecule(mol)
    factory = Gobbi_Pharm2D.factory
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, factory, dMat=Chem.Get3DDistanceMatrix(mol))
    nrg = get_energy(ligand_file)
    row = [fp, nrg]
    learn_csv = str(xp_num)+"_learn_data.csv"
    with open(learn_csv, "a") as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(row)


def process(mol, r_num, l_num, xp_num):
    print("Trying to acquire ligand attributes and binding energy...")
    temp_array = []
    recept_num = 1
    temp_array.append('ligand_'+str(l_num))
    temp_array.append(Descriptors.MolWt(mol)
    temp_array.append(Descriptors.MolLogP(mol)))
    temp_array.append(Chem.MolToSmiles(mol))
    if mol.HasProp('NSC'):
        temp_array.append(mol.GetProp('NSC'))
    else:
        temp_array.append('N/A')
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
    print("Learning molecule data...")
    save_learn_data(mol, ligand_file, xp_num)

def makedir(dir):
    if not os.path.isdir(dir):
        os.makedirs(dir)

def dock(pdb_file, pdbqt_file, num_recept, vina_files_dir, molnum, path_to_py_scripts, pythonsh):
    prep_ligand = os.path.join(path_to_py_scripts, 'prepare_ligand4.py')
    command = pythonsh + " " + prep_ligand
    os.system(command + " -l " + pdb_file + " -o " + pdbqt_file)
    recept_num = 1
    while recept_num <= num_recept:
        # run vina
        print("Docking ligand " + str(molnum) + " with receptor " + str(recept_num) + "...")
        f_out_pdbqt = os.path.join(vina_files_dir, 'ligand_' + str(molnum) + '-r_' + str(recept_num) + '.pdbqt')
        f_out_log = os.path.join(vina_files_dir, 'ligand_' + str(molnum) + '-r_' + str(recept_num) + '.txt')
        try:
            vina_command = "vina --config receptor/r" + str(
                recept_num) + "-vina-config.txt --ligand " + pdbqt_file + " --out " + f_out_pdbqt + " --log " + f_out_log
            os.system(vina_command)
        except Exception:
            pass
        recept_num += 1

def moldocking(xp_num, num_recept, path_to_py_scripts, start_ligand, pythonsh, to_gen):
    data_dir = os.path.join('data', xp_num)
    vina_ligands_dir = os.path.join(data_dir, 'vina_ligands')
    makedir(vina_ligands_dir)
    vina_files_dir = os.path.join(data_dir, 'vina_files')
    makedir(vina_files_dir)
    if to_gen > 0: #run generated PDBs
        molnum = 1
        pdb_file = os.path.join(data_dir, 'output_mols', "ligand_" + str(molnum) + ".pdb")
        pdbqt_file = pdb_file.replace("pdb", "pdbqt")
        pdbqt_file = pdbqt_file.replace("output_mols", "vina_ligands")
        try:
            prep_ligand = os.path.join(path_to_py_scripts, 'prepare_ligand4.py')
            os.system("pythonsh "+prep_ligand+" -l "+str(pdb_file)+" -o "+str(pdbqt_file))
            print("Done.")
        except:
            print("Failed. Skipping ligand...")
            molnum += 1
        dock(pdb_file, pdbqt_file, num_recept, vina_files_dir, molnum, path_to_py_scripts, pythonsh)
        mol = AllChem.MolFromPDBFile(pdb_file)  #EDIT TO RDKIT
        process(mol, num_recept, molnum, xp_num)

    else: #run SDF file
        suppl = Chem.SDMolSupplier(os.path.join(data_dir, str(xp_num)+'.sdf'))
        pdb_file = os.path.join(data_dir, 'ligand.pdb')
        pdbqt_file = pdb_file.replace("pdb", "pdbqt")
        molnum = 1
        for mol in suppl:
            if molnum >= start_ligand:
                try:
                    AllChem.EmbedMolecule(mol)
                    AllChem.UFFOptimizeMolecule(mol, numThreads=0)
                    AllChem.MolToPDBFile(mol, pdb_file)
                    dock(pdb_file, pdbqt_file, num_recept, vina_files_dir, molnum, path_to_py_scripts, pythonsh)
                    os.remove(pdb_file)
                    os.remove(pdbqt_file)
                except Exception:
                    print('Failed to dock ligand '+str(molnum)+'. Skipping ligand...')
                process(mol, num_recept, molnum, xp_num)
            else:
                pass
            molnum += 1
