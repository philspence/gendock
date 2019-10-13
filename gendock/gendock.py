import os
import csv
import rdkit
from rdkit import Chem
import gendock.scripts.molgen as mg
import gendock.scripts.docking as dk
import gendock as gd

def chkDirectory():
    homeDir = os.getcwd()
    if not os.path.exists('receptors'):
        print("No receptors directory found in the current working directory, quitting GenDock.")
        exit()
    else:
        return

def generate(name, target_mass, *args, **kwargs):
    data_path = os.path.join('data', name)
    if not os.path.exists(data_path):
        os.mkdir(data_path)
    nligands = kwargs.get('nligands', 0)
    mol = kwargs.get('mol', None)
    if mol is not None:
        if not type(mol) == rdkit.Chem.rdchem.Mol:
            print("mol is not an RDKit mol type.")
            return
    mg.molgenerate(name, target_mass, nligands, mol, gd.s_list, gd.nt_list, gd.t_list)

def dock(name, *args, **kwargs):
    ligand_num = kwargs.get('ligand_num', 1)
    receptors = kwargs.get('receptors', None)
    rnum = len(receptors)
    if rnum == 0:
        print('No receptor names given, quitting.')
        return
    # check if receptors exist
    for r in receptors:
        rpath = os.path.join('receptors', r + ".pdbqt")
        print(rpath)
        if not os.path.exists(rpath):
            return
    # prepare csv for results
    if ligand_num == 1:
        headers = ["Ligand", "Mw", "LogP", "SMILES", "NSC"]
        csv_fname = str(name) + "_results.csv"
        csv_file = os.path.join('data', name, csv_fname)
        for receptor in receptors:
            r_header = receptor + " Binding Affinity (kcal/mol)"
            headers.append(r_header)
        if not os.path.exists(os.path.join('data', name)):
            os.mkdir(os.path.join('data', name))
        with open(csv_file, "a") as f:
            writer = csv.writer(f, delimiter=",")
            writer.writerow(headers)
    # if not os.path.exists(receptor_pdbqt):
    #     print("Preparing receptors...")
    #     convertMol(receptor_pdb, receptor_pdbqt)
    # else:
    #     print("Receptor PDBQT found.")
    # convertMol(receptor_pdb, receptor_pdbqt)
    #dock
    dk.moldocking(name, ligand_num, receptors)


