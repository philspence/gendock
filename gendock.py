import os
import csv
import rdkit
from rdkit import Chem
import scripts.molgen as gen
import scripts.docking as dk

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
    gen.molgenerate(name, target_mass, nligands, mol)

def dock(name, *args, **kwargs):
    ligand_num = kwargs.get('ligand_num', 1)
    receptor1 = kwargs.get('r1', None)
    receptor2 = kwargs.get('r2', None)
    receptor3 = kwargs.get('r3', None)
    receptors = [receptor1, receptor2, receptor3]
    print(receptors)
    rnum = 0
    for r in receptors:
        if r is not None:
            rnum += 1
    if rnum == 0:
        print('No receptor names given, quitting.')
        return
    receptors = receptors[:rnum]
    #prepare csv for results
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

    #prepare receptors
    rpaths = []
    for r in receptors:
        rpath = os.path.join('receptor', r + ".pdbqt")
        if not os.path.exists(rpath):
            print(rpath + " not found.")
            return
        else:
            rpaths.append(rpath)
    # if not os.path.exists(receptor_pdbqt):
    #     print("Preparing receptors...")
    #     convertMol(receptor_pdb, receptor_pdbqt)
    # else:
    #     print("Receptor PDBQT found.")
    # convertMol(receptor_pdb, receptor_pdbqt)
    #dock
    dk.moldocking(name, ligand_num, receptors)


