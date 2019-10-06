import os
import csv

import scripts.molgen as gen
import scripts.docking as dk

def generate(name, target_mass, *args, **kwargs):
    nligands = kwargs.get('nligands', 0)
    gen.molgenerate(name, target_mass, nligands)

def dock(name, num_recept, *args, **kwargs):
    ligand_num = kwargs.get('ligand_num', 1)
    #prepare csv for results
    if ligand_num == 1:
        headers = ["Ligand", "Mw", "LogP", "SMILES", "NSC"]
        r_count = 1
        csv_fname = str(name) + "_results.csv"
        csv_file = os.path.join('data', name, csv_fname)
        while r_count <= num_recept:
            r_header = str("Receptor " + str(r_count) + " Binding Affinity (kcal/mol)")
            headers.append(r_header)
            r_count += 1
        if not os.path.exists(os.path.join('data', name)):
            os.mkdir(os.path.join('data', name))
        with open(csv_file, "a") as f:
            writer = csv.writer(f, delimiter=",")
            writer.writerow(headers)

    #prepare receptors
    num = 1
    while num <= num_recept:
        receptor = os.path.join('receptor', "receptor"+str(num)+".pdbqt")
        if not os.path.exists(receptor):
            print("Receptor PDBQT not found.")
            exit()
        else:
            print("Receptor PDBQT found")
        # if not os.path.exists(receptor_pdbqt):
        #     print("Preparing receptors...")
        #     convertMol(receptor_pdb, receptor_pdbqt)
        # else:
        #     print("Receptor PDBQT found.")
        # convertMol(receptor_pdb, receptor_pdbqt)
        num += 1
    #dock
    dk.moldocking(name, num_recept, ligand_num)


