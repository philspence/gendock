import pybel as pb
import openbabel as ob
import csv


nmol = 1
results_array = []
for mol in pb.readfile("sdf", "NCI-Open_2012-05-01.sdf"):
    temp_array = []
    lname = "NCIFull/ligand_"+str(nmol)+".pdb"
    ligandname = "ligand_"+str(nmol)
    print("Writing "+str(ligandname))
    temp_array.append(ligandname)
    temp_array.append(mol.write())
    temp_array.append(mol.molwt)
    temp_array.append(mol.formula)
    temp_array.append(mol.data['NSC'])
    mol.addh()
    #mol.make3D(forcefield='uff')
    mol.localopt(forcefield='uff')
    writer = pb.Outputfile('pdb', lname, overwrite=True)
    writer.write(mol)
    print("Finished")
    nmol += 1
    results_array.append(temp_array)
csv_file = "ligand_info.csv"
headers = ["Ligand", "SMILES", "MolWt", "MolForm", "NSC"]
with open(csv_file, "w+") as f:
    writer = csv.writer(f, delimiter=",")
    writer.writerow(headers)
    writer.writerows(results_array)
print("Ligand information saved in ligand_info.csv")


    




    
