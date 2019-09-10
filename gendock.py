import argparse
import os
import sys
import csv

from scripts.molgen import *
from scripts.docking import *

parser = argparse.ArgumentParser()
parser.add_argument('-x', dest='xp_name', metavar ='NAME', help='enter experiment number/name', required=True)
parser.add_argument('-n', dest='num_mols', metavar='X', help='number of molecules you want to generate, enter 0 if you want to skip generating', required=True, type=int)
parser.add_argument('-r', dest='num_recept', metavar='X', help='number of receptors to dock to, enter 0 if you want to skip docking', required=True, type=int)
parser.add_argument('-m', dest='target_mass', metavar='XXX', help='target mass of the generated molecules, default = 400', type=int, default='400')
parser.add_argument('-i', dest='input_smiles', metavar='SMILES String', help='SMILES string of starting molecule, default is to generate from scratch, see readme for more details', default=int('0'))
parser.add_argument('-l', dest='ligand_num', metavar='NUMBER', help='1 (default) will start from the first molecule in the library, enter another ligand number to start from there', type=int, default=int('1'))
parser.add_argument('--version', action='version', version='moldock_v1.00')
args = parser.parse_args()

# orig_stdout = sys.stdout
# log_file = os.path.join('data', args.xp_name, 'log.txt')
# f = open(log_file, 'w')
# sys.stdout = f

if args.ligand_num == 1:
    headers = ["Ligand", "Mw", "LogP", "SMILES", "NSC"]
    r_count = 1
    csv_fname = str(args.xp_name)+"_results.csv"
    csv_file = os.path.join('data', args.xp_name, csv_fname)
    while r_count <= args.num_recept:
        r_header = str("Receptor "+str(r_count)+" Binding Affinity (kcal/mol)")
        headers.append(r_header)
        r_count += 1
    with open(csv_file, "w+") as f:
            writer = csv.writer(f, delimiter=",")
            writer.writerow(headers)

platform = sys.platform
if platform == 'linux' or platform == 'linux2':
    mgltools = os.path.join('tools', 'mgltools_x86_64Linux2_1.5.6')
elif platform == 'darwin':
    mgltools = os.path.join('tools', 'mgltools_i86Darwin9_1.5.6')
else:
    print("You're not running a compatible OS")
    exit()

try:
    if not os.path.isdir(mgltools):
        print("mgltools cannot be found, fatal error")
        exit()
    else:
        pass
except:
    exit()

path_to_py_scripts = os.path.join(mgltools, 'MGLToolsPckgs', 'AutoDockTools', 'Utilities24')
pythonsh = os.path.join(mgltools, 'bin', 'pythonsh')

if args.num_mols > 0:
    molgenerate(args.xp_name, args.num_mols, args.target_mass, args.input_smiles)
else:
    print('Finding SDF file instead of generating molecules')

#prepare receptors
print("Preparing receptors...")
num = 1
while num <= args.num_recept:
    receptor_pdb = os.path.join('receptor', "receptor"+str(num)+".pdb")
    receptor_pdbqt = receptor_pdb.replace("pdb", "pdbqt")
    command = pythonsh+' '+path_to_py_scripts+'/prepare_receptor4.py'
    os.system(command+" -U nphs_lps_waters -r "+receptor_pdb+" -o "+receptor_pdbqt)
    num += 1
print("Finished preparing receptors.")
    
#dock and process or skip all together if the user just wants to generate molecules
if args.num_recept > 0:
    moldocking(args.xp_name, args.num_recept, path_to_py_scripts, args.ligand_num, pythonsh, args.num_mols)
else:
    print("Docking option was to skip...")
    print("Finished.")

# print("Saving log file")
# sys.stdout = orig_stdout
# f.close()


