from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import csv

#def learn(input):
input_learn = "NCIDiv4-fpocket_results.csv"

header = []
infile = open(input_learn, "r")
reader = csv.reader(infile, delimiter=',')
rownum = 0
fps = []
nrgs = []
for row in reader:
    if rownum == 0:
        pass
    else:
        if 'N/A' in row:
            pass
        else:
            mol = Chem.MolFromSmiles(row[3])
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
            fps.append(fp)
            nrgs.append(float(row[4]))
    rownum += 1
infile.close()

np_fps = []
rf = RandomForestRegressor(n_estimators=1000, n_jobs=2)
for fp in fps:
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    np_fps.append(arr)
print('Fitting data...')
rf.fit(np_fps, nrgs)

print('Getting SDF mols')
# input_pred = Chem.SDMolSupplier('NCI-Open_2012-05-01.sdf')
input_pred = ['CC1=CC(=O)C=CC1=O', 'c1ccc2sc(SSc3nc4ccccc4s3)nc2c1', 'O=[N+]([O-])c1cc(Cl)c(O)c([N+](=O)[O-])c1', 'Nc1ncc([N+](=O)[O-])s1', 'Nc1ccc2c(c1)C(=O)c1ccccc1C2=O', 'O=C(O)c1ccccc1-c1c2ccc(=O)c(Br)c-2oc2c(Br)c(O)ccc12', 'CN(C)C1=C(Cl)C(=O)c2ccccc2C1=O', 'Cc1ccc2c(c1[N+](=O)[O-])C(=O)c1ccccc1C2=O', 'CC(=N\O)/C(C)=N/O', 'c1ccc(P(c2ccccc2)c2ccccc2)cc1']
header = []
fps_p = []
nrgs_p = []
maxnum = 10
for m in input_pred:
    mol = Chem.MolFromSmiles(m)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    fps_p.append(arr)

print('Done')
print('Calculating predictions...')
nrg_preds = []
for fp in fps_p:
    nrg_pred = rf.predict(fp)
    nrg_preds.append(nrg_pred)
print(nrg_preds)
print('Done.')