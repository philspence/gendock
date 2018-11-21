from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import csv

#def learn(input):
input = "NCIDiv4-fpocket_results.csv"


testfp = []
testfp.append(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles("c1ccccc1"), 2))
testmol = []
for fp in testfp:
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    testmol.append(arr)

header = []
rawdata = []
infile = open(input, "r")
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
print('Calculating prediction...')
nrg_pred = rf.predict(testmol)
print(nrg_pred)
print('Done.')