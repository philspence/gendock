from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from sklearn.ensemble import RandomForestRegressor
import numpy as np
import csv
import pickle
import os

#def learn(input):
learn_xp = input("Enter data set to learn:")
pred_xp = input("Enter data set to predict:")
input_learn_csv = os.path.join('data', learn_xp, learn_xp+'.csv')
pred_file = os.path.join('data', pred_xp, pred_xp+'.sdf')

filename = 'fitted_data.sav'

if not os.path.exists(filename): #looks to find if the data has been processed already
    infile = open(input_learn_csv, "r")
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
                AllChem.EmbedMolecule(mol)
                factory = Gobbi_Pharm2D.factory
                fp = Generate.Gen2DFingerprint(mol, factory, dMat=Chem.Get3DDistanceMatrix(mol))
                fps.append(fp)
                nrgs.append(float(row[4]))
        rownum += 1
    infile.close()
    np_fps = []
    rf = RandomForestRegressor(n_estimators=10000, n_jobs=4)
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)
    print('Fitting data...')
    rf.fit(np_fps, nrgs)
    pickle.dump(rf, open(filename, 'wb')) #creates file
else:
    print('Found fitted data file')
    rf = pickle.load(open(filename, 'rb')) #loads file

print('Getting SDF mols')
# input_pred = Chem.SDMolSupplier('NCI-Open_2012-05-01.sdf')

input_pred = []
suppl = Chem.SDMolSupplier(pred_file)
for mol in suppl:
    smi = Chem.MolToSmiles(mol)
    input_pred.append(smi)

fps_p = []
nrgs_p = []

for m in input_pred:
    mol = Chem.MolFromSmiles(m)
    factory = Gobbi_Pharm2D.factory
    fp = Generate.Gen2DFingerprint(mol, factory, dMat=Chem.Get3DDistanceMatrix(mol))
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    fps_p.append(arr)

print('Done')
print('Calculating predictions...')
nrg_preds = rf.predict(fps_p)

arraynum = 0
for m in input_pred:
    temparray = []
    temparray.append(m)
    temparray.append(nrg_preds[arraynum])
    csv_file = os.path.join('data', pred_xp, pred_xp+'_predict.csv')
    with open(csv_file, "a") as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(temparray)
    arraynum += 1
print('Done.')
