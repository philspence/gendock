from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
import numpy as np
import csv
import pandas as pd

#def learn(input):
input = "NCIDiv4-fpocket_results.csv"

header = []
rawdata = []
infile = open(input, "r")
reader = csv.reader(infile, delimiter=',')
rownum = 0
for row in reader:
    temp = []
    if rownum == 0:
        pass
    else:
        if 'N/A' in row:
            pass
        else:
            mol = Chem.MolFromSmiles(row[3])
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
            #temp.append(mol)
            temp.append(fp)
            temp.append(float(row[4]))
            rawdata.append(temp)
    rownum += 1
infile.close()


# data = pd.DataFrame({
#     'fingerprint':rawdata[:,0],
#     'energy':rawdata[:,1]
# })
# data.head()