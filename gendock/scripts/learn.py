from rdkit import Chem, DataStructs
from rdkit.ML.Cluster import Butina
from rdkit.Chem import Draw, AllChem
import matplotlib.pyplot as plt
import numpy as np
import csv
import os

def ClusterFps(fps,cutoff):
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1,nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i],fps[:i])
        dists.extend([1-x for x in sims])
    # now cluster the data:
    cs = Butina.ClusterData(dists,nfps,cutoff,isDistData=True, reordering=True)
    return cs

def StandardDev(avg, values, length):
    sum_sqrs = 0
    for val in values:
        sqr = (val - avg) ** 2
        sum_sqrs += sqr
    mean_sqrs = sum_sqrs / length
    sd = mean_sqrs ** 0.5
    return sd


#def learn(input):
learn_xp = input("Enter data set to learn:")
input_learn_csv = os.path.join('data', learn_xp, learn_xp+'_results.csv')

infile = open(input_learn_csv, "r")
reader = csv.reader(infile, delimiter=',')
rownum = 0
mols = []
print('Extracting data...')
for row in reader:
    if rownum == 0:
        pass
    else:
        if row[5] == 'N/A':
            pass
        else:
            np_fp = np.array(0,)
            mol = Chem.MolFromSmiles(row[3])
            mol.SetProp('Energy R1', row[5])
            mol.SetProp('Energy R2', row[6])
            mol.SetProp('NSC', row[4])
            #set other props here for clustering e.g. TO displacement, etc.
            mols.append(mol)
    rownum += 1
infile.close()
print('Done.')

print('Creating clusters...')
fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024) for mol in mols]
clusters = np.array(ClusterFps(fps, 0.5))
print('Done.')

cnt = 0
avg_nrgs = []
sds = []
clustered = 0
clstd_mols = []
for cluster in clusters:
    ncluster = len(cluster)
    sum_nrgs = 0
    nrgs = []
    temp_mols = []
    for i in cluster:
        nrg = float(mols[i].GetProp('Energy R1')) #replace
        sum_nrgs += nrg
        nrgs.append(nrg)
        temp_mols.append(mols[i])
    avg = sum_nrgs / ncluster
    sd = StandardDev(avg, nrgs, ncluster)
    # print('Average for cluster ' + str(cnt) + ': ' + str(avg) + ' p/m ' + str(sd) + '. nmols = ' + str(ncluster))
    if sd > 0:
        clstd_mols.append(temp_mols)
        avg_nrgs.append(avg)
        sds.append(sd)
        clustered += ncluster
    cnt += 1
print('Total clusters: ' + str(cnt))
print('Total compounds in clusters > 1: ' + str(clustered))
print('Best cluster is ' + str(avg_nrgs.index(min(avg_nrgs))) + ' with energy of ' + str(min(avg_nrgs)))

best = clstd_mols[avg_nrgs.index(min(avg_nrgs))]
nbest = len(best)
img = Draw.MolsToGridImage(best[:nbest], molsPerRow=3, legends=[str(i.GetProp('NSC') + ', ' + i.GetProp('Energy R1')) for i in best]) #replace
img.save('best-mols-by-energy-R1.png') #replace

plt.figure()
x = range(0, len(avg_nrgs))
plt.xlabel('Cluster Number')
plt.ylabel('Average Binding Energy (kcal / mol) of pocket 1') #replace
plt.errorbar(x, avg_nrgs, yerr=sds, fmt='s', capsize=3)
plt.savefig('clusters-by-energy-R1.png') #replace
print('Figure saved.')

