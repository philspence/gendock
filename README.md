# GenDock

GenDock is a Python script that can either randomly generate molecules or use a preexisting library of molecules to screen against a macromolecule (or multiple macromolecules) using AutoDock Vina, to aid in finding new ligands for receptor binding. Machine learning alogrithms are also being designed to aid in the identification process.

GenDock is ***not*** designed to give highly detailed computational analysis for the ligand-receptor binding, but ***is*** designed for high throughput screening, the results of which can then be taken on for further analysis. 

## Setup

GenDock is only supported on macOS and Linux. Windows users can install Ubuntu using WSL and use GenDock from there.

I recommend using [Anaconda](https://www.anaconda.com/download) to setup your python environment. GenDock depends upon the installation of *autodock-vina*, *rdkit*, *openbabel*, and  *scikit-learn* - all of which can be handled by Anaconda and should installed to a single environment. You can then install gendock into this environment with:
```
pip install gendock
```

## Preparing Receptors

Save your macromolecules/receptors in the *receptor* directory as **[name].pdbqt** and the Autodock Vina configuration files must then be named **[name]-config.txt**. These should be edited with your vina properties, grid size and location etc. I recommend using AutoDock Tools to find out the suitable grid size and location for your receptor, as well as saving to a *.pdbqt* file format. There is already a default config file named **receptor1-config.txt** and the parameters are set as:
```
receptor = receptor/receptor.pdbqt

center_x =  0
center_y =  0
center_z =  0

size_x = 30
size_y = 30
size_z = 30

num_modes = 10
```
## Running GenDock

Ensure you run gendock from your working directory i.e. one that contains receptors in a directory named **receptors**. The following commands should be followed to run GenDock:

```
import gendock as gd
gd.generate(name, target_mass, nligands=X, mol=rdkit.Chem.rdchem.Mol)
gd.dock(name, ligand_num=Y, receptors=['receptor1', 'receptor2', 'receptor3'])
```
where:
 
**name** is the name of the experiment, all files will be saved the **data** directory. GenDock will create a directory with the name of the experiment that is given by this argument (if one does not already exist).

**target_mass** is the appproximate mass of the molecules that you want to generate.

**receptors** are the names of the receptor files, i.e. 1ELN (do not include the '.pdbqt' in this). This must be a list, even if you are only using one receptor.

***optional***, **nligands** is the number of molecules you want to generate, if an argument is given, all possible molecules will be generated (this can range into the tens of millions of molecules - *be warned*). If a number is given, then that amount of molecules will be randomly chosen from the possible molecules that could be generated.

***optional***, **mol** is an RDKit mol object that contains wildcard atoms i.e '`[*]`' atoms. For example one could generate an RDKit mol by the following:
```
from rdkit import Chem
import gendock as gd
m = Chem.MolFromSmiles('[*]c1cccc[*]c1')
gd.generate('test', 150, mol=m)
```


***optional***, **ligand_num** is the ligand that you want to start with. This is useful if you have generated 1000 ligands and then got cut off after docking 500 of them. Set this to 501 and it will carry on from where it left off.

## Results

GenDock will save a CSV file in the **data/exp_name/** directory that will contain the SMILES string of the ligand as well as other chemical properies, and the best binding energy for each receptor. GenDock stores the results of the AutoDock Vina screening in the **vina_files** directory. Inside each directory is both the PDBQT files and the LOG file. These are named as **ligand_X-r.pdbqt or .txt** where X is the ligand number and r is the receptor name.

## Functional Groups
Functional groups can be found in **scripts/functional_groups.py** as a list of SMILES strings. They are imported as **s_list**, **nt_list** and **t_list** and can be edited just like any another python list before running gendock.

```
import gendock as gd
gd.s_list.append('CCC[*]')
gd.s_list.remove('CCC[*]')
```

## References

AutoDock Vina: *J Comput Chem.* **2010** Jan 30;31(2):455-61. [doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)

OpenBabel: *J. Cheminf.* **2011**, 3, 33. [doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33)

RDKit: Cheminformatics and Machine Learning Software. [RDKit.org](https://www.rdkit.org)
