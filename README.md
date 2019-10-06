# GenDock

GenDock is a Python script that can either randomly generate molecules or use a preexisting library of molecules to screen against a macromolecule (or multiple macromolecules) using AutoDock Vina, to aid in finding new ligands for receptor binding. Machine learning alogrithms are also being designed to aid in the identification process.

GenDock is ***not*** designed to give highly detailed computational analysis for the ligand-receptor binding, but ***is*** designed for high throughput screening, the results of which can then be taken on for further analysis. 

## Setup

GenDock is only supported on macOS and Linux. Windows users can install Ubuntu using WSL and use GenDock from there.

I recommend using [Anaconda](https://www.anaconda.com/download) to setup your python environment. GenDock depends upon the installation of *autodock-vina*, *rdkit*, *openbabel*, and  *scikit-learn* - all of which can be handled by Anaconda and should installed to a single environment.

## Preparing Receptors

Save your macromolecules/receptors in the *receptor* directory as **receptorX.pdbqt** where X is the numerical value, e.g. **receptor1.pdbqt**. The Autodock Vina configuration files are then named **rX-vina-config.txt**. These should be edited with your vina properties, grid size and location etc. I recommend using AutoDock Tools to find out the suitable grid size and location for your receptor, as well as saving to a *.pdbqt* file format. There is already a default config file named **r1-vina-config.txt** and the parameters are set as:
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

The following commands should be followed to run GenDock:

```
import gendock as gd
gd.generate(name, target_mass, nligands=X)
gd.dock(name, num_recept, ligand_num=Y)
```
where:
 
**name** is the name of the experiment, all files will be saved the **data** directory. GenDock will create a directory with the name of the experiment that is given by this argument (if one does not already exist).

**target_mass** is the appproximate mass of the molecules that you want to generate.

**nligands** is the number of molecules you want to generate, if an argument is given, all possible molecules will be generated (this can range into the tens of millions of molecules). If a number is given, then that amount of molecules will be randomly chosen from the possible molecules that could be generated.

**num_recept** is the number of receptors you want to dock against.

**ligand_num** is the ligand that you want to start with. This is useful if you have generated 1000 ligands and then got cut off after docking 500 of them. Set this to 501 and it will carry on from where it left off.

## Results

GenDock will save a CSV file in the **data/exp_name/** directory that will contain the SMILES string of the ligand as well as other chemical properies, and the best binding energy for each receptor. GenDock stores the results of the AutoDock Vina screening in the **vina_files** directory. Inside each directory is both the PDBQT files and the LOG file. These are named as **ligand_X-r_Y.pdbqt or .txt** where X is the ligand number and Y is the receptor number.

## Functional Groups
Functional groups can be found in **scripts/functional_groups.py** as a list of SMILES strings. This file can be edited to add new functional groups.

## References

AutoDock Vina: *J Comput Chem.* **2010** Jan 30;31(2):455-61. [doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)

OpenBabel: *J. Cheminf.* **2011**, 3, 33. [doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33)

RDKit: Cheminformatics and Machine Learning Software. [RDKit.org](https://www.rdkit.org)
