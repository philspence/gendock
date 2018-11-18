# MolDock

MolDock is a Python script that randomly generates molecules and then screens them against a receptor (e.g. a protein) using AutoDock Vina, to aid in finding new ligands for receptor binding. MolDock is ***not*** designed to give highly detailed computational analysis for the ligand-receptor binding, but ***is*** designed for very high throughput screening, the results of which can then be taken on for further analysis. MolDock can either generate molecules from scratch or it can add sidechains to an already established molecule.

## Setup

MolDock can be used on either Windows 10 (using Ubuntu WSL) or Linux (only Ubuntu 16.04 has been tested). MolDock will install all dependencies and set up the environment by running:

````
sudo bash install.sh
````

## Preparing Receptors

Save macromolecules/receptors in the **receptor** directory as **receptorX.pdb** where X is the numerical value, receptor1.pdb or receptor2.pdb, for example. The Autodock Vina configuration files are then named **rX-vina-config.txt**. These should be edited with your vina properties, grid size and location etc. I recommend using AutoDock Tools to find out the suitable grid size and location for your receptor. There is already a default config file named **r1-vina-config.txt** and the parameters are set as:
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
## Running MolDock

MolDock is designed to be very user configurable and runs in 4 stages:
1. Randomly generate molecules (either from nothing or based on a user-inputted molecule)
1. Optimise the geometry of the generated molecules
1. Prepare the molecules for docking as ligands in AutoDock Vina
1. Dock ligands using Autodock Vina

The following commands should be followed to run MolDock

```
source activate moldock
python3.5 moldock.py
```
Using the -h command will bring up the list of required arguments.

#### Experiment name (-x) *required*
This is the name of the experiment, all files will be saved the **data** directory. MolDock will create a directory with the name of the experiment that is given by this argument, e.g. **data/experiment1** if the argument **-x experiment1** is given.

#### Number of molecules to generate or dock (-n) *required*
This is the number of molecules to generate or dock. You can skip generating molecules if necessary (for example in the case where you already have a library of molecules), see the -g argument below.

#### Number of receptors (-r) *required*
MolDock will use as many receptors as you define, comparison of binding is easy using the CSV file generated. Make sure each receptor has its own vina configuration file. MolDock will then prepare the receptors for AutoDock Vina automatically before docking ligands.

#### Mass of molecules to generate (-m) *default = 400*
Use this argument to define the size of molecules that you wish to generate. The default is 400, this value is always an estimate, the mass will have a range of about 100, so for 400 you will most likely get a range of 350-450.

#### SMILES string of input molecule (-i)
MolDock can generate molecules from scratch or it can use an input molecule in the form of a SMILES string. An example is ``` Nc1nc(c2ccc(NC([*])=O)cn2)cs1 ``` where ``` [*] ``` denotes the position of functionalisation. ChemDraw (and other programs) have the ability to copy structures as SMILES strings (Ctr + Alt + C), set the functionalisation site to be labelled as **A** instead of **R** to have ```[*]``` in the SMILES string. Find an example molecule in the **input_mol** directory.

#### To Dock or not to Dock (-d) *default = 1*
Tell MolDock whether you want to dock (default) or to skip docking with (-d 0) if you just want to use MolDock to just generate a library of compounds.

#### Forcefield (-f) *default = MMFF94*
Select the forcefield for the geometry optimization from either the MMFF94 or UFF forcefields.

#### Generate molecules or use existing library (-g) *default = 1*
Use this to tell MolDock to generate its own molecules (default) or to skip this and use a preexisting library (-g 0). Move the PDB files of the library into the **data/exp_name/output_mols** directory before running MolDock.

## Results

MolDock will save a CSV file in the **data/exp_name/** directory that will contain the ligand number, SMILES string and the best binding energy for each receptor. MolDock stores the results of the AutoDock Vina screening in the **vina_files** directory. Inside each directory is both the PDBQT files and the LOG file. These are named as **ligand_X-r_Y.pdbqt or .txt** where X is the ligand number and Y is the receptor number.

## Editing
If you make any significant edits or add new functional groups, please contact me so I can try to incorporate these into the next version of MolDock.

## Extra Tools
MolDock includes some extra tools such as the **smiles-to-pdb.py** sctript that will convert the SMILES strings saved on each line of the **smiles-strings.txt** file to PDBs and store them for use in the **tools/PDBmols** directory. This is launched the same way as MolDock (within the Conda environment):
```
python3.5 smiles-to-pdb.py
```
MolDock can also convert SDF files containing multiple molecules to individual PDBs for screening (useful for something like the NCI library). To tell the script which SDF file to use, use the argument -f to define the filename:
```
python3.5 sdf-to-pdb.py -f filename.sdf
```
These PDBs will be saved in the **tools/PDBmols** directory.

## References

AutoDock Vina: *J Comput Chem.* **2010** Jan 30;31(2):455-61. [doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)

OpenBabel: *J. Cheminf.* **2011**, 3, 33. [doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33)

RDKit: Cheminformatics and Machine Learning Software. [RDKit.org](https://www.rdkit.org)
