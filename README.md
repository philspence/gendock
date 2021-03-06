# GenDock

GenDock is a Python script that will generate molecules to dock against a receptor using a combination of 
autodock-vina and machine learning to improve the speed of ligand identification

## Setup

GenDock requires a few other modules. Create a virtualenv that contains `rdkit tensorflow numpy autodock-vina`

```
git clone https://github.com/philspence/gendock.git
cd gendock
mkdir receptors
cp path/to/receptor.pdbqt receptors/receptor.pdbqt
```

## Running GenDock

```
Arguments:
-r receptor name (no extension)
-c centre of the search grid that vina will use
-s size of the search grid that vina will use
```

`python gendock.py -r receptor -c 0 0 0 -s 20 20 20`

## Results

GenDock will save a CSV file in the **data** directory that will contain the SMILES string of the ligand, the name 
of the receptor, the binding energy calculated from vina and the predicted energy (0 for those that aren't predicted). 

## References

AutoDock Vina: *J Comput Chem.* **2010** Jan 30;31(2):455-61. [doi:10.1002/jcc.21334](https://doi.org/10.1002/jcc.21334)

OpenBabel: *J. Cheminf.* **2011**, 3, 33. [doi:10.1186/1758-2946-3-33](https://doi.org/10.1186/1758-2946-3-33)

RDKit: Cheminformatics and Machine Learning Software. [RDKit.org](https://www.rdkit.org)
