import argparse
from gendock.functional_groups import start_fg, mid_fg
from gendock.docking import DockMol, Receptor
from gendock.learn import Predictor, assess_mols
from rdkit.Chem import MolToSmiles
from pathlib import Path


def main(m, r, c, s):
    model = Predictor()
    model.load_model(Path('model_weights.h5').as_posix())
    receptor = Receptor(r)
    receptor.create_config(center=c, size=s)
    complete = False
    completed_mols = list()
    if not m:
        mols = start_fg.return_all()
    else:
        mols = [DockMol(MolToSmiles(m))]
    while not complete:
        for m in mols:
            receptor.dock_mol(m)
        best_mols, finished_mols = assess_mols(mols)
        if not best_mols:
            complete = True
        next_mols = [m.join(i) for m in best_mols for i in mid_fg.return_all()]
        # predict values for new mols
        mols = next_mols.copy()
        completed_mols.append(finished_mols)
    # add some sort of data output for finished mols


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mol', type=str, required=False, help='SMILES string of mol')
    parser.add_argument('-r', '--receptor', type=str, required=True, help='name of receptor')
    parser.add_argument('-c', '--center', type=tuple, required=True, help='center of receptor grid: (x, y, z)')
    parser.add_argument('-s', '--size', type=tuple, required=True, help='size of receptor grid: (x, y, z)')
    args = parser.parse_args()

    if args.mol:
        mol = args.mol
    else:
        mol = None
    main(mol, args.receptor, args.center, args.size)
