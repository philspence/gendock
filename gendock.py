import argparse
from gendock.functional_groups import load_functional_groups
from gendock.docking import DockMol, Receptor
from gendock.learn import Predictor, assess_mols
from rdkit.Chem import MolToSmiles
from pathlib import Path
from rdkit.Chem import MolToSmiles


def main(r, c, s):
    model = Predictor()
    receptor = Receptor(r)
    receptor.create_config(center=c, size=s)
    func_groups = load_functional_groups()
    if not Path('data').exists():
        Path('data').mkdir()
    Path('data', f'{r}_data.csv').open(mode='w').write(f'SMILES, Receptor, Vina Energy, Predicted Energy\n')
    complete = False
    completed_mols = list()
    mols = func_groups['starter_fg'].return_all()
    counter = 0
    while not complete:
        if counter > 1:  # let it train for 2 rounds before predicting
            model.predict(mols)
            ordered = list(sorted(mols, key=lambda x: x.predicted_energy))[:100]
            mols = ordered
        for m in mols:
            receptor.dock_mol(m)
        finished_mols = assess_mols(mols)
        if len(finished_mols) == len(mols):
            complete = True
        # train model with mols
        model.shape_data(mols)
        if counter == 0:
            model.initial_train()
        else:
            model.update_model()
        next_mols = [m.join(i) for m in mols for i in func_groups['nonterminal_fg'].return_all()]
        mols = next_mols.copy()
        completed_mols.append(finished_mols)
        counter += 1
    # add some sort of data output for finished mols


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('-r', '--receptor', type=str, required=True, help='name of receptor')
    # parser.add_argument('-c', '--center', type=tuple, required=True, help='center of receptor grid: (x, y, z)')
    # parser.add_argument('-s', '--size', type=tuple, required=True, help='size of receptor grid: (x, y, z)')
    # args = parser.parse_args()
    # main(args.receptor, args.center, args.size)
    main('test', (30.515, 69.886, 8.891), (25, 25, 25))
