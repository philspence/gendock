from rdkit.Chem import MolFromSmiles
from gendock.docking import DockMol
from random import choice
import json
from pathlib import Path


class FGLibrary:
    def __init__(self, fgs, position):
        self.position = position
        self.fgs = fgs
        return

    def get_random(self):
        return choice(self.fgs)

    def return_all(self):
        return self.fgs


def load_functional_groups():
    fgs = dict()
    library = json.load(Path('gendock', 'functional_groups.json').open(mode='r'))
    for lib, fg_list in library.items():
        fgs[lib] = FGLibrary([DockMol(MolFromSmiles(i), 0) for i in fg_list], lib)
    return fgs
