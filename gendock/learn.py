from gendock.lipinski import check_lipinski
from keras.models import Sequential
from keras.layers import Embedding, Flatten, Dense, Input, Masking, Concatenate
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
import numpy as np
from pathlib import Path


class Predictor:
    def __init__(self):
        self.model = None
        self.create_model()
        self.data = None
        self.path = Path('model_weights.h5')

    def create_model(self):
        model = Sequential()
        model.add(Dense(64, activation='linear'))
        model.add(Dense(64, activation='linear'))
        model.add(Dense(1, activation='linear'))
        model.compile(optimizer='adam', loss='mse', metrics='mae')
        self.model = model

    def load_model(self):
        self.model.load_weights(self.path.as_posix())

    def predict(self, mols):
        fps = [GetMorganFingerprintAsBitVect(m.rdmol, 2, nBits=2048) for m in mols]
        fps_arr = np.array(fps).reshape(-1, 2048)
        y = self.model.predict(fps_arr)
        for i in range(len(y)):
            mols[i].predicted_energy = y[i].item()
        return y

    def initial_train(self):
        x, y = self.data
        self.model.fit(x, y, epochs=10, verbose=1, batch_size=10, validation_split=0.2)
        self.model.save_weights(self.path.as_posix())

    def update_model(self):
        x, y = self.data
        self.load_model()
        self.model.fit(x, y, epochs=20, verbose=1, batch_size=10, validation_split=0.2)
        self.model.save_weights(self.path.as_posix())

    def shape_data(self, mols):
        fps = [GetMorganFingerprintAsBitVect(m.rdmol, 2, nBits=2048) for m in mols]
        fps_arr = np.array(fps).reshape(-1, 2048)
        energies = np.array([m.docking_energy for m in mols]).reshape(-1, 1)
        self.data = (fps_arr, energies)


def assess_mols(mols):
    finished_mols = list()
    for m in mols:
        criteria, finished = check_lipinski(m)
        if criteria and finished:
            finished_mols.append(m)
    # n_top_mols = int(n_mols * 0.2 // 1)
    # mols_dict = dict(sorted(((mol, mol.docking_energy) for mol in mols), key=lambda x: x[1], reverse=True))
    # top_mols = [mol for mol, energy in mols_dict.items()][:n_top_mols]
    return finished_mols
