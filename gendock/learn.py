from gendock.lipinski import check_lipinski
from keras.models import Model
from keras.layers import Embedding, Flatten, Dense, Input, Masking, Concatenate
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
import numpy as np
from pathlib import Path


class Predictor:
    def __init__(self):
        self.model = None
        self.data = None

    def create_model(self):
        residues_in = Input((2000,))
        res = Masking()
        res = Embedding(24, 4, input_length=2000)(residues_in)
        res = Flatten()(res)
        res = Dense(128, activation='relu')(res)  # essentially, this is "model1"
        res = Dense(128, activation='relu')(res)
        res = Dense(64, activation='relu')(res)
        # FPS layers
        fps_in = Input((2048,))
        fps = Dense(64, activation='linear')(fps_in)
        fps = Dense(64, activation='linear')(fps)  # this is essentially "model2"
        # concat layers
        concat_layer = Concatenate()([res, fps])  # concat last layers of model1 and model2
        # output
        output = Dense(1, activation='linear')(concat_layer)  # output of 6 categories, with softmax
        self.model = Model(inputs=[residues_in, fps_in], outputs=output)
        self.model.compile(loss='categorical_crossentropy', optimizer='adam', metrics=['categorical_accuracy'])

    def load_model(self, path):
        self.create_model()
        self.model.load_weights(path.as_posix())

    def predict(self, x):
        y = self.model.predict(x)
        return y

    def initial_train(self):
        x, y = self.data
        self.model.fit(x, y)
        self.model.save_weights(Path('model_weights.h5').as_posix())

    def update_model(self):
        return

    def shape_data(self, receptor, mols):
        fps = [GetMorganFingerprintAsBitVect(m, 2, nBits=2048) for m in mols]
        fps_arr = np.array(fps).reshape(-1, 2048)
        energies = np.array([m.docking_energy for m in mols]).reshape(-1, 1)
        residues = np.tile(receptor.sequence, len(fps_arr)).reshape(len(fps_arr), -1)
        self.data = ((residues, fps_arr), energies)


def assess_mols(mols):
    okay_mols = list()
    finished_mols = list()
    for m in mols:
        criteria, finished = check_lipinski(m)
        if criteria and not finished:
            okay_mols.append(m)
        if criteria and finished:
            finished_mols.append(m)
    better_mols = [m for m in okay_mols if m.docking_energy < m.previous_docking_energy]
    # n_top_mols = int(n_mols * 0.2 // 1)
    # mols_dict = dict(sorted(((mol, mol.docking_energy) for mol in mols), key=lambda x: x[1], reverse=True))
    # top_mols = [mol for mol, energy in mols_dict.items()][:n_top_mols]
    return better_mols, finished_mols
