from pathlib import Path
from rdkit.Chem.AllChem import EmbedMolecule, MolToPDBFile, SanitizeMol, MolFromPDBFile
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from openbabel import openbabel as ob
from subprocess import Popen, TimeoutExpired
from tensorflow.keras.preprocessing.sequence import pad_sequences


def os_command(command):
    c = Popen(command, shell=True)
    try:
        c.wait(120)
        return c
    except TimeoutExpired:
        print("Command timed out")
        c.kill()


def export_pdbqt(file_name):
    ob_con = ob.OBConversion()
    ob_con.SetInAndOutFormats('pdb', 'pdbqt')
    mol = ob.OBMol()
    ob_con.ReadFile(mol, Path(f'{file_name}.pdb').as_posix())
    ob_con.WriteFile(mol, Path(f'{file_name}.pdbqt').as_posix())
    return Path(f'{file_name}.pdbqt')


class DockMol:
    def __init__(self, mol, prev):
        self.previous_docking_energy = prev
        self.docking_energy = 0
        self.product = None
        self.rdmol = mol
        self.predicted_energy = 0

    def get_energy(self):
        line = Path('temp-out.pdbqt').open(mode='r').readlines()[1]
        result = float(line.split(':')[1].split()[0])
        self.docking_energy = result

    def join(self, next_group):
        rxn = ReactionFromSmarts('[*:1][#0:2].[#0:3][#0:4][*:5]>>[*:1]-[*:5].[#0:2][#0:3][#0:4]')
        product = rxn.RunReactants((self.rdmol, next_group.rdmol))[0][0]
        SanitizeMol(product)
        return DockMol(product, self.previous_docking_energy)

    def cap_with_h(self):
        rxn = ReactionFromSmarts('[*:1]-[#0:2].[H:3]>>[*:1]-[H:3].[#0:2]')
        product = rxn.RunReactant(self.rdmol, 0)[0][0]
        SanitizeMol(product)
        return product


class Receptor:
    def __init__(self, name):
        self.name = name
        self.dir = Path('receptors')
        self.path = Path(self.dir, f"{name}.pdbqt")
        self.sequence = None
        if not self.path.exists():
            print('Receptor not found, place .pdbqt file in the receptors directory.')
            exit()
        self.get_sequence()

    def dock_mol(self, mol):
        dock_mol = mol.cap_with_h()
        EmbedMolecule(dock_mol)
        MolToPDBFile(dock_mol, Path('temp.pdb').as_posix())
        pdbqt_file = export_pdbqt('temp')
        f_out_pdbqt = Path(pdbqt_file.parent, 'temp-out.pdbqt')
        f_out_log = Path(pdbqt_file.parent, 'temp-out.txt')
        vina_command = f"vina --config {self.dir}/{self.name}-config.txt --ligand {pdbqt_file} --out {f_out_pdbqt} " \
                       f"--log {f_out_log}"
        print(vina_command)
        c = os_command(vina_command)
        mol.get_energy()

    def create_config(self, center, size):
        x_center, y_center, z_center = center
        x_size, y_size, z_size = size
        config_txt = f"receptor = {self.path}\n\n" \
                     f"center_x = {x_center}\n" \
                     f"center_y = {y_center}\n" \
                     f"center_z = {z_center}\n\n" \
                     f"size_x = {x_size}\n" \
                     f"size_y = {y_size}\n" \
                     f"size_z = {z_size}"
        Path(self.dir, f"{self.name}-config.txt").open(mode='w').write(config_txt)

    def get_sequence(self):
        three_letter = ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO',
                        'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']
        one_letter = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                      'Y']
        three_to_one = dict(zip(three_letter, one_letter))
        char_to_int = dict((c, i) for i, c in enumerate(one_letter))
        lines = [l for l in self.path.open(mode='r').readlines() if l.startswith('ATOM')]
        res_num = list()
        for line in lines:
            rn = f"{three_to_one[line[17:20]]}_{line[23:26]}"
            if rn not in res_num:
                res_num.append(rn)
        residues_as_ints = list()
        for r in res_num:
            ri = r.split('_')[0]
            residues_as_ints.append(char_to_int[ri])
        self.sequence = pad_sequences([residues_as_ints], padding='post', value=23., dtype='float', maxlen=2000)[0]
