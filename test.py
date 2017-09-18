from yaml import load

from atb_outputsformats import g96
from atb_outputs.mol_data import mol_data_from_mol_data_dict

if __name__ == '__main__':
    with open('data/21.yaml') as fh:
        mol_data = mol_data_from_mol_data_dict(load(fh))

    for united in [True, False]:
        print(g96(mol_data, united=united))
