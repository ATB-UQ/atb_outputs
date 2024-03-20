from yaml import unsafe_load as load

from atb_outputs.formats import pdb, g96, mol_data_dict, lgf
from atb_outputs.itp import itp
from atb_outputs.mol_data import mol_data_from_mol_data_dict, MolData


def finalise_mol_data(m: MolData) -> MolData:
    m.var = {'REV_DATE': '', 'rnme': ''}
    m.completed = lambda x: False
    return m


if __name__ == '__main__':
    with open('data/21.yaml') as fh:
        mol_data = finalise_mol_data(mol_data_from_mol_data_dict(load(fh)))

    print(itp(mol_data))
    print(pdb(mol_data))
    # print(g96(mol_data))
    # print(mol_data_dict(mol_data))
    # print(lgf(mol_data))

