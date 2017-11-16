from yaml import load

from atb_outputs.formats import g96, pdb
from atb_outputs.mol_data import mol_data_from_mol_data_dict, MolData

if __name__ == '__main__':
    with open('data/21.yaml') as fh:
        mol_data = mol_data_from_mol_data_dict(load(fh))

    for united in [True, False]:
        print(g96(mol_data, united=united))

    m = MolData('HETATM    1  C0  UNL     1      -3.254   2.034   1.801  1.00  0.00           C')
    print(pdb(m))
