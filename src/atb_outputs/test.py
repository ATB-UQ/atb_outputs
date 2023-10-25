from yaml import load

from atb_outputs.formats import pdb, g96, mol_data_dict, yml, pickle, template_yml, graph, lgf
from atb_outputs.mol_data import mol_data_from_mol_data_dict, MolData

def finalise_mol_data(m: MolData) -> MolData:
    m.var = {'REV_DATE': '', 'rnme': ''}
    m.completed = lambda x: False
    return m

if __name__ == '__main__':
    with open('data/21.yaml') as fh:
        mol_data = finalise_mol_data(mol_data_from_mol_data_dict(load(fh)))

    for united in [True, False]:
        print(pdb(mol_data, united=united))
        print(g96(mol_data, united=united))

    m = finalise_mol_data(MolData('HETATM    1  C0  UNL     1      -3.254   2.034   1.801  1.00  0.00           C'))

    ALL_OUTPUTS = [
        (pdb, {}),
        (g96, {'optimized': False}),
        #(mol_data_dict, {}),
        #(lgf, {}),
        #(graph, {}),
    ]

    for (output_function, output_kwargs) in ALL_OUTPUTS:
        print(output_function(m, **output_kwargs))
