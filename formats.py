from io import StringIO
from copy import deepcopy
import yaml
import pickle as pickle_module

from atb_outputs.helpers.types_helpers import MolData, Dict, Any, Output_File, Output_Files
import atb_outputs.pdb as PDB
import atb_outputs.yml as YML
import atb_outputs.lgf as LGF
import atb_outputs.graph as molecule_graph
import atb_outputs.lammps as LAMMPS


# This import creates an unfortunate circular dependency, but eventually
# the functions in algorithm.atb._outputs.output_helpers should be
# moved here to the atb_outputs module, which will make things neater
from algorithm.atb._outputs.output_helpers import atb_header


def pdb(mol_data: MolData, optimized: bool = True, united: bool = False, use_rnme: bool = True) -> Output_File:
    '''return a new pdb string reflecting changes of atom order and numbering'''
    io = StringIO()

    atoms = list(
        sorted(
            mol_data.atoms.values(),
            key=lambda x:x['index'],
        ),
    )

    PDB.header(mol_data, io, rev_date=mol_data.var["REV_DATE"], united=united)
    PDB.atoms(mol_data, io, atoms, united=united, use_rnme=use_rnme, optimized=optimized)
    PDB.connectivity(mol_data, io, atoms, united=united)
    PDB.footer(mol_data, io)

    return io.getvalue()

def mol_data_dict(mol_data: MolData) -> Dict[str, Any]:
    return {
        'atoms': mol_data.atoms,
         'bonds': YML.clean_bonds(mol_data.bonds),
         'angles': YML.clean_angles(mol_data.angles),
         'dihedrals': YML.clean_dihedrals(mol_data.dihedrals),
         'impropers': mol_data.impropers,
         'rings': mol_data.rings,
         'var': mol_data.var,
    }

def yml(mol_data: MolData) -> Output_File:
    mol_data = mol_data_dict(mol_data)
    return YML.add_yml_comments(yaml.dump(mol_data))

def pickle(mol_data: MolData) -> Output_File:
    return pickle_module.dumps(mol_data_dict(mol_data))

def template_yml(mol_data: MolData) -> Output_File:
    mol_data = deepcopy(mol_data)
    mol_data = {
        'atoms': YML.clean_atoms(mol_data.atoms, template=True),
         'bonds': YML.clean_bonds(mol_data.bonds, template=True),
         'angles': YML.clean_angles(mol_data.angles, template=True),
         'dihedrals': YML.clean_dihedrals(mol_data.dihedrals, template=True),
         'impropers': YML.clean_impropers(mol_data.impropers, template=True),
         'rings': YML.clean_rings(mol_data.rings, template=True),
         'var': mol_data.var,
    }
    return YML.add_yml_comments(yaml.dump(mol_data))

def graph(mol_data: MolData, **kwargs: Dict[str, Any]) -> Output_Files:
    return list(filter(
        lambda k_v: bool(k_v[1]),
        [
            (
                'svg',
                molecule_graph.graph_img(
                    mol_data,
                    **kwargs,
                ),
            ),
            (
                'gt',
                molecule_graph.graph_gt(
                    mol_data,
                    **kwargs,
                ),
            ),
        ],
    ))

def lgf(mol_data) -> Output_Files:
    try:
        return [
            ('lgf', LGF.graph(mol_data)),
        ]
    except AssertionError:
        return []

def lammps(mol_data: MolData, optimized: bool, united: bool) -> Output_Files:
    header_buffer = StringIO() # IO buffer stores the output string

    # Generate the ATB header
    atb_header(mol_data,
               header_buffer,
               LAMMPS.filename(mol_data.var["rnme"], united, optimized),
               united,
               comment_delimiter = "#")
    header = header_buffer.getvalue()

    return [(
        'lt',
        LAMMPS.moltemplate(mol_data_dict(mol_data), optimized, united, header)
            )]

