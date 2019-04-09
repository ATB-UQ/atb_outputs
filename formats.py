from io import StringIO
from copy import deepcopy
import yaml
import pickle as pickle_module

from atb_outputs.helpers.types_helpers import MolData, Dict, Any, Output_File, Output_Files
import atb_outputs.pdb as PDB
import atb_outputs.yml as YML
import atb_outputs.lgf as LGF
import atb_outputs.ccd_cif as CIF
import atb_outputs.graph as molecule_graph

DOUBLE_BOND_LENGTH_CUTOFF = {
    frozenset(['C', 'C']): 0.139, #nm, Source: phenix.elbow.elbow.quantum.hf_631Gdp.py
    frozenset(['C', 'N']): 0.132, #nm, Source: phenix.elbow.elbow.quantum.hf_631Gdp.py
    frozenset(['N', 'N']): 0.136, #nm, Source: phenix.elbow.elbow.quantum.hf_631Gdp.py
    frozenset(['C', 'O']): 0.134, #nm, Source: phenix.elbow.elbow.quantum.hf_631Gdp.py
}


def ccd_cif(mol_data, comp_id):
    cif_str = CIF.CCD_DISCLAIMER
    cif_str += CIF.MOLECULE_DESCRIPTERS_TEMPLATE.format(
        comp_id=comp_id,
        net_charge=mol_data.var["total_charge"],
    )
    cif_str += CIF.ATOMS_HEADER
    aromatic_ring_atoms = [r["atoms"] for r in list(mol_data.rings.values()) if "aromatic" in r and r["aromatic"]]
    aromatic_atom_ids = [a for r in aromatic_ring_atoms for a in r]
    ocoord_key = "ocoord" if "ocoord" in list(mol_data.atoms.values())[0] else "coord"
    for atom in sorted(mol_data.atoms.values(), key=lambda x:x["index"]):
        line_params = dict(
            comp_id=comp_id,
            name=atom["symbol"],
            type=atom["type"],
            charge=0,
            pdbx_align=1,
            aromatic="Y" if atom["id"] in aromatic_atom_ids else "N",
            terminal_atom="N",
            stereo_config="N",
            x_model=atom["coord"][0]*10,
            y_model=atom["coord"][1]*10,
            z_model=atom["coord"][2]*10,
            x_ideal=atom[ocoord_key][0]*10,
            y_ideal=atom[ocoord_key][1]*10,
            z_ideal=atom[ocoord_key][2]*10,
            index=atom["index"],
        )

        cif_str += CIF.ATOM_LINE_TEMPLATE.format(
            **line_params
        )
    cif_str += CIF.BONDS_HEADER
    for i, bond in enumerate(sorted(mol_data.bonds, key=lambda x:mol_data.atoms[x['atoms'][0]]["index"])):
        bond_length = bond["value"]
        atom1_type = mol_data.atoms[bond["atoms"][0]]["type"]
        atom2_type = mol_data.atoms[bond["atoms"][1]]["type"]
        bond_order_id = frozenset([atom1_type, atom2_type])
        if bond_order_id in DOUBLE_BOND_LENGTH_CUTOFF:
            bond_order = "DOUB" if bond_length < DOUBLE_BOND_LENGTH_CUTOFF[bond_order_id] else "SING"
        else:
            bond_order = "SING"
        cif_str += CIF.BONDS_TEMPLATE.format(
            comp_id=comp_id,
            name1=mol_data.atoms[bond["atoms"][0]]["symbol"],
            name2=mol_data.atoms[bond["atoms"][1]]["symbol"],
            bond_order=bond_order,
            aromatic="Y" if all([atom_id in aromatic_atom_ids for atom_id in bond["atoms"]]) else "N",
            stereo_config="N",
            index=i+1,
        )
    return cif_str


def pdb(mol_data: MolData, optimized: bool = True, united: bool = False, use_rnme: bool = True) -> Output_File:
    '''return a new pdb string reflecting changes of atom order and numbering'''
    io = StringIO()

    atoms = list(
        sorted(
            filter(
                lambda atom: (not united) or (united and 'uindex' in atom),
                mol_data.atoms.values(),
            ),
            key=lambda atom: atom['uindex' if united else 'index'],
        ),
    )

    PDB.header(mol_data, io, rev_date=mol_data.var["REV_DATE"], united=united)
    PDB.atoms(mol_data, io, atoms, united=united, use_rnme=use_rnme, optimized=optimized)
    PDB.connectivity(mol_data, io, atoms, united=united)
    PDB.footer(mol_data, io)

    return io.getvalue()


def g96(mol_data: MolData, optimized: bool = True, united: bool = False) -> Output_File:
    io = StringIO()
    print_to_io = lambda *args: print(*args, file=io)

    atoms = list(
        sorted(
            filter(
                lambda atom: (not united) or (united and 'uindex' in atom),
                mol_data.atoms.values(),
            ),
            key=lambda atom: atom['uindex' if united else 'index'],
        ),
    )

    print_to_io('TITLE')
    print_to_io('')
    print_to_io('END')

    print_to_io('POSITION')
    for atom in atoms:
        # Source: GROMOS96 Manual (ISBN 3 7281 2422 2), page III-41
        print_to_io(
            '{residue_index:>5d}{X:1s}{residue_name:>5s}{X:1s}{atom_name:>5s}{atom_index:>7d}{x:15.9f}{y:15.9f}{z:15.9f}'.format(
                X=' ',
                residue_index=1,
                residue_name=mol_data.var['rnme'],
                atom_name=atom['symbol'],
                atom_index=atom['uindex' if united else 'index'],
                **dict(zip(('x', 'y', 'z'), atom['ocoord' if optimized else 'coord'])),
            ),
        )

    print_to_io('END')

    return io.getvalue()


def mol_data_dict(mol_data: MolData) -> Dict[str, Any]:
    def clean_up_flavours(atom: Dict[str, Any]) -> Dict[str, Any]:
        return {k: v for (k, v) in atom.items() if k != 'flavour'}

    return {
        'atoms': {atom_id: clean_up_flavours(atom) for (atom_id, atom) in mol_data.atoms.items()},
         'bonds': YML.clean_bonds(mol_data.bonds),
         'angles': YML.clean_angles(mol_data.angles),
         'dihedrals': YML.clean_dihedrals(mol_data.dihedrals),
         '_dihedrals': mol_data._dihedrals if hasattr(mol_data, '_dihedral') else None,
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


STORE_GRAPH_GT = False


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
        ]
        +
        (
            [
                (
                    'gt',
                    molecule_graph.graph_gt(
                        mol_data,
                        **kwargs,
                    ),
                ),
            ]
            if STORE_GRAPH_GT
            else
            []
        )
    ))


def lgf(mol_data, **kwargs: Dict[str, Any]) -> Output_Files:
    try:
        return [
            ('lgf', LGF.graph(mol_data, **kwargs)),
        ]
    except AssertionError:
        return []
