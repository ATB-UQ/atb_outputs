
CCD_DISCLAIMER = '''#
# Note:
#  - Stereo config flags for atoms and bonds have not been set.
#  - The double bond flag is currently only set for the following cases: C=C, C=N, N=N, and C=O;
#    the bond-length cutoffs used to assign double bonds (DOUB) were derived from those used by eLBOW (Phenix).
#
'''

MOLECULE_DESCRIPTERS_TEMPLATE = '''data_{comp_id}
# 
_chem_comp.id                                    {comp_id}
_chem_comp.name                                  ?
_chem_comp.type                                  ?
_chem_comp.pdbx_type                             ATOMP
_chem_comp.formula                               ?
_chem_comp.mon_nstd_parent_comp_id               ?
_chem_comp.pdbx_synonyms                         ?
_chem_comp.pdbx_formal_charge                    {net_charge}
_chem_comp.pdbx_initial_date                     ?
_chem_comp.pdbx_modified_date                    ?
_chem_comp.pdbx_ambiguous_flag                   ?
_chem_comp.pdbx_release_status                   ?
_chem_comp.pdbx_replaced_by                      ?
_chem_comp.pdbx_replaces                         ?
_chem_comp.formula_weight                        ?
_chem_comp.one_letter_code                       ?
_chem_comp.three_letter_code                     {comp_id_3char}
_chem_comp.pdbx_model_coordinates_details        ?
_chem_comp.pdbx_model_coordinates_missing_flag   N
_chem_comp.pdbx_ideal_coordinates_details        ATB
_chem_comp.pdbx_ideal_coordinates_missing_flag   N
_chem_comp.pdbx_model_coordinates_db_code        ?
_chem_comp.pdbx_subcomponent_list                ?
_chem_comp.pdbx_processing_site                  ?
'''

ATOMS_HEADER = '''#
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.alt_atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.pdbx_align
_chem_comp_atom.pdbx_aromatic_flag
_chem_comp_atom.pdbx_leaving_atom_flag
_chem_comp_atom.pdbx_stereo_config
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
_chem_comp_atom.pdbx_component_atom_id
_chem_comp_atom.pdbx_component_comp_id
_chem_comp_atom.pdbx_ordinal
'''
ATOM_LINE_TEMPLATE = "{comp_id:<4} {name:<3} {name:<3} {type:<2} {charge:1} {pdbx_align:1} {aromatic:1} " \
                     "{terminal_atom} {stereo_config} {x_model:>6.3f} {y_model:>6.3f} {z_model:>6.3f} {x_ideal:>6.3f} " \
                     "{y_ideal:>6.3f} {z_ideal:>6.3f} {name:<3} {comp_id:<4} {index:<}\n" \

BONDS_HEADER = '''#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
_chem_comp_bond.pdbx_aromatic_flag
_chem_comp_bond.pdbx_stereo_config
_chem_comp_bond.pdbx_ordinal
'''

BONDS_TEMPLATE = "{comp_id:<4} {name1:<3} {name2:<3} {bond_order:4} {aromatic:1} {stereo_config:1} {index:<}\n"
