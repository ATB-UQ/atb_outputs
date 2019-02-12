

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
_chem_comp.three_letter_code                     {comp_id} 
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
ATOM_LINE_TEMPLATE = "{comp_id:<3} {name:<3} {name:<3} {type:<2} {charge:1} {pdbx_align:1} {aromatic:1} " \
                     "{terminal_atom} {stereo_config} {x_model:>6.3f} {y_model:>6.3f} {z_model:>6.3f} {x_ideal:>6.3f} " \
                     "{y_ideal:>6.3f} {z_ideal:>6.3f} {name:<3} {comp_id:<3} {index:<}\n" \

'''
NLE N   N   N 0 1 N N N 16.557 39.518 17.898 0.720  1.773  0.288  N   NLE 1  
NLE CA  CA  C 0 1 N N S 15.812 40.611 17.285 0.763  0.319  0.492  CA  NLE 2  
NLE C   C   C 0 1 N N N 16.773 41.690 16.789 2.084  -0.218 0.003  C   NLE 3  
NLE O   O   O 0 1 N N N 16.479 42.322 15.753 2.747  0.426  -0.776 O   NLE 4  
NLE OXT OXT O 0 1 N Y N 17.818 41.883 17.441 2.524  -1.411 0.433  OXT NLE 5  
NLE CB  CB  C 0 1 N N N 14.816 41.205 18.283 -0.375 -0.340 -0.289 CB  NLE 6  
NLE CG  CG  C 0 1 N N N 13.697 40.254 18.678 -1.718 0.110  0.290  CG  NLE 7  
NLE CD  CD  C 0 1 N N N 12.730 40.911 19.645 -2.857 -0.549 -0.491 CD  NLE 8  
NLE CE  CE  C 0 1 N N N 11.636 39.956 20.071 -4.200 -0.099 0.087  CE  NLE 9  
NLE H   1HN H 0 1 N N N 16.728 38.807 17.216 0.822  2.004  -0.689 H   NLE 10 
NLE HN2 2HN H 0 1 N Y N 17.429 39.863 18.245 -0.129 2.166  0.666  HN2 NLE 11 
NLE HA  HA  H 0 1 N N N 15.250 40.215 16.426 0.652  0.097  1.553  HA  NLE 12 
NLE HXT HXT H 0 1 N Y N 18.329 42.568 17.026 3.377  -1.713 0.092  HXT NLE 13 
NLE HB2 1HB H 0 1 N N N 15.369 41.477 19.194 -0.315 -0.046 -1.337 HB2 NLE 14 
NLE HB3 2HB H 0 1 N N N 14.345 42.069 17.792 -0.290 -1.424 -0.211 HB3 NLE 15 
NLE HG2 1HG H 0 1 N N N 13.147 39.956 17.773 -1.779 -0.184 1.338  HG2 NLE 16 
NLE HG3 2HG H 0 1 N N N 14.143 39.379 19.173 -1.803 1.194  0.211  HG3 NLE 17 
NLE HD2 1HD H 0 1 N N N 13.286 41.234 20.538 -2.796 -0.255 -1.539 HD2 NLE 18 
NLE HD3 2HD H 0 1 N N N 12.263 41.768 19.138 -2.772 -1.633 -0.413 HD3 NLE 19 
NLE HE1 1HE H 0 1 N N N 11.747 39.724 21.141 -4.284 0.985  0.009  HE1 NLE 20 
NLE HE2 2HE H 0 1 N N N 10.655 40.422 19.897 -5.011 -0.568 -0.469 HE2 NLE 21 
NLE HE3 3HE H 0 1 N N N 11.711 39.028 19.485 -4.260 -0.393 1.135  HE3 NLE 22 '''

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

BONDS_TEMPLATE = "{comp_id:<3} {name1:<3} {name2:<3} {bond_order:4} {aromatic:1} {stereo_config:1} {index:<}\n"
'''
NLE N   CA  SING N N 1  
NLE N   H   SING N N 2  
NLE N   HN2 SING N N 3  
NLE CA  C   SING N N 4  
NLE CA  CB  SING N N 5  
NLE CA  HA  SING N N 6  
NLE C   O   DOUB N N 7  
NLE C   OXT SING N N 8  
NLE OXT HXT SING N N 9  
NLE CB  CG  SING N N 10 
NLE CB  HB2 SING N N 11 
NLE CB  HB3 SING N N 12 
NLE CG  CD  SING N N 13 
NLE CG  HG2 SING N N 14 
NLE CG  HG3 SING N N 15 
NLE CD  CE  SING N N 16 
NLE CD  HD2 SING N N 17 
NLE CD  HD3 SING N N 18 
NLE CE  HE1 SING N N 19 
NLE CE  HE2 SING N N 20 
NLE CE  HE3 SING N N 21 '''