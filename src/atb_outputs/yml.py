from copy import deepcopy
import re

def clean_atoms(atoms, template=False):
    if not template:
        return deepcopy(atoms)
    atoms_= {}
    for atom_id, atom in list(atoms.items()):
        atoms_[atom_id] = _make_atom_template(atom)
    return atoms_

def _make_atom_template(a):
    a["ocoord"] = []
    a["charge"] = None
    a["ucharge"] = None
    a["cgroup"] = None
    return a

def clean_bonds(bonds, template=False):
    bonds_ = []
    for b in bonds:
        if template:
            _make_bond_template(b)
        else:
            b["code"] = deepcopy(b["code"])
            b["value"] = float(b["value"])

        bonds_.append(b)
    return bonds_

def _make_bond_template(b):
    b["code"] = []
    b['fc'] = None
    if "hfc" in b:
        del b["hfc"]
    b["value"] = None
    b["order_qm"] = None


def clean_angles(angles, template=False):
    angles_ = []
    for a in angles:
        if template:
            _make_angle_template(a)
        else:
            a["code"] = deepcopy(a["code"])
            a["value"] = float(a["value"])
            
        angles_.append(a)
    return angles_

def _make_angle_template(a):
    a["code"] = []
    a['fc'] = None
    if "hfc" in a:
        del a["hfc"]
    a["value"] = None

def clean_dihedrals(dihedrals, template=False):
    dihedrals_ = []
    for d in dihedrals:
        if template:
            _make_dihedral_template(d)
        else:
            if "code" in d:
                d["code"] = deepcopy(d["code"])  
            else:
                print(d)  
            d["value"] = float(d["value"])
        dihedrals_.append(d)
    
    return dihedrals_

def _make_dihedral_template(d):
    d["code"] = []
    d['fc'] = None
    d["value"] = None
    d["essential"] = None
    
def clean_impropers(impropers, template=False):
    return []

def clean_rings(rings, template=False):
    rings_ = {}
    for r_id, r in list(rings.items()):
        rings_[r_id] = {"atoms": r["atoms"]}
    return rings_

def add_yml_comments(yml_str):
    atoms_comments = '''
#
# atom data structure
#    atom_id:
#        'index'  : index after rearrange
#        'uindex' : united atom index
#        'symbol' : atom symbol from pdb file
#        'type'   : atom type code
#        'group'  : residue name in pdb
#        'charge' : atom charge from QM
#        'coord'  : atom coordinates from pdb
#        'ocoord' : optimized coordinates from QM
#        'excl'   : exclusion list, use sorted atomed index rather than atomid, BE CAREFUL
#        'uexcl'  : united atom exclusion list
#        'equivalenceGroup'  : symmetry group that this atoms belongs to
'''.strip()

    bonds_comments = '''
#
# force_constant    Kb: force constant
# value_dist        b0: equilibrium bond length
#                    b: length of bond
#                    V: potential energy of bond
#
# functional form
#   fc: 
#       V = 1/2*Kb*(b - b0)^2       Harmonic
#   hfc:
#       V = 1/4*Kb*(b^2 - b0^2)^2   Quartic
#
# bonds data structure
#    'atoms' : atom IDs of this bond 
#    'value' : optimized length (nm)
#    'fc'    : force constant
#    'hfc'   : harmonic force constant
#    'code'  : [[mcb, fc, len] ....] possible parameters
#    'order' : bond order information from QM
'''.strip()
    angles_comments = '''
#
# force_constant      Ka: force constant
# value_angle       phi0: equilibrium angle
#                    phi: actual angle
#                      V: potential energy of angle
#
# functional form
#   fc: 
#       V = 1/2*Ka*( cos(phi) - cos(phi0) )^2
#   hfc:
#       V = 1/2*Ka*(phi - phi0)^2
#
# angles data structure
#    'atoms' : atom IDs of this angle
#    'value' : optimized bond angle (degree)
#    'code'  : [[mcb, fc, len] ....] possible parameters
#    'fc'    : force constant
#    'hfc'   : harmonic force constant
'''.strip()
    dihedrals_comments = '''
#
# force_constant      Kd: force constant for dihedral potential
# value_angle     theta0: equilibrium angle between planes 
# period               m: multiplicity 
# phase_shift      delta: phase shift
#                  theta: actual angle between planes
#                      V: potential energy of dihedral
#
# functional form
#       V = Kd*(1 + cos(delta)cos(m*theta))
#
# dihedrals data structure
#    'atoms'     : atom IDs of this dihedral
#    'value'     : optimized dihedral angle, in degree
#    'code'      : [[mcb, fc, phase shift, multiplicity, [possible atom type list]] ...
#    'fc'        : force constant
#    'phase'     : phase shift
#    'mul'       : multiplicity
#    'essential' : redudant or not
'''.strip()
    impropers_comments = '''
#
# force_constant      Kd: force constant for dihedral potential
# value_angle     theta0: equilibrium angle between planes 
# period               m: multiplicity 
# phase_shift      delta: phase shift
#                  theta: actual angle between planes
#                      V: potential energy of dihedral
#
# functional form
#       V = 1/2*Kd*(theta - theta0)^2
#
# improper dihedrals data structure
#    'atoms' : atom IDs of this improper dihedral
#    'ring'  : in which ring if in a ring
#    'code'  : improper type code (because there's no possible alternatives this is not a list)
'''.strip()
    rings_comments = '''
#
# rings data structure
#   ring_id:
#       'aromatic' : aromatic ring?
#       'atoms'    : list of atoms on the ring
'''.strip()
    yml_str = re.sub(re.compile(r"^atoms:", re.MULTILINE), "{0}\natoms:".format(atoms_comments), yml_str)
    yml_str = re.sub(re.compile(r"^angles:", re.MULTILINE), "{0}\nangles:".format(angles_comments), yml_str)
    yml_str = re.sub(re.compile(r"^bonds:", re.MULTILINE), "{0}\nbonds:".format(bonds_comments), yml_str)
    yml_str = re.sub(re.compile(r"^dihedrals:", re.MULTILINE), "{0}\ndihedrals:".format(dihedrals_comments), yml_str)
    yml_str = re.sub(re.compile(r"^impropers:", re.MULTILINE), "{0}\nimpropers:".format(impropers_comments), yml_str)
    yml_str = re.sub(re.compile(r"^rings:", re.MULTILINE), "{0}\nrings:".format(rings_comments), yml_str)
    return yml_str
