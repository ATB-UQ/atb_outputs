from io import StringIO

def graph(molecule_data, enforce_has_charges: bool = True, enforce_has_ocoords: bool = True):
    io = StringIO()

    print('''@nodes
partial_charge  label   label2  atomType    coordX  coordY  coordZ  initColor''', file=io)

    atoms = list(molecule_data.atoms.values())

    if enforce_has_charges:
        assert all([('charge' in atom) for atom in atoms]), '''ERROR: Can't output .lgf because of missing charges.'''

    if enforce_has_ocoords:
        assert all([('ocoord' in atom) for atom in atoms]), '''ERROR: Can't output .lgf because of missing optimised coordinates.'''
        coordinate_key = 'ocoord'
    else:
        coordinate_key = 'coord'

    for atom in atoms:
        print('''{0:4.3f} {1} {2} {3} {4:4.3f} {5:4.3f} {6:4.3f} {7}'''.format(
            atom['charge'] if 'charge' in atom else 0.0,
            atom['id'],
            atom['symbol'],
            atom['iacm'],
            atom[coordinate_key][0],
            atom[coordinate_key][1],
            atom[coordinate_key][2],
            atom['cgroup'],
        ), file=io)

    print('''@edges
        label''', file=io)

    for (i, bond_atoms) in enumerate([bond['atoms'] for bond in molecule_data.bonds]):
        print('''{0} {1} {2}'''.format(
            bond_atoms[0],
            bond_atoms[1],
            i,
        ), file=io)

    return io.getvalue()

if __name__ == '__main__':
    example_file = '''
@nodes
partial_charge  label   label2  atomType    coordX  coordY  coordZ  initColor   
0.129   1   H1  20  0.73    2.376   -0.001  11  
0.129   10  H2  20  -1.692  1.82    -0.001  2   
-0.129  11  C3  12  -1.363  -0.312  -0.001  1   
-0.129  12  C2  12  -0.952  1.024   0.001   0   
-0.129  2   C6  12  1.363   0.312   0.001   10  
0.129   3   H6  20  2.423   0.554   -0.001  9   
-0.129  4   C1  12  0.411   1.337   0.003   8   
-0.129  5   C5  12  0.952   -1.024  -0.001  7   
-0.129  6   C4  12  -0.411  -1.337  -0.003  6   
0.129   7   H5  20  1.692   -1.82   0.001   5   
0.129   8   H4  20  -0.73   -2.376  0.001   4   
0.129   9   H3  20  -2.423  -0.554  0.001   3   
@edges
        label   
1   4   0   
2   3   1   
2   4   2   
2   5   3   
4   12  4   
5   6   5   
5   7   6   
6   8   7   
6   11  8   
9   11  9   
10  12  10  
11  12  11  '''
