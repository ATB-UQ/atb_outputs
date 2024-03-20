from functools import reduce
from io import StringIO

GROMOS_IMPROPER_DIHEDRALS = {
1: {'fc': 0.0510, 'value': 0.0       },
2: {'fc': 0.102,  'value': 35.26439  },
3: {'fc': 0.204,  'value': 0.0       },
4: {'fc': 0.0510, 'value': 180.0     },
5: {'fc': 0.102,  'value': -35.26439 },
}


def has_all_charges(data):
    return all( ["charge" in atom for atom in list(data.atoms.values())] )


def first_code(obj):
    return obj['code'][0]


def to_aaindex(atom_id, data):
    return data.atoms[atom_id]['index']


def to_uindex(atom_id, data):
    return data.atoms[atom_id]['uindex']


def to_atom_id(atom_index, data, united):

    def uindex_to_id_filter(item):
        ref_index = item[1]["uindex"] if "uindex" in item[1] else item[1]["index"]
        return ref_index == atom_index

    def index_to_id_filter(item):
        return item[1]["index"] == atom_index

    index_filter = uindex_to_id_filter if united else index_to_id_filter

    return list(filter(index_filter, list(data.atoms.items())))[0][0]


def a_ljsym(atom, united_atom_prefix):
    return atom[united_atom_prefix + 'ljsym'] if (united_atom_prefix + 'ljsym') in atom else atom['ljsym']


def a_charge(atom, united_atom_prefix):
    return atom[united_atom_prefix + 'charge'] if (united_atom_prefix + 'charge') in atom else (atom['charge'] if 'charge' in atom else None )


def header(data, io):
    # Write header
    print(";", file=io)


def title(data, io):
    # Write 'title'
    #print >> io, '\n'.join([';\t'+l for l in data['title']])
    print(';', file=io)


def moleculetype(data, io):
    # Write 'moleculetype' block
    print('[ moleculetype ]', file=io)
    print('; Name   nrexcl', file=io)
    print(data.var['rnme'] + (9-len(data.var['rnme']))*' ' + '3', file=io)


def atoms(data, io, united):
    united_atom_prefix = 'u' if united else ''
    # Write 'atoms' block
    print('[ atoms ]', file=io)
    print(';  nr  type  resnr  resid  atom  cgnr  charge    mass', file=io)
    a_rnme = data.var['rnme']
    totalcharge = 0.0

    has_charges = has_all_charges(data)

    atom_format = '%8.3f' if has_charges else '%8s'
    if united:
        sorted_atoms = sorted(
            [a for a in data.atoms.values() if 'uindex' in a],
            key=lambda atom: atom['uindex'],
        )
    else:
        sorted_atoms = sorted(
            data.atoms.values(),
            key=lambda atom: atom['index'],
         )

    atom_lines = []

    for atom in sorted_atoms:
        a_symbol = atom['symbol']
        a_index = atom[united_atom_prefix + 'index']
        a_cgroup = atom['cgroup']
        a_mass = atom[united_atom_prefix + 'mass'] if (united_atom_prefix + 'mass') in atom else atom['mass']

        if has_charges:
            atom_charge = a_charge(atom, united_atom_prefix)
        else:
            atom_charge = '%%'

        atom_lines.append(
            '%5d %5s %4s %7s %6s %4d {0} %8.4f'.format(atom_format) % (
                a_index,
                a_ljsym(atom, united_atom_prefix),
                '1',
                a_rnme,
                a_symbol,
                a_index,
                atom_charge,
                a_mass,
            ),
        )
        if has_charges:
            totalcharge += atom_charge
    atom_lines.append('; total charge of the molecule: %7.3f' % totalcharge)

    print("\n".join(atom_lines), file=io)


def bonds(data, io, united):
    # Write 'bonds' block
    print('[ bonds ]', file=io)
    print(';  ai   aj  funct   c0         c1', file=io)
    to_index = to_uindex if united else to_aaindex
    for b in data.bonds:
        if united and 'united' in b: continue
        if 'value' in b and len(b['code']) > 0 :
            print('%5d %4d %4s %8.4f %12.4e' \
                    %( to_index(b['atoms'][0], data),
                       to_index(b['atoms'][1], data),
                       '2',
                       first_code(b)['value'],
                       first_code(b)['fc']), file=io)
        else:
            print('%5d %4d %4s %8s %12s' \
                    %( to_index(b['atoms'][0], data),
                       to_index(b['atoms'][1], data),
                       '2',
                       '%%',
                       '%%'), file=io)


def dummy_blocks(data, io):
    print('[ pairs ]', file=io)
    print(';  ai   aj  funct  ;  all 1-4 pairs but the ones excluded in GROMOS itp', file=io)
    print('[ angles ]', file=io)
    print(';  ai   aj   ak  funct   angle     fc', file=io)
    print('[ dihedrals ]', file=io)
    print('; GROMOS improper dihedrals', file=io)
    print(';  ai   aj   ak   al  funct   angle     fc', file=io)
    print('[ dihedrals ]', file=io)
    print(';  ai   aj   ak   al  funct    ph0      cp     mult', file=io)


def pairs_1_4(data, io, nbr3, united):

    # Print 1-4 pairs
    print('[ pairs ]', file=io)
    print(';  ai   aj  funct  ;  all 1-4 pairs but the ones excluded in GROMOS itp', file=io)
    united_atom_prefix = 'u' if united else ''
    for k, v in sorted(list(nbr3.items())):
        atom = data.atoms[to_atom_id(k, data, united)]
        for n in v:
            should_be_displayed = n not in atom[united_atom_prefix + 'excl']
            if should_be_displayed:
                print('%5d %4d %4s' %(
                          k,
                          n,
                          '1'), file=io)


def angles(data, io, united):
    # Print angle block
    print('[ angles ]', file=io)
    print(';  ai   aj   ak  funct   angle     fc', file=io)
    to_index = to_uindex if united else to_aaindex
    for angle in data.angles:
        if united and 'united' in angle: continue
        if 'value' in angle and 'code' in angle and len(angle['code']) > 0:
            print('%5d %4d %4d %4s %9.2f %8.2f' \
                    %( to_index(angle['atoms'][0], data),
                       to_index(angle['atoms'][1], data),
                       to_index(angle['atoms'][2], data),
                       '2',
                       first_code(angle)['value'],
                       first_code(angle)['fc']), file=io)
        else:
            print('%5d %4d %4d %4s %9s %8s' \
                    %( to_index(angle['atoms'][0], data),
                       to_index(angle['atoms'][1], data),
                       to_index(angle['atoms'][2], data),
                       '2',
                       '%%',
                       '%%'), file=io)


def impropers(data, io, united):
    # Print improper dihedral block
    print('[ dihedrals ]', file=io)
    print('; GROMOS improper dihedrals', file=io)
    print(';  ai   aj   ak   al  funct   angle     fc', file=io)
    to_index = to_uindex if united else to_aaindex
    for i in data.impropers:
        if united and 'united' in i: continue
        # if this is an all atom output, skip all type 2 impropers
        if not united and i['code'] == 2:
            continue
        #force constant has to be converted from kJ/mol/deg^2 to kJ/mol/rad^2
        print('%5d %4d %4d %4d %4s %9.2f %8.2f' \
                %( to_index(i['atoms'][0], data),
                   to_index(i['atoms'][1], data),
                   to_index(i['atoms'][2], data),
                   to_index(i['atoms'][3], data),
                   '2',
                   GROMOS_IMPROPER_DIHEDRALS[ i['code'] ]['value'],
                   GROMOS_IMPROPER_DIHEDRALS[ i['code'] ]['fc']*3281.5686
                  ), file=io)


def dihedrals(data, io, united):
    # Print dihedral block
    print('[ dihedrals ]', file=io)
    print(';  ai   aj   ak   al  funct    ph0      cp     mult', file=io)
    to_index = to_uindex if united else to_aaindex
    for i in data.dihedrals:
        if united and 'united' in i:
            # print("skipping united")
            continue
        # Only print essential dihedrals
        if 'essential' in i and not i['essential']:
            # print("skipping essential")
            continue
        if 'code' in i and len(i['code']) > 0:
            print('%5d %4d %4d %4d %4s %9.2f %8.2f %4d' \
                    %( to_index(i['atoms'][0], data),
                       to_index(i['atoms'][1], data),
                       to_index(i['atoms'][2], data),
                       to_index(i['atoms'][3], data),
                       '1',
                       first_code(i)['value'],
                       first_code(i)['fc'],
                       first_code(i)['mul']
                     ), file=io)
        else:
            print('%5d %4d %4d %4d %4s %9s %8s %4s' \
                    %( to_index(i['atoms'][0], data),
                       to_index(i['atoms'][1], data),
                       to_index(i['atoms'][2], data),
                       to_index(i['atoms'][3], data),
                       '1',
                       '%%',
                       '%%',
                       '%%'
                     ), file=io)


def graph_dihedrals(data, io, united):
    # Print dihedral block
    print('[ dihedrals ]', file=io)
    print(';  ai   aj   ak   al  funct    ph0      cp     mult', file=io)
    to_index = to_uindex if united else to_aaindex
    if united:
        dih_key = "united_atom"
    else:
        dih_key = "all_atom"
    for i in data.graph_dihedrals[dih_key]:
        # print(i)
        if united and 'united' in i:
            # print("skipping united")
            continue
        # Only print essential dihedrals
        if 'essential' in i and not i['essential']:
            # print("skipping essential")
            continue
        if 'code' in i and len(i['code']) > 0:
            print('%5d %4d %4d %4d %4s %9.2f %8.2f %4d' \
                    %( to_index(i['atoms'][0], data),
                       to_index(i['atoms'][1], data),
                       to_index(i['atoms'][2], data),
                       to_index(i['atoms'][3], data),
                       '1',
                       first_code(i)['value'],  # phase
                       first_code(i)['fc'],     # force constant
                       first_code(i)['mul']     # multiplicity
                     ), file=io)
        else:
            print('%5d %4d %4d %4d %4s %9s %8s %4s' \
                    %( to_index(i['atoms'][0], data),
                       to_index(i['atoms'][1], data),
                       to_index(i['atoms'][2], data),
                       to_index(i['atoms'][3], data),
                       '1',
                       '%%',
                       '%%',
                       '%%'
                     ), file=io)


def exclusions_1_4(data, io, nbr3, united):
    # Print 1-4 exclusions
    print('[ exclusions ]', file=io)
    print(';  ai   aj  funct  ;  GROMOS 1-4 exclusions', file=io)
    united_atom_prefix = 'u' if united else ''
    for k, v in sorted(list(nbr3.items())):
        atom = data.atoms[to_atom_id(k, data, united)]
        for n in v:
            should_be_displayed = n in atom[united_atom_prefix + 'excl']
            if should_be_displayed:
                print('%5d %4d' %(k, n), file=io)


def calculate_1_4_neighbours(data, united):
    # Calculate all 1-4 neighbourship
    if not united:
        conn_index = dict([(a["index"], _conv_to_index(a["conn"], data, united)) for a in list(data.atoms.values())])
    else:
        conn_index = dict([(a["uindex"], _conv_to_index(a["uconn"], data, united)) for a in list(data.atoms.values()) if
                           'uindex' in a])
    nbr3 = {}
    for index, conn in list(conn_index.items()):
        # 1st neighbor
        first_neighbours = sorted(conn)
        # 2nd neighbor
        second_neighbours = [x for x in _get_neighbour_ids(first_neighbours, conn_index) if x not in first_neighbours]
        # 3rd neighbor
        third_neighbours = [x for x in _get_neighbour_ids(second_neighbours, conn_index) if
                            x not in first_neighbours + second_neighbours]

        nbr3[index] = [x for x in third_neighbours if x > index]
    return nbr3


def _get_neighbour_ids(atom_ids, conn_index):
    return sorted( list( set( reduce(lambda x,y:x+y, [conn_index[atom_id] for atom_id in atom_ids], []) ) ) )


def _conv_to_index(atom_ids, data, united):
    to_index = to_uindex if united else to_aaindex
    return [to_index(x, data) for x in atom_ids]


def itp(data, united=False):
    nbr3 = calculate_1_4_neighbours(data, united)

    io = StringIO()
    header(data, io)
    title(data, io)
    moleculetype(data, io)
    atoms(data, io, united)
    bonds(data, io, united)
    pairs_1_4(data, io, nbr3, united)
    angles(data, io, united)
    impropers(data, io, united)
    # graph_dihedrals(data, io, united)
    dihedrals(data, io, united)
    exclusions_1_4(data, io, nbr3, united)

    return io.getvalue()
