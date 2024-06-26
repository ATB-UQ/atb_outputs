from logging import Logger
from sys import stderr
from copy import deepcopy
from math import sqrt
from typing import Optional, Dict, Any, List

from atb_outputs.helpers.types_helpers import Atom, Tuple, Ring, Coordinate
from atb_outputs.helpers.dijkstra import shortestPath


class MolDataFailure(Exception):
    pass


class MolData(object):
    def __init__(self,
                 initialiser_object: Optional[str],
                 log: Optional[Logger] = None,
                 build_ring: bool = True,
                 enforce_single_molecule: bool = True) -> None:
        """
        TODO: MolData docs
        :param initialiser_object: One of a Molecule3D, FDB_Molecule or a pdb_string, each has a different initialiser
        :param log: ?
        :param build_ring: ?
        :param enforce_single_molecule: ?
        :param atom_index_name: If using a Molecule3D to initialise, which atom ID to use for the mapping. Needs to be
            and ID of type int
        """

        #### Added method for instanciating from FDBMolecule objects

        self.atoms = {}
        self.bonds = []
        self.equivalenceGroups = {}

        self._atom_index_map = {}

        if type(initialiser_object).__name__ == 'FDBMolecule':
            self._readFDBMolecule(initialiser_object)
        elif type(initialiser_object).__name__ == 'Molecule3D':
            self._readMolecule3D(initialiser_object)
        else:
            self._readPDB(initialiser_object, enforce_single_molecule=enforce_single_molecule)

        if build_ring:
            self.rings = build_rings(self, log)

        self.united_hydrogens = []

    def unite_atoms(self) -> None:
        self.united_hydrogens = []
        for atom in self.atoms.values():
            connected_hydrogens = [a_id for a_id in atom["conn"] if self.atoms[a_id]["type"] == "H"]
            if atom["type"] == "C" and len(connected_hydrogens) > 1:
                self.united_hydrogens.extend([a_id for a_id in connected_hydrogens])

        united_atoms = []
        for atom_id, atom in sorted(self.atoms.items()):
            if atom_id not in self.united_hydrogens:
                united_atoms.append(atom_id)
                self.atoms[atom_id]["uindex"] = len(united_atoms)
        self._unite_bonds()

    def _unite_bonds(self) -> None:
        for bond in self.bonds:
            if any([a in self.united_hydrogens for a in bond["atoms"]]):
                continue
            else:
                bond["united"] = True

    def get_id(self, index: int) -> int:
        '''return id of the atom with specified index number'''
        return [k for (k, v) in self.atoms.items() if v['index'] == index][0]

    def __getitem__(self, atom_id: int) -> Atom:
        '''return an atom with atom_id. '''
        assert type(atom_id) == int, 'Atom identifiers are integers'
        try:
            return self.atoms[atom_id]
        except:
            raise Exception('atom with id number %d not found.' % atom_id)

    def _addBondData(self, atm1: int, atm2: int) -> None:
        if atm1 == atm2:
            return
        if set([atm1, atm2]) in [set(b["atoms"]) for b in self.bonds]:
            return
        self.bonds.append({"atoms": [int(atm1), int(atm2)]})

    def _readFDBMolecule(self, Molecule: 'FDBMolecule') -> None:
        fdb_atom_info = Molecule.atom_info
        for a in Molecule.atoms:
            index_dict = Molecule.atoms[a]._index

            bond_info = Molecule.bonds(a)
            connectivity = [Molecule.atoms[a_bonded]._index['id'] for _, a_bonded in bond_info]
            i = index_dict['id']
            name = index_dict['name']
            fdb_info = [a_fdb for a_fdb in fdb_atom_info if a_fdb['elementID'] == name]
            assert len(fdb_info) == 1, fdb_info
            fdb_info = fdb_info[0]

            self.atoms[i] = {'id': i,
                             'index': i,
                             'symbol': name,
                             'type': Molecule.atoms[a].element,
                             'conn': connectivity,
                             'coord': [fdb_info['x3d'], fdb_info['y3d'], fdb_info['z3d']]
                             }

            self.bonds = [
                {'atoms': [Molecule.atoms[a1]._index['id'], Molecule.atoms[a2]._index['id']]}
                for a1, a2 in Molecule.bonds]

    def _readMolecule3D(self, Molecule: 'Molecule3D') -> None:
        if tuple(sorted(Molecule.atoms)) != sorted(range(1,len(Molecule.atoms)+1)):
            map_id_flag = True
            # I don't fully understand the indexing approach, but the naughty interface expects objects
            # Initialised from the ATB to have indexes starting from 1 and offsets them to start at 0, so
            # I'm mapping the atom names to a 1 index in this instance
            self._atom_index_map = {
                a: i+1 for i, a in enumerate(Molecule.atoms)
            }
        else:
            map_id_flag = False

        for a in Molecule.atoms:
            atom_obj = Molecule.atoms[a]
            # index_dict = Molecule.atoms[a]._index
            name = atom_obj.name
            bond_info = Molecule.bonds(a)

            if map_id_flag:
                # if the ids are getting mapped by iter, then add the original atom id to the dict
                extra_id = {'Molecule3D_id': atom_obj.get_index()}

                mold_data_index = self._atom_index_map[name]

            else:
                extra_id = {}
                mold_data_index = atom_obj.get_index()

            connectivity = [self._atom_index_map[Molecule.atoms[a_bonded].get_index()] if map_id_flag else Molecule.atoms[a_bonded].get_index()
                            for _, a_bonded in bond_info]

            self.atoms[mold_data_index] = {**{'id': mold_data_index,
                             'index': mold_data_index,
                             'symbol': name,
                             'type': atom_obj.element,
                             'conn': connectivity,
                             'coord': atom_obj.coordinates
                             }, **extra_id}

        self.bonds = [
            {'atoms': [
                self._atom_index_map[Molecule.atoms[a1].get_index()] if map_id_flag else Molecule.atoms[a1].get_index(),
                self._atom_index_map[Molecule.atoms[a2].get_index()] if map_id_flag else Molecule.atoms[a2].get_index()]}
            for a1, a2 in Molecule.bonds]
        # if map_id_flag:
        #     for a1, a2 in Molecule.bonds
        #     assert self._atom_index_map[Molecule.atoms[a1].get_index()] if map_id_flag

    def _readPDB(self, pdb_str: Optional[str], enforce_single_molecule: bool = True) -> None:
        '''Read lines of PDB files'''
        assert pdb_str is None or isinstance(pdb_str, str), type(pdb_str)

        pdbDict = {}
        for line in ([] if pdb_str is None else pdb_str.splitlines()):
            # lines with atom coordinates
            if line.startswith("ATOM") or line.startswith('HETATM'):
                # split line for different fields
                # it = line.split()
                # in case of some fields are missing, read according to pdb standard
                # if len(it) < 11:
                it = [line[0:6], line[6:11], line[12:16], line[17:20],
                      line[22:27], line[30:38], line[38:46], line[46:54], line[54:60],
                      line[60:66], line[76:78]]
                it = [i.strip() for i in it]
                # store information that we are interested in. Refer to MoleculeData
                # class for more details, also here we have a angstrom to nm conversion
                pdbDict[int(it[1])] = {
                    'index': int(it[1]),
                    'symbol': it[2],
                    'group': it[3],
                    'coord': [float(it[5]) / 10., float(it[6]) / 10., float(it[7]) / 10.],
                    'pdb': line.strip(),
                    'type': it[-1].upper(),
                }
            # connectivity records
            if line.startswith('CONECT'):
                it = line.split()
                for key, item in pdbDict.items():
                    if key == int(it[1]):
                        conn_list = []
                        for i in it[2:6]:
                            try:
                                num = int(i)
                                if num == 0:
                                    break
                                conn_list.append(num)
                            except Exception:
                                break
                        if 'conn' in item:
                            item['conn'].extend(conn_list)
                        else:
                            item['conn'] = conn_list
                        break

        # make sure connectivity information is symmetric by simply mirroring connections
        # flag any orphan connectivities for removal
        orphanAtomReference = {}
        for key in pdbDict:  # for each atom in pdb
            if 'conn' in pdbDict[key]:  # if this atom lists connections to others
                for conn in pdbDict[key]['conn']:  # for each connection to others
                    if conn not in pdbDict:  # if it connects to a non-existent atom
                        stderr.write("connectivity made from atom %s to non-existent atom %s!" % (key, conn))
                        if key not in orphanAtomReference:
                            orphanAtomReference[key] = [conn]
                        else:
                            orphanAtomReference[key].append(conn)
                        continue
                    if 'conn' in pdbDict[conn]:  # if the other atom has a list of connections
                        pdbDict[conn]['conn'].append(key)  # append this atom to the end of the other atom's list
                    else:
                        pdbDict[conn]['conn'] = [key]  # create new list containing this atom

        # remove any orphan connection records
        for (k, connList) in orphanAtomReference.items():
            for c in connList:
                pdbDict[k]["conn"].remove(c)
                error_msg = "connectivity made from %s to non-existent atom %s removed" % (k, c)

                stderr.write(error_msg)
                raise MolDataFailure(error_msg)

        has_connects = lambda atom: 'conn' in atom and atom['conn']
        if not all(has_connects(atom) for atom in pdbDict.values()):
            if len(pdbDict) == 1:
                # Only single atom molecules are allowed to have no bonds
                list(pdbDict.values())[0]['conn'] = []
            else:
                if enforce_single_molecule:
                    raise MolDataFailure(
                        'Mol_Data Error: Missing connectivities for atoms {0}'.format(
                            [atom['index'] for atom in pdbDict.values() if not has_connects(atom)],
                        ),
                    )
                else:
                    pass

        # sort and unique connectivities
        for (ID, atom) in pdbDict.items():
            atom['conn'] = sorted(set(atom['conn']))
            for neighbour in atom['conn']:
                self._addBondData(ID, neighbour)

        self.atoms = pdbDict
        for (atom_id, atom) in self.atoms.items():
            atom['id'] = atom_id

    def __str__(self) -> str:
        return 'MolData(atoms={0}, bonds={1})'.format(
            self.atoms,
            self.bonds,
        )


def mol_data_from_mol_data_dict(mol_data_dict: Dict[str, Any]) -> MolData:
    mol_data = MolData(None)

    for key in mol_data_dict.keys():
        setattr(mol_data, key, mol_data_dict[key])

    return mol_data


def build_rings(data: MolData, log: Optional[Logger] = None) -> Dict[int, Ring]:
    def _is_ring_in_all_rings(ring: Any, all_rings: List[Any]) -> bool:
        for existing_ring in all_rings.values():
            if frozenset(existing_ring["atoms"]) == frozenset(ring):
                return True
        return False

    all_rings = {}
    ring_count = 1
    mol_graph = _get_graph_dict(data.atoms)

    # serialized_graph = []
    # _serialize_weighted_graph(mol_graph,output=serialized_graph)
    # log.debug("connectivity_graph:\n" + "\n".join(serialized_graph))

    for bond in data.bonds:
        rings = _get_all_rings_for_bond(deepcopy(mol_graph), bond["atoms"])
        for ring in rings:
            if len(ring) == 2:
                continue
            ring = list(map(int, ring))
            if not _is_ring_in_all_rings(ring, all_rings):
                ring_dict = {"atoms": ring, "aromatic": False}
                ring_dict["aromatic"] = is_ring_aromatic(data, ring_dict, log)
                all_rings[ring_count] = ring_dict
                ring_count += 1
    return all_rings


ACCEPTED_PLANAR_VALENCE_PER_ATOM_TYPE = {
    'C': [3],
    'N': [2, 3],  # For pyridine
    'O': [2],
    'S': [2],
}
PLANAR_DISTANCE_TOL = 0.025


def is_ring_aromatic(data: MolData, ring: Ring, log: Logger) -> bool:
    return has_ring_planar_geometry(data, ring, log) and has_ring_planar_valences(data, ring, log)


def has_ring_planar_geometry(data: MolData, ring: Ring, log: Logger) -> bool:
    if len(ring["atoms"]) < 4:
        return False
    else:
        # Get dihedral atoms
        return is_ring_planar(data, ring, log)


def has_ring_planar_valences(data: MolData, ring: Ring, log: Logger) -> bool:
    ''' Prevents rings with unual valences to be assigned as planar'''

    ring_atoms = [data.atoms[atom_id] for atom_id in ring['atoms']]
    has_aromatic_valences = True
    for atom in ring_atoms:
        for atom_type, accepted_valences in ACCEPTED_PLANAR_VALENCE_PER_ATOM_TYPE.items():
            if atom['type'].upper() == atom_type:
                if not len(atom['conn']) in accepted_valences: has_aromatic_valences = False
                break
        if not has_aromatic_valences: break
    if has_aromatic_valences:
        if log: log.debug("{0} has the valences expected for a planar ring and will be treated as such".format(
            [data[x]["symbol"] for x in ring['atoms']]))
    else:
        if log: log.debug(
            "{0} DOES NOT have the valences expected for a planar ring and will not be treated as such".format(
                [data[x]["symbol"] for x in ring['atoms']]))
    return has_aromatic_valences


def _serialize_weighted_graph(G: Any, indent: str = "", output: List[Any] = []) -> None:
    for node, branches in G.items():
        line = indent + str(node)
        if type(branches) is dict:
            output.append(line + "--")
            _serialize_weighted_graph(branches, indent=" " * 4, output=output)
        else:
            line += ": {0}".format(branches)
            output.append(line)


def is_ring_planar(data: MolData, ring: Ring, log: Logger) -> bool:
    atoms = data.atoms
    coord_type = "ocoord" if "ocoord" in list(atoms.values())[0] else "coord"

    ring_atoms = [atoms[atom_id] for atom_id in ring["atoms"]]
    A, B, C, D = equation_of_plane(
        ring_atoms[0][coord_type],
        ring_atoms[1][coord_type],
        ring_atoms[2][coord_type],
    )
    max_distance = 0
    for atom in ring_atoms[3:]:
        distance = _distance_from_plane(A, B, C, D, atom[coord_type])
        max_distance = max(abs(distance), max_distance)
        if abs(distance) > PLANAR_DISTANCE_TOL:
            return False
    if log:
        log.debug(
            "Maximum distance to plane is {0:.3f}nm ({1})".format(
                max_distance,
                [data[x]["symbol"] for x in ring["atoms"]],
            ),
        )
    return True


def equation_of_plane(a0: Coordinate, a1: Coordinate, a2: Coordinate) -> Tuple[float, float, float, float]:
    det1 = ((a1[1] - a0[1]) * (a2[2] - a0[2])) - ((a2[1] - a0[1]) * (a1[2] - a0[2]))
    det2 = ((a1[2] - a0[2]) * (a2[0] - a0[0])) - ((a2[2] - a0[2]) * (a1[0] - a0[0]))
    det3 = ((a1[0] - a0[0]) * (a2[1] - a0[1])) - ((a2[0] - a0[0]) * (a1[1] - a0[1]))
    D = det1 * a0[0] + det2 * a0[1] + det3 * a0[2]
    D = -D
    return (det1, det2, det3, D)


def _distance_from_plane(A: float, B: float, C: float, D: float, pt: Coordinate) -> float:
    x, y, z = pt
    denom = sqrt(A ** 2 + B ** 2 + C ** 2)
    if not denom: return 0
    numer = A * x + B * y + C * z + D
    distance = numer / denom
    return distance


def _get_all_rings_for_bond(mol_graph: Any, bond_atom_ids: Any) -> List[Any]:
    i0 = str(bond_atom_ids[0])
    i1 = str(bond_atom_ids[1])
    all_rings = []
    mol_graph[i0][i1] = 999
    found = True
    while found:
        found = False
        ring = shortestPath(
            mol_graph,
            i0,
            i1,
        )
        if not ring in all_rings:
            all_rings.append(ring)
            found = True
            for i in range(len(ring) - 1):
                mol_graph[ring[i]][ring[i + 1]] = 2
    mol_graph[i0][i1] = 1
    return all_rings


def _get_graph_dict(atoms: Dict[int, Atom]) -> Dict[str, Atom]:
    G = {}
    for (atom_id, atom) in atoms.items():
        tmp = {}
        for i in atom["conn"]:
            tmp[str(i)] = 1
        G[str(atom_id)] = tmp
    return G
