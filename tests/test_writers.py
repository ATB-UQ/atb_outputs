"""Characterisation tests of the LIVE writers (review 2026-10 wave 2).

Goldens in tests/golden/ were produced by the writers as they stood at the
time (tests/make_golden.py) -- they pin current behaviour, they do not claim it
is right (see the review findings C-xx).
"""
import datetime
import io

import pytest
import yaml

import atb_outputs.pdb as PDB
from atb_outputs.formats import pdb, g96, yml, lgf, ccd_cif, mol_data_dict
from atb_outputs.mol_data import MolData, MolDataFailure, mol_data_from_mol_data_dict

from conftest import golden, mask_date


def test_pdb_golden(mol_data):
    assert mask_date(pdb(mol_data)) == golden("toluene.pdb")


def test_g96_golden(mol_data):
    assert g96(mol_data) == golden("toluene.g96")


def test_lgf_golden(mol_data):
    out = lgf(mol_data)
    assert [name for name, _ in out] == ["lgf"]
    assert out[0][1] == golden("toluene.lgf")


def test_ccd_cif_golden(mol_data):
    assert ccd_cif(mol_data, "TOL", "TOL") == golden("toluene.cif")


def test_yml_round_trips_through_safe_load(mol_data):
    text = yml(mol_data)
    loaded = yaml.safe_load(text)
    assert set(loaded) == set(mol_data_dict(mol_data))
    assert len(loaded["atoms"]) == len(mol_data.atoms) == 15
    assert loaded["_dihedrals"] is None
    assert loaded["var"]["molid"] == "21"
    # the explanatory comment block is spliced in after the first line
    assert "# force_constant      Ka: force constant" in text


def test_yml_prints_nothing_to_stdout(mol_data, capsys):
    yml(mol_data)
    assert capsys.readouterr().out == ""


def test_mol_data_dict_exact_key_set(mol_data):
    """core hashes this dict for the topology cache key (core outputs.py
    ``_hashed_mol_data_dict``), so ANY change to the key set -- or to the value
    of ``_dihedrals`` -- invalidates every cached topology.

    ``'_dihedrals'`` is written ``mol_data._dihedrals if hasattr(mol_data,
    '_dihedral') else None`` (typo: ``_dihedral``), so it is ALWAYS None.
    That typo is load-bearing: do NOT "fix" it without rehashing the whole
    cache in the same change (review 2026-10 atb_outputs C-04).
    """
    d = mol_data_dict(mol_data)
    assert set(d) == {"atoms", "bonds", "angles", "dihedrals", "_dihedrals",
                      "impropers", "rings", "var"}
    assert d["_dihedrals"] is None
    # even when the attribute the typo should have looked at exists
    mol_data._dihedrals = ["something"]
    assert mol_data_dict(mol_data)["_dihedrals"] is None
    assert all("flavour" not in a for a in d["atoms"].values())


# ---- _readPDB ---------------------------------------------------------------

def _atom(serial, name, el, x, y, z):
    return ("HETATM%5d %-4s %3s %5d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
            % (serial, name, "MOL", 1, x, y, z, el))


WATER = "\n".join([
    _atom(1, "O", "O", 0.0, 0.0, 0.0),
    _atom(2, "H1", "H", 0.96, 0.0, 0.0),
    _atom(3, "H2", "H", -0.24, 0.93, 0.0),
    "CONECT    1    2    3",   # one-sided: 2 and 3 have no CONECT of their own
    "END",
])


def test_readpdb_conect_is_made_symmetric():
    m = MolData(WATER, build_ring=False)
    assert m.atoms[1]["conn"] == [2, 3]
    assert m.atoms[2]["conn"] == [1]
    assert m.atoms[3]["conn"] == [1]
    assert m.atoms[1]["coord"] == pytest.approx([0.0, 0.0, 0.0])
    assert m.atoms[2]["coord"] == pytest.approx([0.096, 0.0, 0.0])  # Angstrom -> nm
    assert m.atoms[2]["type"] == "H"
    assert [a["id"] for a in m.atoms.values()] == [1, 2, 3]
    assert len(m.bonds) == 2


def test_readpdb_single_atom_allowed_without_conect():
    m = MolData(_atom(1, "NA", "NA", 0, 0, 0), build_ring=False)
    assert m.atoms[1]["conn"] == []
    assert m.bonds == []


@pytest.mark.parametrize("text", ["", "END\n", "garbage", None])
def test_readpdb_no_atoms_raises(text):
    with pytest.raises(MolDataFailure, match="contains no atoms"):
        MolData(text, build_ring=False)


def test_readpdb_unbonded_atom_in_multi_atom_raises():
    pdb_str = "\n".join([_atom(1, "O", "O", 0, 0, 0), _atom(2, "H1", "H", 1, 0, 0),
                         _atom(3, "H2", "H", 0, 1, 0), "CONECT    1    2"])
    with pytest.raises(MolDataFailure, match="Missing connectivities"):
        MolData(pdb_str, build_ring=False)


def test_readpdb_orphan_conect_raises():
    pdb_str = "\n".join([_atom(1, "O", "O", 0, 0, 0), _atom(2, "H1", "H", 1, 0, 0),
                         "CONECT    1    2   99"])
    with pytest.raises(MolDataFailure, match="non-existent atom"):
        MolData(pdb_str, build_ring=False)


def test_readpdb_unbonded_atom_allowed_without_single_molecule_enforcement():
    # C-09: enforce_single_molecule=False used to raise KeyError 'conn'.
    pdb_str = "\n".join([_atom(1, "O", "O", 0, 0, 0), _atom(2, "H1", "H", 1, 0, 0),
                         _atom(3, "H2", "H", 5, 5, 5), "CONECT    1    2"])
    md = MolData(pdb_str, build_ring=False, enforce_single_molecule=False)
    assert md.atoms[3]["conn"] == []
    assert md.atoms[1]["conn"] == [2] and md.atoms[2]["conn"] == [1]
    assert len(md.bonds) == 1


# ---- formerly strict xfails (review 2026-10 C-02, C-06), fixed in wave 2 -----

def test_mol_data_from_mol_data_dict_loads(mol_data):
    # C-02: the dict loader built MolData(None), which the empty-PDB guard rejects.
    md = mol_data_from_mol_data_dict(mol_data_dict(mol_data))
    assert len(md.atoms) == 15
    md.var = dict(md.var)
    md.completed = mol_data.completed
    assert g96(md) == golden("toluene.g96")


def test_pdb_header_year_is_calendar_year(mol_data, monkeypatch):
    # C-06: the header used %g (ISO-week year): 2027-01-01 came out as '01-Jan-26'.
    class FakeDate(datetime.date):
        @classmethod
        def today(cls):
            return cls(2027, 1, 1)

    monkeypatch.setattr(PDB, "date", FakeDate)
    first = pdb(mol_data).splitlines()[0]
    assert first.endswith("01-Jan-27")
