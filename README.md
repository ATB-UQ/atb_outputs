# atb_outputs

The ATB molecule data object (`MolData`) and the writers that turn it into output
files (PDB, G96, YML, pickle, LGF, graph image, CCD CIF, MOL2, GROMACS ITP).
Pure Python, installable on its own (`pip install -e .`), but note that
`setup.cfg` declares **no dependencies**: it needs `numpy`, `pyyaml` and, for
`mol2()`, the `chemistry_helpers` sibling (see Configuration).

**Status:** live-support. Imported by `core` (the topology generator), `website`,
`chemical_equivalence`, `NMR_Interface`, `dihedral_fragments` and `atb_protein_ff`
(`fragment_merger` is only mentioned in its README; it has no `atb_outputs` import).

## Capabilities

* `atb_outputs.mol_data.MolData(initialiser)` -- builds atoms, bonds, rings and
  equivalence-group slots from a PDB string, a `chemistry_data_structure`
  `Molecule3D` (`_readMolecule3D`) or an `FDBMolecule`. `mol_data_from_mol_data_dict`
  rebuilds one from a stored dict (e.g. the cached `yml_moldata`).
  `MolData.unite_atoms()` produces the united-atom view.
* `atb_outputs.formats` -- one function per output, each taking a `MolData`:
  `pdb`, `g96`, `yml`, `template_yml`, `pickle`, `mol_data_dict`, `lgf`, `graph`,
  `ccd_cif`, `mol2` (via babel). The `yml` writer sanitises the data and dumps with
  libyaml's `CSafeDumper` (`yml.py`).
* `atb_outputs.itp.itp(mol_data, united=False)` -- GROMACS `.itp` text (not
  re-exported by `formats`). **Not the served writer**: nothing outside
  `testing/test.py` imports it; the `itp_allatom` the platform serves is
  `core/_outputs/itp.py`. Treat this one as an example (review 2026-10 C-01).
* `atb_outputs.graph` -- molecule graph rendering; needs the optional
  `graph_tool` package, which is **not installed** in the platform venv, so
  `formats.graph` returns an empty result (message on stderr). Unused (C-03).

## Usage

```python
from yaml import unsafe_load
from atb_outputs.formats import pdb, g96
from atb_outputs.itp import itp
from atb_outputs.mol_data import MolData, mol_data_from_mol_data_dict

def finalise(m):                      # writers read these two attributes
    m.var = {'REV_DATE': '', 'rnme': ''}
    m.completed = lambda x: False
    return m

m = finalise(MolData(open('mol.pdb').read()))           # from a PDB string
md = finalise(mol_data_from_mol_data_dict(unsafe_load(open('testing/data/21.yaml'))))
print(pdb(md, united=False)); print(g96(md)); print(itp(md))
# KNOWN BREAKAGE (review 2026-10 C-02): mol_data_from_mol_data_dict always raises
# MolDataFailure ('PDB string contains no atoms'), so the `md = ...` line above
# does not work. tests/conftest.py shows a working construction
# (MolData.__new__ + setattr of each dict key).

# from a Molecule3D (needs enough populated fields):
#   from chemistry_data_structure.parsing.input_parsers import GAMESS_to_Molecule3D
#   MolData(GAMESS_to_Molecule3D(open('qm.log').read()))
```

## How it fits the platform

`core.atb.outputs.Output` wraps these writers for every file the ATB serves;
`chemical_equivalence` computes the `equivalenceGroup` values on a `MolData`;
`website` uses it when ingesting submissions. Atom-order permutation for users'
structures is done in `core` (`Output.with_atom_order`), not here.

## Configuration

No env vars, no binaries. `mol2()` shells out to OpenBabel through
`chemistry_helpers.babel` (`atb_server_settings.babel_path`).
`formats.STORE_GRAPH_GT` / `graph.RAISE_IF_MISSING_GRAPH_TOOL` are module flags.

## Tests

`python -m pytest -p no:cacheprovider -q` from the repo root runs `tests/`
(golden-file characterisation of pdb/g96/yml/lgf/ccd_cif, `_readPDB` edge cases,
the exact `mol_data_dict` key set, strict xfails for C-02 and C-06; goldens are
regenerated with `tests/make_golden.py`). `testing/test.py` is the old print-only
smoke script and currently fails on its first line (C-02). The `Makefile` is a leftover (Python 3.5 paths, `git submodule` targets) and does not work.
