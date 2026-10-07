# atb_outputs

The ATB molecule data object (`MolData`) and the writers that turn it into output
files (PDB, G96, YML, pickle, LGF, graph image, CCD CIF, MOL2, GROMACS ITP).
Pure Python, installable on its own (`pip install -e .`), but note that
`setup.cfg` declares **no dependencies**: it needs `numpy`, `pyyaml` and, for
`mol2()`, the `chemistry_helpers` sibling (see Configuration).

**Status:** live-support. Imported by `core` (the topology generator), `website`,
`chemical_equivalence`, `NMR_Interface`, `fragment_merger`, `dihedral_fragments`
and `atb_protein_ff`.

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
  re-exported by `formats`). Improper force constants are converted with `(180/pi)^2`.
* `atb_outputs.graph` -- molecule graph rendering; needs the optional
  `graph_tool` package and degrades to an empty result (message on stderr)
  without it.

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

`testing/test.py` is a print-only smoke script (run from `testing/`:
`python test.py`), not a pytest suite -- there are no assertions. The
`Makefile` is a leftover (Python 3.5 paths, `git submodule` targets) and does not work.
