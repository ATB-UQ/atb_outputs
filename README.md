# Molecule Data (MolData) object

The main class is the `MolData` class defined in `atb_outputs.mol_data`.
It takes a `PDB` string as input and returns a `MolData` object.

The MolData object can then be fed to any of the functions in `atb_outputs.formats`:

* `pdb()`
* `yml()`
* `pickle()`
* `template_yml()`
* `graph()`
* `lgf()`
