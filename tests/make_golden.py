"""Regenerate tests/golden/* from the CURRENT writers (characterisation, not
specification). Run deliberately, review the diff."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from conftest import build_mol_data, mask_date, GOLDEN
from atb_outputs.formats import pdb, g96, lgf, ccd_cif

m = build_mol_data()
(GOLDEN / "toluene.pdb").write_text(mask_date(pdb(m)))
(GOLDEN / "toluene.g96").write_text(g96(m))
(GOLDEN / "toluene.lgf").write_text(lgf(m)[0][1])
(GOLDEN / "toluene.cif").write_text(ccd_cif(m, "TOL", "TOL"))
