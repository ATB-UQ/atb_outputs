"""Fixtures: a MolData for toluene (testing/data/21.yaml), built by hand so the
writers' goldens do not depend on ``mol_data_from_mol_data_dict`` (review 2026-10
C-02, fixed in wave 2 and tested on its own)."""
import re
from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from atb_outputs.mol_data import MolData

HERE = Path(__file__).parent
GOLDEN = HERE / "golden"
YAML_21 = HERE.parent / "testing" / "data" / "21.yaml"


def build_mol_data():
    d = yaml.unsafe_load(YAML_21.read_text())
    m = MolData.__new__(MolData)  # bypass __init__ (no PDB to parse)
    for k, v in d.items():
        setattr(m, k, deepcopy(v))
    m.var = dict(m.var)
    m.var.setdefault("rnme", "MOL")
    m.var["total_charge"] = 0
    m.completed = lambda key: True
    return m


@pytest.fixture
def mol_data():
    return build_mol_data()


def mask_date(text):
    """The PDB HEADER carries today's date (dd-Mon-yy)."""
    return re.sub(r"\d{2}-[A-Za-z]{3}-\d{2}", "DD-Mon-YY", text)


def golden(name):
    return (GOLDEN / name).read_text()
