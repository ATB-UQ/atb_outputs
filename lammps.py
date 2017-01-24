
from atb_outputs.atb_lammps.lib import Moltemplate 
from atb_outputs.atb_lammps.lib.LammpsGromosForcefield import LammpsGromosForcefield as Forcefield
from atb_outputs.helpers.types_helpers import Dict, Output_File
from io import StringIO
from os import path
import yaml

SCALE14 = 0.5

def moltemplate(
        mol_data: Dict,
        optimized: bool,
        united: bool,
        header: str,
        ) -> Output_File:

    # use residue name as molecule identifier
    name = mol_data["var"]["rnme"]

    # the name of the forcefield is used in the moltemplate namespace
    ffname = "GROMOS_{}_ATB".format(mol_data["var"]["ff_version"])

    forcefield = None # A forcefield object can optionally be supplied to
                      # cross-check atom names and masses

    pdb_precision = False # Only used for testing to match precision

    output_buffer = StringIO() # IO buffer stores the output string
    output_buffer.write(header)

    # Write instructions to the buffer
    instructions_file = "{}/atb_lammps/molecule_instructions.md".format(
        path.dirname(path.abspath(__file__))
        )
    with open(instructions_file) as readme:
        instructions = Moltemplate.make_header(readme, "#| ",
                ffname, name,
                filename(name, united, optimized))
    output_buffer.write(instructions)

    # Generate the molecule data and write it to the buffer
    Moltemplate.write_molecule(
        output_buffer,
        mol_data,
        name,
        forcefield,
        ffname,
        united,
        optimized,
        pdb_precision,
        check_forcefield = False,
        )
    # Note: check_forcefield = True requires the force field information to cross-
    # check atom types and masses. This check has never found a problem during
    # tests

    # return the contents of the buffer
    return output_buffer.getvalue()

def filename(resname, united, optimized):
    return "{}_{}atom_{}_geometry.lt".format(
            resname,
            "uni" if united else "all",
            "optimized" if optimized else "original")

def ifp2moltemplate(ifp: str, forcefield_name: str):

    output_buffer = StringIO() # IO buffer stores the output string

    atblammpspath = path.dirname(path.abspath(__file__))+"/atb_lammps"
    
    with open("{}/masses.yml".format(atblammpspath)) as f:
        masses = yaml.load(f)

    notice = """
####
# For this LAMMPS moltemplate version of the forcefield, the Lennard-Jones
# 1-4 interactions are uniformly scaled by 0.5 instead of using the usual
# GROMOS 1-4 Lennard-Jones parameters.
####
"""
    output_buffer.write(notice)

    instructions_file = "{}/molecule_instructions.md".format(atblammpspath)
    with open(instructions_file) as readme:
        instructions = Moltemplate.make_header(readme, "#| ",
                "GROMOS_{}_ATB".format(forcefield_name),
                "my_molecule",
                "my_molecule")
    output_buffer.write(instructions)

    ifp_buffer = StringIO(ifp)

    forcefield = Forcefield(StringIO(ifp), masses, SCALE14)
    Moltemplate.write_forcefield(output_buffer, forcefield_name, forcefield)

    return output_buffer.getvalue()
