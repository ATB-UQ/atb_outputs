from chemistry_helpers.babel import babel_output


def use_babel(pdb):
    return babel_output(pdb, in_format='pdb', out_format='mol2', babel_executable="/Users/uqmstroe/miniconda3/bin/obabel")
