#from distutils.core import setup
from setuptools import setup, find_packages

setup(
    name='atb_outputs',
    version='1.0',
    author='Bertrand Caron',
    author_email='b.caron@uq.edu.au',
    url='https://atb.uq.edu.au',
    packages=find_packages(''),
    install_requires=['pyyaml', 'numpy'],
    description='Automated Topology Builder (ATB) molecule data object and output modules (PDB, YML, PICKLE, LGF, GRAPH) ',
    keywords='ATB PDB graph',
    python_requires='>=3.5',
)
