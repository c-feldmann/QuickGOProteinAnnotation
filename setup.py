from setuptools import setup

setup(
    name='QuckGOProteinAnnotation',
    version='0.1',
    author='Christian Feldmann',
    license="BSD",
    packages=['go_protein_annotation', ],
    author_email='cfeldmann@bit.uni-bonn.de',
    description='Classes and functions useful for chemoinformatics',
    install_requires=['pandas', 'networkx', 'requests', 'tqdm']
)
