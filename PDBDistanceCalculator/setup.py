'''
@author: dansaper
'''
from setuptools import setup, find_packages

setup (
    name = "PDBDistanceCalculator",
    version = "0.1",
    description = "A program for finding the distance of atoms to isoalloxazine atoms",
    author = "Daniel Saper",
    packages = find_packages(),
    install_requires = [
        "numpy",
        "biopython"
    ],
    entry_points = {
        'console_scripts': [
            'flavinDistanceCalculator = pdbmanip.main:flavin_funct'
        ]
    }      
)