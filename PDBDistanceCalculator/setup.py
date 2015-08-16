'''
@author: dansaper
'''
from cx_Freeze import setup, Executable

executables = [
    Executable('pdbmanip/main.py')
]

setup(name='flavinDistanceCalculator',
      version='0.1',
      description='test of flavinDistanceCalculator',
      executables=executables
)