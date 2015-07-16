'''
@author: dansaper
'''

from Bio.PDB import PDBList, PDBParser

class PDBRetriever:
    '''
    This class retrieves the requested PDB files and parses them
    '''

    def __init__(self, codes):
        '''
        Constructor
        '''
        pdbl = PDBList()
        #Download the requested PDB files from wwpdb
        self.files = [pdbl.retrieve_pdb_file(code) for code in codes]
        
        
            