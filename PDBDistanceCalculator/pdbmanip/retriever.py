'''
@author: dansaper
'''

from Bio.PDB import PDBList

class PDBRetriever:
    '''
    This class retrieves the requested PDB files
    '''

    def __init__(self, codes):
        '''
        Constructor
        '''
        pdbl = PDBList()
        #Download the requested PDB files from wwpdb
        self.files = [pdbl.retrieve_pdb_file(code) for code in codes]
        
        
            