'''
@author: dansaper
'''

class NeighborListing:
    '''
    A holder for an atom and its neighbors
    '''
    
    def __init__(self, atom, neighbors):
        self.atom = atom
        self.neighbors = neighbors 


class DistListing:
    '''
    A holder for an atom with distance information
    '''
    def __init__(self, atom, distance):
        self.atom = atom
        self.distance = distance