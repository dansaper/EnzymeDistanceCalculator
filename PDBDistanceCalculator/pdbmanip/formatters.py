'''
@author: dansaper
'''

class JSONNeighborFormatter(object):
    '''
    abstract interface for formatting neightobr listings
    '''
    def format_DistListing(self, distListing):
        return {
            'atom': self.format_atom(distListing.atom),
            'distance': distListing.distance
        }
    
    def format_atom(self, atom):
        return {
            'kind': 'ATOM',
            'serial': atom.get_serial_number(),
            'name': atom.get_name(),
            'altLoc': atom.get_altloc(),
            'resName': atom.parent.get_resname(),
            'chainID': atom.parent.parent.get_id(), #ask about this
            'resSeq': atom.parent.get_id()[1], #ask about this
            'coords': {
                'x': atom.get_coord()[0],
                'y': atom.get_coord()[1],
                'z': atom.get_coord()[2]
            },
            'occupancy': atom.get_occupancy()    
        }
    
    def format_NeighborListing(self, listing):
        return {
            'atom': self.format_atom(listing.atom),
            'neighbors': [self.format_DistListing(neighbor) for neighbor in listing.neighbors]
        }

    def write_listing(self, file):

    def __init__(self, params):
        '''
        Constructor
        '''
        