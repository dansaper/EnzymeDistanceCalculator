'''
@author: dansaper
'''

import abc
import json
import numpy as np

class NeighborFormatter(object):
    '''    
    abstract interface for formatting neighbor listings
    '''
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, neighborListings):
        self.neighborListings = neighborListings

    def write_listing(self, filename):
        with open(filename, 'w+') as f:
            f.write(self.format_listings())
    
    @abc.abstractmethod
    def format_listings(self):
        pass
    
    
class DictionaryFormatter(NeighborFormatter):
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
            'coords': {#We do need to convert from numpy type to normal in order to use JSON
                'x': atom.get_coord()[0].item(),
                'y': atom.get_coord()[1].item(),
                'z': atom.get_coord()[2].item()
            },
            'occupancy': atom.get_occupancy()    
        }
    
    def format_NeighborListing(self, listing):
        return {
            'atom': self.format_atom(listing.atom),
            'neighbors': sorted([self.format_DistListing(neighbor) for neighbor in listing.neighbors], key=lambda listing: listing['atom']['serial'])
        }
        
    def format_listings(self):
        return sorted([self.format_NeighborListing(listing) for listing in self.neighborListings], key=lambda listing: listing['atom']['serial'])


class JSONNeighborFormatter(DictionaryFormatter):
    '''
    Neighbor Formatter using JSON format
    '''
    def format_listings(self):
        #t = self.format_NeighborListing(self.neighborListings[0])
        #json.dumps(t)
        return json.dumps(DictionaryFormatter.format_listings(self), sort_keys=True)    
        