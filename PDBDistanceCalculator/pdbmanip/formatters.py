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

    @classmethod
    def write_listing(cls, listings, output_file):
        with open(output_file, 'w+') as f:
            f.write(cls.format_listings(listings))
    
    @classmethod
    @abc.abstractmethod
    def format_listings(cls, listings):
        pass
    
    
class DictionaryFormatter(NeighborFormatter):
    '''
    Converts listing to a dictionary of with plain strings and numbers
    '''
    @classmethod
    def format_DistListing(cls, distListing):
        return {
            'atom': cls.format_atom(distListing.atom),
            'distance': distListing.distance
        }
    
    @classmethod
    def format_atom(cls, atom):
        return {
            'kind': 'ATOM',
            'serial': atom.get_serial_number(),
            'name': atom.get_name(),
            'altLoc': atom.get_altloc(),
            'resName': atom.parent.get_resname(),
            'modelID': atom.parent.parent.parent.get_id(),
            'chainID': atom.parent.parent.get_id(),
            'resSeq': atom.parent.get_id()[1], #ask about this
            'coords': {#We do need to convert from numpy type to normal in order to use JSON
                'x': atom.get_coord()[0].item(),
                'y': atom.get_coord()[1].item(),
                'z': atom.get_coord()[2].item()
            },
            'occupancy': atom.get_occupancy()    
        }
    
    @classmethod
    def format_NeighborListing(cls, listing):
        return {
            'atom': cls.format_atom(listing.atom),
            'neighbors': sorted([cls.format_DistListing(neighbor) for neighbor in listing.neighbors], key=lambda listing: listing['atom']['serial'])
        }
    
    @classmethod
    def format_listings(cls, listings):
        return sorted([cls.format_NeighborListing(listing) for listing in listings], key=lambda listing: listing['atom']['serial'])


class JSONNeighborFormatter(DictionaryFormatter):
    '''
    Neighbor Formatter using JSON format
    '''
    
    @classmethod
    def format_listings(cls, listings):
        #t = self.format_NeighborListing(self.neighborListings[0])
        #json.dumps(t)
        return json.dumps(super().format_listings(listings), sort_keys=True)
    
class PDBStyleNeighborFormatter(DictionaryFormatter):
    '''
    Neighbor Formatter that mimics the style of a PDB file
    '''