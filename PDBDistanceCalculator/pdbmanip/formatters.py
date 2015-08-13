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
    def write_listing(cls, structure, listings, output_file):
        with open(output_file, 'w+') as f:
            f.write(cls.format_listings(structure, listings))
    
    @classmethod
    @abc.abstractmethod
    def format_listings(cls, structure, listings):
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
            'resSeq': atom.parent.get_id()[1],
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
    def format_listings(cls, structure, listings):
        header = structure.header
        return {
            'head': header['head'],
            'title': header['name'],
            'author': header['author'],
            'listings': sorted([cls.format_NeighborListing(listing) for listing in listings], key=lambda listing: listing['atom']['serial'])
        }


class JSONNeighborFormatter(DictionaryFormatter):
    '''
    Neighbor Formatter using JSON format
    '''
    
    @classmethod
    def format_listings(cls, header, listings):
        #t = self.format_NeighborListing(self.neighborListings[0])
        #json.dumps(t)
        return json.dumps(super().format_listings(header, listings), sort_keys=True)
    
class PDBStyleNeighborFormatter(NeighborFormatter):
    '''
    Neighbor Formatter that mimics the style of a PDB file
    '''
    
    @classmethod
    def format_listings(cls, structure, listings):
        header = structure.header
        rows = []
        rows.append("HEADER    {0}".format(header["head"]))
        rows.append("TITLE    {0}".format(header["name"]))
        rows.append("AUTHOR    {0}".format(header["author"]))
        for listing in listings:
            def atom_info_format(atom):
                return "{name:<4} {resName:<3} {resSeq:<4} {chainID:<1} {modelID}".format(
                    name = atom.get_name(),
                    resName = atom.parent.get_resname(),
                    resSeq = atom.parent.get_id()[1],
                    chainID = atom.parent.parent.get_id(),
                    modelID = atom.parent.parent.parent.get_id()
                )
            
            def pdb_line_format(atom):
                pdb_format = "{record_name} {serial} {name} {altLoc} {resName} {chainID} {resSeq} {icode} {x} {y} {z} {occupancy} {bFactor} {element}"
                return pdb_format.format(
                    record_name = ("ATOM" if atom.parent.get_id()[0] == " " else "HETATM"),
                    serial = atom.get_serial_number(),
                    name = atom.get_fullname(),
                    altLoc = atom.get_altloc(),
                    resName = atom.parent.get_resname(),
                    chainID = atom.parent.parent.get_id(),
                    resSeq = atom.parent.get_id()[1],
                    icode = atom.parent.get_id()[2],
                    x = atom.get_coord()[0].item(),
                    y = atom.get_coord()[1].item(),
                    z = atom.get_coord()[2].item(),
                    occupancy = atom.get_occupancy(),
                    bFactor = atom.get_bfactor(),
                    element = atom.element
                )
             
            if not listing.neighbors:
                rows.append(atom_info_format(listing.atom) + " - No neighbors")
            else:
                for neighbor in listing.neighbors:
                    rows.append("{0}    {1}    {2}    {3}".format(
                        atom_info_format(listing.atom),
                        neighbor.distance,
                        atom_info_format(neighbor.atom),
                        pdb_line_format(neighbor.atom)
                    ))
        return "\n".join(rows)
    
    
    