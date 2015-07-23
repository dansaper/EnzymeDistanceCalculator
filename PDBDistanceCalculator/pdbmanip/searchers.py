'''
@author: dansaper
'''

import abc
import math
from pdbmanip.structure_tuples import DistListing, NeighborListing

def coord_distance_sq(start, end):
    return math.fsum([(a - b)**2 for a, b in zip(start, end)])

class StructureSearcher(object):
    '''
    An abstract class for generating a list of distances from certain atoms in in a structure
    '''    
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_search_atoms(self):
        '''
        Should return a list of atoms in the search residues
            even those in different chains
        '''
        pass
    
    @abc.abstractmethod
    def get_search_residues(self):
        '''
        Should return a list of the requested residue names (typically passed in)
        '''
        pass
    
    @abc.abstractmethod
    def search_for_atoms(self, max_dist):
        '''
        Should return a list of NeighborListings, with one element for each of the Search Atoms
        '''
        pass
    
class ListSearcher(StructureSearcher):
    '''
    A StructureSearcher that uses an list for Basic O(n^2) searching
    '''
    
    def __init__(self, allAtoms, targetResidues):
        self.targetResidues = targetResidues[:]
        self.targetAtoms = [atom for atom in allAtoms 
                            if atom.parent.get_resname() in self.targetResidues]
        self.atoms = allAtoms

    def get_search_atoms(self):
        return self.targetAtoms
    
    def get_search_residues(self):
        return self.targetResidues
    
    def search_in_list(self, coords, max_dist):
        max_dist_sq = max_dist**2
        found = []
        for atom in self.atoms:
            distance_sq = coord_distance_sq(coords, atom.get_coord())
            if distance_sq <= max_dist_sq:
                found.append(DistListing(atom, math.sqrt(distance_sq)))
        return found
     
    def search_for_atoms(self, max_dist):
        return [NeighborListing(atom, self.search_in_list(atom.get_coord(), max_dist)) for atom in self.targetAtoms]


class KdNode:
    '''
    Private node class for kd-tree
    ''' 
    def __init__(self, axis, atom, left_child, right_child):
        self.axis = axis
        self.atom = atom
        self.left_child = left_child
        self.right_child = right_child
        
class KdTreeSearcher(StructureSearcher):
    '''
    A StructureSearcher that uses a kd-tree to lower the search space (see https://en.wikipedia.org/wiki/K-d_tree)
    '''
    
    def __init__(self, allAtoms, targetResidues):
        self.targetResidues = targetResidues[:]
        self.targetAtoms = []
        
        def kdTree(atoms, depth):
            if not atoms:
                return None
            
            axis = depth % 3
            
            atoms.sort(key=lambda atom: atom.get_coord()[axis])
            median = len(atoms) // 2
            
            #NEEDS PROFILING
            #While creating tree, we also pick out the atoms in residues we are interested in, 
            #    in order to save searching through each chain for the correct residues
            #    Speed benefit will probably depend on the percentage of chains that have each residue
            if atoms[median].parent.get_resname() in self.targetResidues:
                self.targetAtoms.append(atoms[median])
            
            return KdNode(
                axis = axis,
                atom = atoms[median],
                left_child = kdTree(atoms[:median], depth + 1),
                right_child = kdTree(atoms[median+1:], depth + 1)
            )
        
        self.tree = kdTree(allAtoms, 0)
        
    def get_search_atoms(self):
        return self.targetAtoms
    
    def get_search_residues(self):
        return self.targetResidues
    
    def search_in_tree(self, coords, max_dist):
        '''
        Returns a list of DistListing's within max_dist of coords
        '''
        max_dist_sq = max_dist**2
        def search_tree(node, within_range, coords):
            if node is not None:
                node_coords = node.atom.get_coord()
                axis = node.axis
                
                if math.fabs(coords[axis] - node_coords[axis]) <= max_dist:
                    distance_sq = coord_distance_sq(coords, node_coords)
                    if distance_sq <= max_dist_sq:
                        within_range.append(DistListing(node.atom, math.sqrt(distance_sq))) #store with actual distance
                    for tree in (node.left_child, node.right_child):
                        search_tree(tree, within_range, coords)
                else:
                    search_tree(
                        #Equal coords are always in range
                        node.left_child if coords[axis] < node_coords[axis] else node.right_child,
                        within_range, coords
                    )
            return within_range
        return search_tree(self.tree, [], coords)
    
    def search_for_atoms(self, max_dist):
        return [NeighborListing(atom, self.search_in_tree(atom.get_coord(), max_dist))
                for atom in self.targetAtoms]
