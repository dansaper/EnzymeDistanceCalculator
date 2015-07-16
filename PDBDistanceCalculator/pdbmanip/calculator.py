'''
@author: dansaper
'''

import math

class DistListing:
    '''
    A holder for an atom with distance information
    '''
    def __init__(self, atom, distance):
        self.atom = atom
        self.distance = distance

class Node:
    '''
    Private node class for kd-tree
    ''' 
    def __init__(self, axis, atom, left_child, right_child):
        self.axis = axis
        self.atom = atom
        self.left_child = left_child
        self.right_child = right_child
    
class AtomKdTree:
    '''
    An encapsulation of a kd-tree for atoms (see https://en.wikipedia.org/wiki/K-d_tree)
    '''
    
    def __init__(self, allAtoms, targetResidues, targetAtomHolder):
        targetAtomHolder.clear()
        def kdTree(atoms, depth):
            if not atoms:
                return None
            
            axis = depth % 3
            
            atoms.sort(key=lambda atom: atom.get_coord()[axis])
            median = len(atoms) // 2
            
            #While creating tree, we also pick out the atoms in residues we are interested in, 
            #    in order to save an iteration
            if atoms[median].parent.get_resname() in targetResidues:
                targetAtomHolder.append(atoms[median])
            
            return Node(
                axis = axis,
                location = atoms[median],
                left_child = kdTree(atoms[:median], depth + 1),
                right_child = kdTree(atoms[median+1:], depth + 1)
            )
        
        self.tree = kdTree(allAtoms, 0)
    
    def search_within_range(self, coords, range_):
        range_sq = range_**2
        def search_tree(node, within_range, coords, range_):
            node_coords = node.atom.get_coord()
            axis = node.axix
            
            if math.fabs(coords[axis] - node_coords[axis]) <= range:
                distance_sq = math.fsum([(a - b)**2 for a, b in zip(coords, node_coords)])
                if distance_sq <= range_sq:
                    within_range.append(DistListing(node.atom, math.sqrt(distance_sq))) #store with actual distance
                for tree in (node.left_child, node.right_child):
                    search_tree(tree, within_range, coords, range_)
            else:
                search_tree(
                    node.left_child if coords[axis] < node_coords[axis] else node.right_child, #coords won't be equal
                    within_range, coords, range_
                )
            
            return within_range
    
        return search_tree(self.tree, [], coords, range_)

class NeighborListing:
    '''
    A holder for an atom and its neighbors
    '''
    
    def __init__(self, atom, neighbors):
        self.atom = atom
        self.neighbors = neighbors 

class DistanceCalculator:
    '''
    A datatype suitable for use in calculating distances within structure  (instead of Bio.PDB.PDBParser and Bio.PDB.Structure)
    '''

    def get_searcher_class(self):
        return AtomKdTree

    def __init__(self, header, targetResidues, allAtoms):
        '''
        targetResidues is the types of residues we actually care about
        allAtoms is a list of Bio.PDB.Atom objects
        '''
        
        self.header = header
        self.targetResidues = targetResidues
        self.targetAtoms = []
        self.searcher = self.get_searcher_class()(allAtoms, self.targetResidues, self.targetAtoms)
    
    def neighbor_search(self, max_range):
        for atom in self.targetAtoms:
            yield NeighborListing(atom, self.searcher.search_within_range(atom.get_coord(), max_range))
        
        