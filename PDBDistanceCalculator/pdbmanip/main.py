'''
Created on Jul 23, 2015

@author: dansaper
'''

from pdbmanip.retriever import PDBRetriever
from Bio.PDB.PDBParser import PDBParser
import os
import pdbmanip.searchers as searchers
import pdbmanip.formatters as formatters
import pdbmanip.filters.target_filters as target_filters
import pdbmanip.filters.compare_filters as compare_filters

if __name__ == '__main__':
    codes = ['2dor']
    max_distance = '3'
    
    subunit_dir = os.path.join(os.getcwd(), "subunit_files")
    retriever = PDBRetriever(subunit_dir)
    
    parser = PDBParser()    
    distance_dir = os.path.join(os.getcwd(), "distance_files")
    if not os.path.isdir(distance_dir):
        os.makedirs(distance_dir)
            
    for code in codes:
        for file_path in retriever.get_subunit_files(code):
            structure_id = os.path.basename(file_path)
            structure = parser.get_structure(structure_id, file_path)
            searcher = searchers.KdTreeSearcher(
                list(structure.get_atoms()), 
                target_filters.FlavinFilter(),
                compare_filters.FlavinCompareFilter()
            )
            listings = searcher.search_for_atoms(int(max_distance))
            output_path = os.path.join(distance_dir, "{0}_dist_{1}".format(structure_id, max_distance))
            formatters.JSONNeighborFormatter.write_listing(listings, output_path)
