'''
Created on Jul 23, 2015

@author: dansaper
'''

from pdbmanip.retriever import PDBRetriever
from Bio.PDB.PDBParser import PDBParser
import pdbmanip.searchers as searchers
import pdbmanip.formatters as formatters
import pdbmanip.filters.target_filters as target_filters
import pdbmanip.filters.compare_filters as compare_filters

if __name__ == '__main__':
    retriever = PDBRetriever(['2dor'])
    file = retriever.files[0]
    parser = PDBParser()
    structure = parser.get_structure('2dor', file)
    searcher = searchers.KdTreeSearcher(
                                            list(structure.get_atoms()), 
                                            target_filters.FlavinFilter(),
                                            compare_filters.FlavinCompareFilter()
                                        )
    listings = searcher.search_for_atoms(3)
    formatter = formatters.JSONNeighborFormatter(listings)
    output_filename = "outputTree.json"
    formatter.write_listing(output_filename)
