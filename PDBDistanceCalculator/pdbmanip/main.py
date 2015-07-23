'''
Created on Jul 23, 2015

@author: dansaper
'''

from pdbmanip.retriever import PDBRetriever
from Bio.PDB.PDBParser import PDBParser
from pdbmanip.searchers import ListSearcher, KdTreeSearcher
from pdbmanip.formatters import JSONNeighborFormatter, DictionaryFormatter

if __name__ == '__main__':
    retriever = PDBRetriever(['2dor'])
    file = retriever.files[0]
    print(retriever.files)
    parser = PDBParser()
    structure = parser.get_structure('2dor', file)
    searcher = KdTreeSearcher(list(structure.get_atoms()), ["FMN"])
    listings = searcher.search_for_atoms(3)
    formatter = JSONNeighborFormatter(listings)
    output_filename = "outputTree.json"
    formatter.write_listing(output_filename)
    