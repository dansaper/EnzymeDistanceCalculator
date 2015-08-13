'''
Created on Jul 23, 2015

@author: dansaper
'''

from pdbmanip.retriever import PDBRetriever
from Bio.PDB.PDBParser import PDBParser
import argparse
import os
import re
import pdbmanip.searchers as searchers
import pdbmanip.formatters as formatters
import pdbmanip.filters.target_filters as target_filters
import pdbmanip.filters.compare_filters as compare_filters

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bio_dir",
                            default=os.path.join(os.getcwd(), "subunit_files"),
                            help="the directory that contains biological subunit files (default 'CWD/subunit_files')")
    parser.add_argument("--dist_dir",
                            default=os.path.join(os.getcwd(), "distance_files"),
                            help="the directory where the distance files will be created (default 'CWD/distance_files')")
    parser.add_argument("-f", "--input_file", help="""a file that provides a comma separated list of PDB codes.
                                                    If set, positional code arguments are ignored""")
    parser.add_argument("--json", action="store_true", help="Outputs in a JSON format, instead of psuedo-PDB format")
    parser.add_argument("-d", "--max_distance",
                            default='3.5',
                            help="The upper bound for the distance from checked atoms (default %(default)s)")
    parser.add_argument("codes", nargs=argparse.REMAINDER, help="any number of pdb code arguments - required if -f is not set")
    args = parser.parse_args()
    
    if (not args.input_file) and len(args.codes) == 0:
        parser.error("No PDB codes specified - see -h for details")
        
    
    return args

def generate_distance_files(subunits_dir, distance_dir, codes, max_distance, formatter):
    if not os.path.isdir(distance_dir):
        os.makedirs(distance_dir)
        
    retriever = PDBRetriever(subunits_dir)
    parser = PDBParser(QUIET=True)
    
    for code in codes:
        for file_path in retriever.get_subunit_files(code):
            structure_id = os.path.basename(file_path)
            structure = parser.get_structure(structure_id, file_path)
            searcher = searchers.KdTreeSearcher(
                list(structure.get_atoms()), 
                target_filters.FlavinFilter(),
                compare_filters.FlavinCompareFilter()
            )
            listings = searcher.search_for_atoms(float(max_distance))
            output_path = os.path.join(distance_dir, "{0}_dist_{1}".format(structure_id, max_distance))
            if os.path.isfile(output_path):
                print("Distance file already exists for {0}".format(structure_id))
            else:
                formatter.write_listing(structure, listings, output_path)


def flavin_funct():
    args = get_arguments()
    
    if args.input_file:
        with open(os.path.expanduser(args.input_file), 'r') as input_file:
            codes = re.findall('[a-zA-Z0-9_]{4}', input_file.read())
    elif args.codes:
        codes = args.codes
    else:
        print("")
        
    max_distance = args.max_distance
    subunits_dir = os.path.expanduser(args.bio_dir)
    distance_dir = os.path.expanduser(args.dist_dir)
    
    if args.json:
        formatter = formatters.JSONNeighborFormatter
    else:
        formatter = formatters.PDBStyleNeighborFormatter
    
    generate_distance_files(subunits_dir, distance_dir, codes, max_distance, formatter)
    