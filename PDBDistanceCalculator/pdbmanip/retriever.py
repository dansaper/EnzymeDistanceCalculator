'''
@author: dansaper
'''

from ftplib import FTP
import gzip
import urllib.request
import os

pdb_site = "ftp.wwpdb.org"
subunits_dir = "/pub/pdb/data/biounit/coordinates/divided/"

class PDBRetriever:
    '''
    This class retrieves the requested PDB files
    '''

    def __init__(self, codes, target_dir = os.path.join(os.getcwd(), "subunit_files")):
        '''
        We have to roll our own downloading, since Bio.PDB.PDBList only downloads structures, not biounits
        '''
        
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
            
        ftp = FTP(pdb_site)
        ftp.login()
        
        to_unpack = []
        
        for code in codes:
            code = code.lower()
            subdir = code[1:3] + "/"
            files = ftp.mlsd(subunits_dir + subdir)
            found_files = [f for f in files if f[0][0:4] == code]
            
            if not found_files:
                print("No biological assembly files found for {0}".format(code))
                return
            else:
                for f in found_files:
                    request = urllib.request.urlopen("ftp://" + pdb_site + subunits_dir)
                    data = request.data
                    file_path = os.path.join(target_dir, f[0])
                    with open(file_path, "x+b") as target_file:
                        target_file.write(data)
                        to_unpack.append(file_path)
        
        unzipped_files = []
        for zipped in to_unpack:
            new_file_path = os.path.splitext(zipped)
            #if it already exists, we don't need to unzip it,
            #    but we pretend we did for the rest of the program
            unzipped_files.append(new_file_path)
            
            with open(new_file_path, "x+") as new_file:
                with gzip.open(zipped, "rt") as unzipped:
                    new_file.write(unzipped.read())
        
        self.files = unzipped_files
        
        
            