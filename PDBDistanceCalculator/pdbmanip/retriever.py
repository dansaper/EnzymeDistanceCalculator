'''
@author: dansaper
'''

from ftplib import FTP
import gzip
import urllib
import os

pdb_site = "ftp.wwpdb.org"
subunits_dir = "/pub/pdb/data/biounit/coordinates/divided/"

class PDBRetriever:
    '''
    This class retrieves the requested PDB files
    '''

    def __init__(self, target_dir):
        '''
        We have to roll our own downloading, since Bio.PDB.PDBList only downloads structures, not biounits
        '''
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)
            
        self.target_dir = target_dir
        self.subunit_files_dict = {}

    def get_subunit_files(self, code):
        if code not in self.subunit_files_dict:
            self.subunit_files_dict[code] = self.download_subunit_files(code)
        return self.subunit_files_dict[code]
        
    def find_subunit_urls(self, code):
        ftp = FTP(pdb_site)
        ftp.login()
        subdir = subunits_dir + code[1:3] + "/"
        found_files = [f for f in ftp.mlsd(subdir) if f[0][0:4] == code]
        return ["ftp://" + pdb_site + subdir + f[0] for f in found_files]
    
    def download_subunit_files(self, code):
        return [self.download_subunit_file(url) for url in self.find_subunit_urls(code)]
        
    def download_subunit_file(self, url):
        #parse for filename and remove .gz extension
        new_file_name = os.path.splitext(urllib.parse.urlparse(url).path.split("/")[-1])[0]
        new_file_path = os.path.join(self.target_dir, new_file_name)
        
        #if file does not already exist, download it
        if not os.path.isfile(new_file_path):
            with urllib.request.urlopen(url) as response:
                data = response.read()
                with open(new_file_path, "x+b") as target_file:
                    target_file.write(gzip.decompress(data))
        return new_file_path