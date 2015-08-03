'''
@author: dansaper
'''

isoalloxazine_atoms = {'N1', 'C2', 'O2',
                                    'N3', 'C4', 'O4',
                                    'C4A', 'N5', 'C5A',
                                    'C6', 'C7', 'C7M',
                                    'C8', 'C8M', 'C9',
                                    'C9A', 'N10', 'C10'}

def on_flavin(atom):
    '''
    Checks if atom is on a flavin
    '''
    return atom.parent.get_resname() in {'FMN', 'FAD'}

def on_same_residue(atom_a, atom_b):
    a_id = atom_a.get_full_id()
    b_id = atom_b.get_full_id()
    
    #[3][1] == resSeq, [2] == chain id
    return (a_id[3][1] == b_id[3][1]) and (a_id[2] == b_id[2])

def in_isoalloxazine(atom):
    return on_flavin(atom) and atom.get_name() in isoalloxazine_atoms
