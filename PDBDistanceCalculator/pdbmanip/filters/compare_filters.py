'''
@author: dansaper
'''

import abc
import pdbmanip.flavin_info as fl

class AtomCompareFilter(object):
    '''
    Abstract base class for checking if atoms should be compared
    '''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def should_compare(self, atom_a, atom_b):
        '''
        Checks if atom_a and atom_b are identical
        '''
        return atom_a.get_full_id() != atom_b.get_full_id()

class FlavinCompareFilter(AtomCompareFilter):
    '''
    We don't compare atoms on the same isoalloxazine
    '''
    def should_compare(self, atom_a, atom_b):
        return not (fl.on_same_residue(atom_a, atom_b)
                and fl.in_isoalloxazine(atom_a)
                and fl.in_isoalloxazine(atom_b)
            ) and super().should_compare(atom_a, atom_b)