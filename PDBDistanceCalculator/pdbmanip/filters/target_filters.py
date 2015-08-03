'''
@author: dansaper
'''

import abc
import pdbmanip.flavin_info as fl


class AtomFilter(object):
    '''
    Abstract base class for filtering atoms for use in search
    '''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def fits_filter(self, atom, target_list):
        '''
        Checks if atom should be added to target_list
            Does not modify target_list
        '''
        pass



class FlavinFilter(AtomFilter):
    '''
        Filters for isoalloxazine atoms
    '''

    def fits_filter(self, atom):
        return fl.in_isoalloxazine(atom)