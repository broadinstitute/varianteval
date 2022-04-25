from bitarray import bitarray
from collections import defaultdict
import random
import numpy as np


class NestedDict(defaultdict):
    """Nested defaultdict"""
    def __call__(self):
        return NestedDict(self.default_factory)

def chr2num(chr_name):
    """
    Returns a numerical mapping for a primary chromosome name
    Assumes names are prefixed with "chr"
    """
    chr_id = chr_name[3:]
    if chr_id == 'X':
        return 23
    elif chr_id == 'Y':
        return 24
    elif chr_id == 'M' or chr_id == 'MT':
        return 25
    return int(chr_id)

def seq_to_num(seq):
    base2num = {'A': bitarray('00'), 'C': bitarray('01'), 'G': bitarray('10'), 'T': bitarray('11')}
    seq_bits = bitarray()
    seq_bits.encode(base2num, seq)
    return int(seq_bits.to01(), 2)

def shuffle_and_split(records, n_chunks):
    random.shuffle(records)
    return np.array_split(np.array(records), n_chunks)


def is_chromosome(string: str) -> bool:
    lower=string.lower()
    if lower[0].isdigit() or \
       lower == 'x' or lower == 'chrx' or \
       lower == 'y' or lower == 'chry' or \
       lower == 'm' or lower == 'chrm' or lower == 'mt' or lower == 'chrmt' or \
       (lower[0:3] == 'chr' and len(lower) <= 5):
       return True
     return False

