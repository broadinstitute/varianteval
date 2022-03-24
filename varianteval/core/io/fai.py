from collections import namedtuple

Chr = namedtuple('Chr', 'name len')

class ChrFAIndex:
    def __init__(self):
        self.tid2chr = {}
        self.chr2tid = {}

    def add(self, tid, chr):
        self.tid2chr[tid] = chr
        self.chr2tid[chr.name] = tid

    def tids(self):
        return self.tid2chr.keys()

    def chr(self, tid):
        return self.tid2chr[tid]

    def chr_from_name(self, chr_name):
        return self.tid2chr[self.tid(chr_name)]

    def tid(self, chr_name):
        return self.chr2tid[chr_name]

    def contigs(self):
        return self.tid2chr.values()

    def has(self, chr_name):
        return chr_name in self.chr2tid

    def chr_names(self):
        return self.chr2tid.keys()


def load_faidx(fai_fname, primary=False):
    chr_index = ChrFAIndex()
    with open(fai_fname, "r") as faidx:
        for tid, line in enumerate(faidx):
            if primary and tid > 23:
                break  # only keep the autosomes
            # TODO: tid + 1 to correspond to chr names
            name, length, _, _, _ = line[:-1].split()
            chr_index.add(tid+1, Chr(name, int(length)))
    return chr_index
