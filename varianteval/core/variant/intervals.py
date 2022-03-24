from collections import defaultdict
from intervaltree import IntervalTree
from varianteval.core.utils import chr2num

class GenomeInterval(tuple):
    """
    Represents an interval on the genome defined by the 3-tuple: chromosome name, start, and end positions
    Assumes a 0-based start, includes both endpoints
    """
    def __new__(cls, chr_name, start, end):
        return tuple.__new__(GenomeInterval, (chr_name, start, end))

    def __init__(self, chr_name, start, end):
        self.chr_name = chr_name
        self.start = start
        self.end = end

    def __len__(self):
        """ Returns the length of the interval """
        return self.end - self.start + 1

    def __str__(self):
        """ Returns a string representation of the interval """
        return "%s:%d-%d" % (self.chr_name, self.start, self.end)

    def __lt__(self, interval):
        """ Less-than interval comparator """
        if self.chr_name == interval.chr_name:
            return self.start < interval.start
        return chr2num(self.chr_name) < chr2num(interval.chr_name)

    @staticmethod
    def pos_to_interval(chr_name, pos):
        """ Converts a single position to an interval """
        return GenomeInterval(chr_name, pos, pos)

    def overlap(self, interval):
        """ Returns the size of the overlap between the two intervals """
        return max(0, min(interval.end, self.end) - max(interval.start, self.start) + 1)

class GenomeIntervalPair:
    """
    Represents a pair of genome intervals
    """
    def __init__(self, interval_a, interval_b):
        self.interval_a = interval_a
        self.interval_b = interval_b

    def __str__(self):
        return "%s_&_%s" % (str(self.interval_a), str(self.interval_b))


class GenomeIntervalTree:
    """
    Stores genome intervals as a set of interval trees indexed by chromosome name
    Supports range overlap queries
    """
    def __init__(self, genome_intervals):
        self.chr2tree = defaultdict(IntervalTree)
        for interval in genome_intervals:
            self.add(interval)

    def add(self, interval):
        self.chr2tree[interval.chr_name].addi(interval.start, interval.end, interval)

    def get_overlap(self, interval, overlap_frac=0.2):
        """
        Returns the list of overlapped intervals in the tree s.t. their overlap is >= overlap_frac
        Interval overlap is computed as IoM (intersection/min)
        """
        if not self.chr2tree[interval.chr_name].overlaps(interval.start, interval.end):
            return []
        candidates = self.chr2tree[interval.chr_name].overlap(interval.start, interval.end)
        output_intervals = []
        for c in candidates:
            candidate_interval = c.data
            overlap = candidate_interval.overlap(interval)
            if float(overlap/min(len(interval), len(candidate_interval))) >= overlap_frac:
                output_intervals.append(candidate_interval)
        return output_intervals

    def overlaps(self, interval, overlap_frac=0.2):
        """
        Returns true if the given interval overlaps at least one interval in the tree
        s.t. their overlap is >= overlap_frac
        Interval overlap fraction is computed as intersection/min
        """
        # if no interval overlaps at all
        if not self.chr2tree[interval.chr_name].overlaps(interval.start, interval.end):
            return False
        if not overlap_frac:
            return True
        return len(self.get_overlap(interval, overlap_frac))
