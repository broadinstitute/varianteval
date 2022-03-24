from intervaltree import IntervalTree, Interval
from varianteval.core.utils import NestedDict
import varianteval.core.constants as constants
from varianteval.core.io.repeat_tracks import repeat_track_iter

class GenomeContext:
    """
    Stores genome sequence annotation tracks (e.g. repeats)
    Internally each track is represented as an IntervalTree, indexed by track type and chromosome
    """
    def __init__(self, reference_genome):
        self.reference_genome = reference_genome
        self.tracks = NestedDict(NestedDict(IntervalTree))

    def load_track(self, repeat_track_fname, track_type):
        """Loads a genome track from the specified file"""
        if isinstance(track_type, constants.REPEAT_TYPE):
            for chr, start, end, repeat_name in repeat_track_iter(repeat_track_fname):
                self.tracks[track_type][chr].add(Interval(start, end, repeat_name))
        # TODO (loaders for other annotation types)

    def find_overlap(self, genome_interval, track_type):
        """ Return the set of track intervals overlapping the given genome interval """
        return self.tracks[track_type][genome_interval.chr_name].overlap(genome_interval.start, genome_interval.end)
