"""
The main data structures for representing variants.

.. image:: ../../../docs/sv.png

Remark: phasing is not currently modeled and it is left to future work.
"""

import copy
import math
import pickle

from varianteval.core.constants import *
import varianteval.core.utils as utils
import varianteval.core.variant.intervals as intervals


class Track:
    """
    Interval annotation of a sequence
    """
    def __init__(self):
        self.name: str = None

        # List of ``[first..last]`` coordinates, inclusive, zero-based.
        self.intervals: list[intervals.GenomeInterval] = []

    def __eq__(self, other):
        if isinstance(other, Track): return self.name == other.name and self.intervals == other.intervals
        return False

    def __hash__(self):
        return hash(self.name) ^ hash(tuple(self.intervals))

    def __str__(self):
        return ("%s \n" % self.name), "\n".join([str(i) for i in self.intervals])


class Sequence:
    """
    A distinct chromosome, contig, or inserted sequence. Only sequences with some breakpoint should be loaded in memory.
    """
    def __init__(self):
        # General properties
        self.name: str = None  # Can be None if it is not a known sequence.
        self.is_circular: bool = False
        self.sequence: str = None  # The actual nucleotides. Can be None if not needed.
        self.tracks: list[Track] = []

        # The length of the sequence may be uncertain (e.g. if it is an insertion).
        self.length_min: int = 0
        self.length_max: int = 0
        self.length_avg: int = 0
        self.length_std: float = 0.0
        self.length_probability_function: int = PROBABILITY_FUNCTION.UNKNOWN  # ID of a function

    def __eq__(self, other):
        if isinstance(other, Sequence):
            return self.name == other.name and \
                   self.is_circular == other.is_circular and \
                   self.length_min == other.length_min and \
                   self.length_max == other.length_max and \
                   self.length_avg == other.length_avg and \
                   self.length_std == other.length_std and \
                   self.length_probability_function == other.length_probability_function and \
                   self.sequence == other.sequence
        return False

    def __hash__(self):
        return hash(self.name) ^ hash(self.is_circular) ^ hash(self.sequence) ^ hash(self.length_min) ^ \
               hash(self.length_max) ^ hash(self.length_avg) ^ hash(self.length_std) ^ \
               hash(self.length_probability_function)

    def __lt__(self, other):
        if self.name < other.name: return True
        elif self.name > other.name: return False
        if self.sequence < other.sequence: return True
        elif self.sequence > other.sequence: return False
        return False

    def __le__(self, other):
        if self.name < other.name: return True
        elif self.name > other.name: return False
        if self.sequence < other.sequence: return True
        elif self.sequence > other.sequence: return False
        return True

    def __gt__(self, other):
        if self.name < other.name: return False
        elif self.name > other.name: return True
        if self.sequence < other.sequence: return False
        elif self.sequence > other.sequence: return True
        return False

    def __ge__(self, other):
        if self.name < other.name: return False
        elif self.name > other.name: return True
        if self.sequence < other.sequence: return False
        elif self.sequence > other.sequence: return True
        return True

    def __str__(self):
        return ("%s(%s) %s[%d..(%d,%f)..%d]\n" % (
            self.name, self.is_circular, self.length_probability_function, self.length_min, self.length_avg,
            self.length_std, self.length_max)), "Tracks: \n", ("\n".join([str(i) for i in self.tracks]))

    def clear(self):
        self.name = None
        self.is_circular = False
        self.sequence = None
        self.tracks = None
        self.length_min = 0
        self.length_max = 0
        self.length_avg = 0
        self.length_std = 0
        self.length_probability_function = PROBABILITY_FUNCTION.UNKNOWN

    def get_type(self) -> int:
        """
        :return: 0=chromosome, 1=contig, 2=inserted string.
        """
        if utils.is_chromosome(self.name): return 0
        elif self.name == VCF_INSERTION_STRING_NAME: return 2
        return 1


class Breakpoint:
    """
    Assume that every position of the reference is a strip of paper with a left side and a right side. A breakpoint is
    an uncertain position with certain side.
    """
    def __init__(self):
        self.sequence: Sequence = None

        # If we cut a sequence of length ``n`` at position ``p`` (zero-based), we need to know whether we are using
        # ``[0..p]`` (true) or ``[p..n-1]`` (false) in an adjacency. These are defined in ``constants.py``.
        self.side: int = BREAKPOINT_SIDE.UNKNOWN

        # To model uncertainty, a breakpoint position is a probability distribution over an interval (the interval has
        # size one iff there is no uncertainty).
        #
        # Remark: some callers might not return an interval, but just a probability distribution. In this case we set
        # ``position_first = position_last = -1``.
        #
        # Remark: if ``sequence`` is circular, then it can happen that ``position_first > position_last``.
        self.position_first: int = -1  # Inclusive, zero-based.
        self.position_last: int = -1  # Inclusive, zero-based.
        self.position_avg: int = -1
        self.position_std: float = 0.0
        self.position_probability_function: int = PROBABILITY_FUNCTION.UNKNOWN  # Defined in ``constants.py``.

        # One count for each individual. The count tells how many haplotypes in the individual (0,1,2,...) use an
        # adjacency that involves the breakpoint. Use -1 for unknown.
        #
        # Remark: we give a genotype to a breakpoint, as well, since the same breakpoint might be involved in multiple,
        # distinct events in the haplotypes of the same individual.
        self.genotype: list[int] = []

    def __eq__(self, other):
        if isinstance(other, Breakpoint):
            return self.sequence == other.sequence and \
                   self.side == other.side and \
                   self.position_first == other.position_first and \
                   self.position_last == other.position_last and \
                   self.position_avg == other.position_avg and \
                   self.position_std == other.position_std and \
                   self.position_probability_function == other.position_probability_function
        return False

    def __hash__(self):
        return hash(self.sequence) ^ hash(self.side) ^ hash(self.position_first) ^ hash(self.position_last) ^ \
               hash(self.position_avg) ^ hash(self.position_std) ^ hash(self.position_probability_function)

    def __lt__(self, other):
        if self.sequence < other.sequence: return True
        elif self.sequence > other.sequence: return False
        if self.position_first < other.position_first: return True
        elif self.position_first > other.position_first: return False
        if self.position_last < other.position_last: return True
        elif self.position_last > other.position_last: return False
        if self.side < other.side: return True
        elif self.side > other.side: return False
        return False

    def __le__(self, other):
        if self.sequence < other.sequence: return True
        elif self.sequence > other.sequence: return False
        if self.position_first < other.position_first: return True
        elif self.position_first > other.position_first: return False
        if self.position_last < other.position_last: return True
        elif self.position_last > other.position_last: return False
        if self.side < other.side: return True
        elif self.side > other.side: return False
        return True

    def __gt__(self, other):
        if self.sequence < other.sequence: return False
        elif self.sequence > other.sequence: return True
        if self.position_first < other.position_first: return False
        elif self.position_first > other.position_first: return True
        if self.position_last < other.position_last: return False
        elif self.position_last > other.position_last: return True
        if self.side < other.side: return False
        elif self.side > other.side: return True
        return False

    def __ge__(self, other):
        if self.sequence < other.sequence: return False
        elif self.sequence > other.sequence: return True
        if self.position_first < other.position_first: return False
        elif self.position_first > other.position_first: return True
        if self.position_last < other.position_last: return False
        elif self.position_last > other.position_last: return True
        if self.side < other.side: return False
        elif self.side > other.side: return True
        return True

    def __str__(self):
        if self.side == BREAKPOINT_SIDE.LEFT:
            self.left = "<"
            self.right = "]"
        else:
            self.left = "["
            self.right = ">"
        return ("%s %s%s%d..(%d,%f)..%d%s" % (self.sequence.name, self.position_probability_function, self.left,
                self.position_first, self.position_avg, self.position_std, self.position_last, self.right)), \
            "Genotype: ", (" ".join([str(i) for i in self.genotype]))

    def clear(self):
        self.sequence = None
        self.side = BREAKPOINT_SIDE.UNKNOWN
        self.position_first = -1
        self.position_last = -1
        self.position_avg = -1
        self.position_std = -1
        self.position_probability_function = PROBABILITY_FUNCTION.UNKNOWN

    def is_similar(self, other, identity_threshold: int, jaccard_threshold: float) -> bool:
        """
        :other: of type Breakpoint;
        :return: TRUE iff the two breakpoints are intervals of size one and are at most ``identity_threshold`` bps
        apart, or if at least one breakpoint is an interval of size greater than one and the two intervals have Jaccard
        similarity at least ``jaccard_threshold``.
        """
        if self.sequence != other.sequence or self.side != other.side: return False
        if self.position_first == self.position_last:
            if other.position_first == other.position_last:
                return abs(self.position_avg - other.position_avg) <= identity_threshold
            else:
                jaccard = float(min(other.position_last, self.position_last) - max(other.position_first, self.position_first) + 1) / (max(other.position_last, self.position_last) - min(other.position_first, self.position_first) + 1)
                return jaccard >= jaccard_threshold
        else:
            jaccard = float(min(other.position_last, self.position_last) - max(other.position_first, self.position_first) + 1) / (max(other.position_last, self.position_last) - min(other.position_first, self.position_first) + 1)
            return jaccard >= jaccard_threshold

    def update_position(self, update_avg: bool):
        """
        Given ``position_{first,last}``, resets all other ``position_*`` fields.

        :param update_avg: if FALSE, ``position_avg`` is not updated.
        """
        if self.position_first == self.position_last:
            self.position_probability_function = PROBABILITY_FUNCTION.DELTA
            if update_avg: self.position_avg = self.position_first
            self.position_std = 0
        else:
            if self.position_probability_function in [PROBABILITY_FUNCTION.DELTA, PROBABILITY_FUNCTION.UNIFORM]:
                self.position_probability_function = PROBABILITY_FUNCTION.UNIFORM
                if update_avg: self.position_avg = (self.position_first + self.position_last) // 2
                self.position_std = math.sqrt(float((self.position_last - self.position_first + 1) ** 2) / 12)
            else:
                pass  # We keep the STD and AVG of the original function

    def move_right(self, bps: int):
        self.position_avg = min(self.position_avg + bps, self.sequence.length_max - 1)
        self.position_first = min(self.position_first + bps, self.sequence.length_max - 1)
        self.position_last = min(self.position_last + bps, self.sequence.length_max - 1)
        self.update_position(False)

    def move_left(self, bps: int):
        self.position_avg = max(self.position_avg - bps, 0)
        self.position_first = max(self.position_first - bps, 0)
        self.position_last = max(self.position_last - bps, 0)
        self.update_position(False)

    def relax_by_radius(self, radius: int):
        """
        Ensures that ``[position_avg-radius..position_avg+radius]`` belongs to the uncertainty interval.
        """
        self.position_first = max(0, min(self.position_avg - radius, self.position_first))
        self.position_last = min(max(self.position_avg + radius, self.position_last), self.sequence.length_max)
        if self.position_probability_function == PROBABILITY_FUNCTION.DELTA:
            self.position_probability_function = PROBABILITY_FUNCTION.UNIFORM
            self.position_std = math.sqrt(float((self.position_last - self.position_first + 1) ** 2) / 12)
        elif self.position_probability_function == PROBABILITY_FUNCTION.UNIFORM:
            self.position_std = math.sqrt(float((self.position_last - self.position_first + 1) ** 2) / 12)
        else:
            pass  # We keep the STD of the original function


class Reference_Interval:
    """
    An interval of the reference, with anomalous copy number.

    This class models the case in which the caller has evidence for a clear change in coverage, but it doesn't have
    evidence for new context.
    """
    def __init__(self, breakpoint1: Breakpoint, breakpoint2: Breakpoint, copy_number: float):
        # Cannot be None.
        self.breakpoint1: Breakpoint = breakpoint1
        self.breakpoint2: Breakpoint = breakpoint2
        self.copy_number: float = copy_number

        # One count for each individual. The count tells how many haplotypes in the individual (0,1,2,...) contain this
        # interval with this specific copy number. Use -1 for unknown.
        #
        # Remark: we give a genotype to an interval, as well, just for generality. I.e. we allow the same interval to be
        # involved in multiple, distinct events in the haplotypes of the same individual.
        self.genotype: list[int] = []

    def __eq__(self, other):
        if isinstance(other, Reference_Interval):
            return self.breakpoint1 == other.breakpoint1 and \
                   self.breakpoint2 == other.breakpoint2 and \
                   self.copy_number == other.copy_number
        return False

    def __hash__(self):
        return hash(self.breakpoint1) ^ hash(self.breakpoint2) ^ hash(self.copy_number)

    def __lt__(self, other):
        if self.breakpoint1 < other.breakpoint1: return True
        elif self.breakpoint1 > other.breakpoint1: return False
        if self.breakpoint2 < other.breakpoint2: return True
        elif self.breakpoint2 > other.breakpoint2: return False
        return False

    def __le__(self, other):
        if self.breakpoint1 < other.breakpoint1: return True
        elif self.breakpoint1 > other.breakpoint1: return False
        if self.breakpoint2 < other.breakpoint2: return True
        elif self.breakpoint2 > other.breakpoint2: return False
        return True

    def __gt__(self, other):
        if self.breakpoint1 < other.breakpoint1: return False
        elif self.breakpoint1 > other.breakpoint1: return True
        if self.breakpoint2 < other.breakpoint2: return False
        elif self.breakpoint2 > other.breakpoint2: return True
        return False

    def __ge__(self, other):
        if self.breakpoint1 < other.breakpoint1: return False
        elif self.breakpoint1 > other.breakpoint1: return True
        if self.breakpoint2 < other.breakpoint2: return False
        elif self.breakpoint2 > other.breakpoint2: return True
        return True

    def __str__(self):
        return "From: ", str(self.breakpoint1), "\n To: ", str(self.breakpoint2), "\n", \
               ("Copy number: %f \n" % self.copy_number), "Genotype: ", (" ".join([str(i) for i in self.genotype]))

    def canonize(self):
        if self.breakpoint2 < self.breakpoint1:
            tmp_breakpoint = self.breakpoint1
            self.breakpoint1 = self.breakpoint2
            self.breakpoint2 = tmp_breakpoint


class Adjacency:
    """
    An edge of the variation graph that does not belong to the reference.

    A *new adjacency* is the observation that two nucleotides are close to one another in a chromosome, despite them not
    being close to one another in the reference: this creates a new edge in the variation graph. This class represents
    only new adjacencies: adjacencies in the reference do not need to be represented explicitly.

    Remark: an adjacency might be involved in no event, i.e. we might have evidence for an adjacency without being able
    to assign it to an event. This is likely to happen with complex events.
    """
    def __init__(self, breakpoint1: Breakpoint, breakpoint2: Breakpoint):
        # One of the two breakpoints might be None, to represent e.g. a
        # chromosomal break that is not rejoined.
        self.breakpoint1: Breakpoint = breakpoint1
        self.breakpoint2: Breakpoint = breakpoint2

        # One count for each individual. The count tells how many haplotypes in the individual (0,1,2,...) use the
        # adjacency at least once. Use -1 for unknown.
        #
        # Remark: we give a genotype to an adjacency, as well, since the same adjacency might be involved in multiple,
        # distinct events in the haplotypes of the same individual.
        self.genotype: list[int] = []

    def __eq__(self, other):
        if isinstance(other, Adjacency):
            return self.breakpoint1 == other.breakpoint1 and self.breakpoint2 == other.breakpoint2
        return False

    def __hash__(self):
        return hash(self.breakpoint1) ^ hash(self.breakpoint2)

    def __lt__(self, other):
        if self.breakpoint1 < other.breakpoint1: return True
        elif self.breakpoint1 > other.breakpoint1: return False
        if self.breakpoint2 < other.breakpoint2: return True
        elif self.breakpoint2 > other.breakpoint2: return False
        return False

    def __le__(self, other):
        if self.breakpoint1 < other.breakpoint1: return True
        elif self.breakpoint1 > other.breakpoint1: return False
        if self.breakpoint2 < other.breakpoint2: return True
        elif self.breakpoint2 > other.breakpoint2: return False
        return True

    def __gt__(self, other):
        if self.breakpoint1 < other.breakpoint1: return False
        elif self.breakpoint1 > other.breakpoint1: return True
        if self.breakpoint2 < other.breakpoint2: return False
        elif self.breakpoint2 > other.breakpoint2: return True
        return False

    def __ge__(self, other):
        if self.breakpoint1 < other.breakpoint1: return False
        elif self.breakpoint1 > other.breakpoint1: return True
        if self.breakpoint2 < other.breakpoint2: return False
        elif self.breakpoint2 > other.breakpoint2: return True
        return True

    def __str__(self):
        return str(self.breakpoint1), " -- ", \
               str(self.breakpoint2), "\n", "Genotype: ", (" ".join([str(i) for i in self.genotype]))

    def canonize(self):
        if self.breakpoint2 < self.breakpoint1:
            tmp_breakpoint = self.breakpoint1
            self.breakpoint1 = self.breakpoint2
            self.breakpoint2 = tmp_breakpoint


class Event:
    """
    A variant of arbitrary complexity. This is either a walk in the variation graph (if the variant is supported by
    evidence of new context), or a set of intervals of the reference with anomalous copy number (if the variant is
    supported just by evidence of anomalous coverage), or both.

    Remark: this is similar to the notion of event in the VCF specification.
    """
    def __init__(self):
        # One of the known SV types defined in ``constants.py``.
        self.event_type: int = SV_TYPE.UNK

        # Sequence of new adjacencies that describes a walk in the variation graph, in arbitrary orientation. The order
        # of this list is important. Consecutive adjacencies are assumed to be connected by intervals of the reference,
        # taken in the orientation determined by the adjacencies. The sequence can be None.
        #
        # Remark: an event might contain just one new adjacency (e.g. in a deletion or in a telomere-telomere fusion).
        #
        # Remark: several adjacencies in the same event might use the same breakpoint.
        self.adjacencies: list[Adjacency] = []

        # List of intervals of the reference with anomalous copy number but with no evidence of new context (might be
        # None).
        self.reference_intervals: list[Reference_Interval] = []

        # One count for every individual. The count specifies how many haplotypes in the individual contain the event
        # (0,1,2,...). Use -1 for unknown.
        #
        # Remark: the reference is always assumed not to contain the event.
        self.genotype: list[int] = []

    def __eq__(self, other):
        if isinstance(other, Event):
            return self.event_type == other.event_type and \
                   self.adjacencies == other.adjacencies and \
                   self.reference_intervals == other.reference_intervals
        return False

    def __hash__(self):
        return hash(self.event_type) ^ hash(self.adjacencies) ^ hash(self.reference_intervals)

    def __str__(self):
        return ("Type: %d \n" % (self.event_type)), \
               "Adjacencies: \n", ('\n'.join([str(i) for i in self.adjacencies])), \
               "Intervals: \n", ('\n'.join([str(i) for i in self.reference_intervals])), \
               "Genotype: ", (" ".join([str(i) for i in self.genotype]))

    def clear(self):
        self.event_type = SV_TYPE.UNK
        self.adjacencies = []
        self.reference_intervals = []
        self.genotype = []

    def canonize(self):
        for a in self.adjacencies: a.canonize()
        reversed_list = copy.copy(self.adjacencies)
        reversed_list.reverse()
        if reversed_list < self.adjacencies: self.adjacencies = reversed_list
        self.reference_intervals.sort()


class Callset:
    """
    Collects all the distinct objects of a given type. The same object might be pointed to by several other objects.
    Identical insertion strings are collapsed into a single object.
    """
    def __init__(self):
        self.events: dict[int, list[Event]] = {}
        self.reference_intervals: dict[int, list[Reference_Interval]] = {}
        self.adjacencies: dict[int, list[Adjacency]] = {}
        self.breakpoints: dict[int, list[Breakpoint]] = {}
        self.sequences: dict[int, list[Sequence]] = {}

    def clear(self):
        self.events = {}
        self.reference_intervals = {}
        self.adjacencies = {}
        self.breakpoints = {}
        self.sequences = {}

    def serialize(self, output_file: str):
        """
        All core data structures. ``pickle`` ensures that the same object is stored only once when it is pointed to by
        several other objects.
        """
        with open(output_file, "wb") as f:
            pickle.dump(self.events, f)
            pickle.dump(self.reference_intervals, f)
            pickle.dump(self.adjacencies, f)
            pickle.dump(self.breakpoints, f)
            pickle.dump(self.sequences, f)

    def deserialize(self, input_file: str):
        """
        All core data structures
        """
        with open(input_file, "rb") as f:
            self.events = pickle.load(f)
            self.reference_intervals = pickle.load(f)
            self.adjacencies = pickle.load(f)
            self.breakpoints = pickle.load(f)
            self.sequences = pickle.load(f)

    def check_consistency(self):
        """
        TODO: makes sure that all the data structures are in a valid state...
        """
        pass

    def sort_and_compact_breakpoints(self):
        """
        TODO: Merges similar breakpoints...
        """
        pass




# ---------------------- SIMPLIFICATION PROCEDURES -----------------------

def relax_breakpoints_with_track(breakpoints, track: Track, min_intersection: int):
    """
    If the uncertainty interval ``A`` of a breakpoint overlaps with an interval ``B`` in ``track`` by
    ``>=min_intersection`` bps, then ``A`` is reset to ``A union B``.

    Remark: the procedure assumes ``breakpoints`` to be sorted by ``position_first``.

    :param track: procedure ``merge_track_intervals()`` is assumed to have already been executed.
    """
    if (breakpoints is None) or (len(breakpoints) == 0) or (track is None) or (len(track.intervals) == 0): return
    i = 0
    j = 0
    n_breakpoints = len(breakpoints)
    n_tracks = len(track.intervals)
    while i < n_breakpoints and j < n_tracks:
        if breakpoints[i].position_last < track.intervals[j].start + min_intersection - 1:
            i += 1
            continue
        if breakpoints[i].position_first > track.intervals[j].end - min_intersection + 1:
            j += 1
            continue
        breakpoints[i].position_first = min(breakpoints[i].position_first, track.intervals[j].start)
        breakpoints[i].position_last = max(breakpoints[i].position_last, track.intervals[j].end)
        breakpoints[i].update_position(True)
        i += 1


def merge_track_intervals(track: Track):
    """
    Merges every pair of overlapping intervals in ``track``.
    """
    j = 0
    n_intervals = len(track.intervals)
    for i in range(1,n_intervals):
        if track.intervals[i].start>track.intervals[j].end:
            j += 1
            track.intervals[j].chr_name = track.intervals[i].chr_name
            track.intervals[j].start = track.intervals[i].start
            track.intervals[j].end = track.intervals[i].end
        else:
            track.intervals[j].end = max(track.intervals[j].end,track.intervals[i].end)
    del track.intervals[j+1:]
