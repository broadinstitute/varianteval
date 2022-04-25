"""
Basic data structures for representing variants.

Remark: phasing is outside the scope of this work. We do not model phased 
genotypes on purpose.

Remark: the VCF format allows calls and genotypes to be annotated with several 
fields (see e.g. the FILTER and INFO fields). None of these annotations is 
crucial for our purposes. If needed, we could handle them by storing pointers to
unstructured strings or to external files.
"""

import math
import pickle
import pysam

import varianteval.core.constants as constants
import varianteval.core.variant.intervals as intervals
import varianteval.core.utils as utils


# All distinct objects of a given type. The same object might be pointed to by
# several other objects.
events: set[Event]
reference_intervals: set[Reference_Interval]
adjacencies: set[Adjacency]
breakpoints: set[Breakpoint]
sequences: set[Sequence]


def serialize(output_file: str):
    """
    All core data structures. ``pickle`` ensures that the same object is stored
    only once when it is pointed to by several other objects.
    """
    pickler = pickle.Pickler(output_file,pickle.HIGHEST_PROTOCOL)
    pickler.dump(events)
    pickler.dump(reference_intervals)
    pickler.dump(adjacencies)
    pickler.dump(breakpoints)
    pickler.dump(sequences)


def deserialize(input_file: str):
    """
    All core data structures
    """
    unpickler = pickle.Unpickler(input_file)
    unpickler.load(events)
    unpickler.load(reference_intervals)
    unpickler.load(adjacencies)
    unpickler.load(breakpoints)
    unpickler.load(sequences)


def new_precise_beakpoint(sequence: Sequence, position: int, side: int) -> Breakpoint:
    """ Builds a breakpoint with no uncertainty """
    
    out = Breakpoint()
    out.sequence = sequence
    out.side = side
    out.position_first = position
    out.position_last = position
    out.position_avg = position
    out.position_std = 0
    out.position_probability_function = constants.PROBABILITY_FUNCTION.DELTA
    return out


def new_uniform_breakpoint(sequence: Sequence, position_first: int, position_last: int, side: int) -> Breakpoint:
    """ Builds a breakpoint with uniform uncertainty """
    
    out = Breakpoint()
    out.sequence = sequence
    out.side = side
    out.position_first = position_first
    out.position_last = position_last
    out.position_avg = (position_last+position_first)//2;
    out.position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
    out.position_probability_function = constants.PROBABILITY_FUNCTION.UNIFORM
    return out


def new_precise_sequence(name: str, is_circular: bool, length: int, sequence: str) -> Sequence:
    """
    Builds a sequence with no length uncertainty.
    
    :param sequence: can be None.
    """
    
    out = Sequence()
    out.name = name
    out.is_circular = is_circular
    out.sequence = sequence
    out.length_min = length
    out.length_max = length
    out.length_avg = length
    out.length_std = 0
    out.length_probability_function = constants.PROBABILITY_FUNCTION.DELTA
    return out

   
def new_uniform_sequence(name: str, is_circular: bool, length_min: int, length_max: int, sequence: str) -> Sequence:
    """
    Builds a sequence with uniform length uncertainty.
    
    :param sequence: can be None.
    """
    
    out = Sequence()
    out.name = name
    out.is_circular = is_circular
    out.sequence = sequence
    out.length_first = length_min
    out.length_max = length_max
    out.length_avg = (length_max+length_min)//2;
    out.length_std = math.sqrt(float((length_max-length_min+1)**2)/12)
    out.length_probability_function = constants.PROBABILITY_FUNCTION.UNIFORM
    return out


def relax_breakpoints_with_track(breakpoints: list[Breakpoint], track: Track, min_intersection: int):
    """
    If the uncertainty interval ``A`` of a breakpoint overlaps with an interval 
    ``B`` in ``track`` by ``>=min_intersection`` bps, then ``A`` is reset to 
    ``A \union B``.
    
    :param breakpoints: assumed to be sorted by ``position_first``;
    :param track: procedure ``merge_intervals()`` is assumed to have already
         been executed.
    """
    if (breakpoints is None) or len(breakpoints) == 0 or (track is None) or len(track) == 0: return
    i = 0; j = 0
    n_breakpoints = len(breakpoints)
    n_tracks = len(tracks)
    while i<n_breakpoints and j<n_tracks:
        if breakpoints[i].position_last < track[j][1]+min_intersection-1:
            i += 1
            continue
        if breakpoints[i].position_first > track[j][2]-min_intersection+1:
            j += 1
            continue
        breakpoints[i].position_first = min(breakpoints[i].position_first,track[j][1])
        breakpoints[i].position_last = max(breakpoints[i].position_last,track[j][2])
        breakpoints[i].update_position()
        i += 1
    



# ----------------------------- DATA STRUCTURES --------------------------------

class Event:
    """
    A human-understandable variant.
    
    This is either a walk in the variation graph (if the variant is supported by
    evidence of new context), or a set of intervals of the reference with 
    anomalous copy number (if the variant is supported just by evidence of
    anomalous coverage).
    """
    
    # One of the known SV types defined in ``constants.py``.
    event_type: int
    
    # Sequence of new adjacencies that descibes a walk in the variation graph
    # (might be None). Consecutive adjacencies are assumed to be connected by
    # intervals of the reference, taken in the orientation determined by the
    # adjacencies.
    #
    # Remark: an event might contain just one new adjacency (e.g. in a deletion
    # or in a telomere-telomere fusion).
    #
    # Remark: several adjacencies in the same event might use the same
    # breakpoint.
    adjacencies: list[Adjacency]
    
    # List of intervals of the reference with anomalous copy number (might be
    # None).
    reference_intervals: list[Reference_Interval]
    
    # One count for every individual. The count specifies how many haplotypes
    # in the individual contain the event (0,1,2,...). Use -1 for unknown.
    # 
    # Remark: the reference is assumed not to contain the event.
    genotype: list[int]
    
    
    def __eq__(self, other):
        if isinstance(other,Event):
            return self.event_type == other.event_type and \
                   self.adjacencies == other.adjacencies and \
                   self.reference_intervals == other.reference_intervals
        return False
    
    
    def __hash__(self):
        return hash(self.event_type) ^ hash(self.adjacencies) ^ hash(self.reference_intervals)
    
    
    def __str__(self):
        return ("Type: %d \n" % (self.event_type)), "Adjacencies: \n", ('\n'.join([str(i) for i in self.adjacencies])), "Intervals: \n", ('\n'.join([str(i) for i in self.reference_intervals])), "Genotype: ", (" ".join([str(i) for i in self.genotype]))



class Adjacency:
    """
    An edge of the variation graph that does not belong to the reference. 
    
    A *new adjacency* is the observation that two nucleotides are close to 
    one another in a chromosome, despite them not being close to one another in 
    the reference: this creates a new edge in the variation graph. This class 
    represents only new adjacencies: adjacencies in the reference do not need to
    be represented explicitly.
    
    Remark: an adjacency might be involved in no event, i.e. we might have 
    evidence for an adjacency without being able to assign it to an event. This
    is likely to happen with complex events.
    """
    
    # One of the two breakpoints might be None, to represent e.g. a chromosomal
    # break that is not rejoined.
    position1, position2: Breakpoint
    
    # One count for each individual. The count tells how many haplotypes in the
    # individual (0,1,2,...) use the adjacency at least once. Use -1 for
    # unknown.
    #
    # Remark: we give a genotype to an adjacency, as well, since the same
    # adjacency might be involved in multiple, distinct events in the 
    # haplotypes of the same individual.
    genotype: list[int]
    
    
    def __init__(self, position1: Breakpoint, position2: Breakpoint):
        self.position1 = position1
        self.position2 = position2
    
    
    def __eq__(self, other):
        if isinstance(other,Adjacency):
            return self.position1 == other.position1 and self.position2 == other.position2
        return False
    
    
    def __hash__(self):
        return hash(self.position1) ^ hash(self.position2)
    
    
    def __str__(self):
        return str(self.position1), " -- ", str(self.position2), "\n", "Genotype: ", (" ".join([str(i) for i in self.genotype]))



class Reference_Interval:
    """
    An interval of the reference, with anomalous copy number.
    
    This class models the case in which we have evidence for a clear change in 
    coverage, but we don't have evidence for new contexts.
    """
    
    # Cannot be None.
    first_position, last_position: Breakpoint
    copy_number: float
    
    # One count for each individual. The count tells how many haplotypes in the
    # individual (0,1,2,...) contain this interval with this specific copy
    # number. Use -1 for unknown.
    #
    # Remark: we give a genotype to an interval, as well, just for generality.
    # I.e. we allow the same interval to be involved in multiple, distinct 
    # events in the haplotypes of the same individual.
    genotype: list[int]
    
    
    def __init__(self, first_position: Breakpoint, last_position: Breakpoint, copy_number: float):
        self.first_position = first_position
        self.last_position = last_position
        self.copy_number = copy_number
    
    
    def __eq__(self, other):
        if isinstance(other,Reference_Interval):
            return self.first_position == other.first_position and \
                   self.last_position == other.last_position
        return False
    
    
    def __hash__(self):
        return hash(self.first_position) ^ hash(self.last_position)
    
    
    def __str__(self):
        return "From: ", str(self.first_position), "\n To: ", str(self.last_position), "\n", ("Copy number: %f \n" % self.copy_number), "Genotype: ", (" ".join([str(i) for i in self.genotype]))



class Breakpoint:
    """
    An uncertain position in a sequence.
    
    Every breakpoint must be involved in an adjacency or in a reference 
    interval.
    """
    
    sequence: Sequence
    
    # If we cut a sequence of length ``n`` at position ``p`` (zero-based), we
    # need to know whether we are using ``[0..p]`` (true) or ``[p..n-1]``
    # (false) in the adjacency. These are defined in ``constants.py``.
    side: int
    
    # To model uncertainty, a breakpoint position is a probability distribution
    # over an interval (the interval has size one iff there is no uncertainty).
    #
    # Remark: some callers might not return an interval, but just a probability
    # distribution. In this case we set ``position_first = position_last = -1``.
    #
    # Remark: if ``sequence`` is circular, then it can happen that 
    # ``position_first > position_last``.
    position_first: int  # Inclusive, zero-based.
    position_last: int  # Inclusive, zero-based.
    position_avg: int
    position_std: float
    position_probability_function: int  # Defined in ``constants.py``.
    
    # One count for each individual. The count tells how many haplotypes in the
    # individual (0,1,2,...) use an adjacency that involves the breakpoint. Use
    # -1 for unknown.
    #
    # Remark: we give a genotype to a breakpoint, as well, since the same
    # breakpoint might be involved in multiple, distinct events in the 
    # haplotypes of the same individual.
    genotype: list[int]

    
    def __eq__(self, other):
        if isinstance(other,Breakpoint):
            return self.sequence == other.sequence and \
                   self.side == other.side and \
                   self.position_first == other.position_first and \
                   self.position_last == other.position_last and \
                   self.position_avg == other.position_avg and \
                   self.position_std == other.position_std and \
                   self.position_probability_function == other.position_probability_function
        return False
    
    
    def __hash__(self):
        return hash(self.sequence) ^ hash(self.side) ^ hash(self.position_first) ^ hash(self.position_last) ^ hash(self.position_avg) ^ hash(self.position_std) ^ hash(self.position_probability_function)
    
    
    def __str__(self):
        if self.side == constants.BREAKPOINT_SIDE.LEFT:
            left = "<"
            right = "]"
        else:
            left: "["
            right: ">"
        return ("%s %s%s%d..(%d,%f)..%d%s" % (self.sequence.name, self.position_probability_function, left, self.position_first, self.position_avg, self.position_std, self.position_last, right)), "Genotype: ", (" ".join([str(i) for i in self.genotype]))
    
    
    def update_position():
        """ 
        Given ``position_{first,last}``, resets all other ``position_*`` fields.
        """
        if position_first == position_last:
            position_probability_function == constants.PROBABILITY_FUNCTION.DELTA
            position_avg == position_first
            position_std = 0
        else:
            if position_probability_function == constants.PROBABILITY_FUNCTION.DELTA || position_probability_function == constants.PROBABILITY_FUNCTION.UNIFORM:
                position_probability_function == constants.PROBABILITY_FUNCTION.UNIFORM
                position_avg = (position_first+position_last)//2
                position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
            else:
                pass  # We keep the STD and AVG of the original function
    
    
    def relax_by_radius(radius: int):
        """
        Ensures that ``[position_avg-radius..position_avg+radius]`` belongs to 
        the uncertainty interval.
        """
        position_first = max(0,min(position_avg-radius,position_first))
        position_last = min(max(position_avg+radius,position_last),sequence.length_max)
        if position_probability_function == constants.PROBABILITY_FUNCTION.DELTA:
            position_probability_function == constants.PROBABILITY_FUNCTION.UNIFORM
            position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
        elif position_probability_function == constants.PROBABILITY_FUNCTION.UNIFORM:
            position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
        else:
            pass  # We keep the STD of the original function



class Sequence:
    """
    A distinct chromosome, contig, or insertion sequence.
    
    Remark: a sequence might have no breakpoint associated with it.
    """
    
    # General properties
    name: str  # Can be None if it is not a known sequence.
    is_circular: bool
    sequence: str  # The actual nucleotides. Can be None if not needed.
    tracks: list[Track]
    
    # The length of the sequence may be uncertain (e.g. if it is an insertion).
    length_min: int
    length_max: int
    length_avg: int
    length_std: float
    length_probability_function: int  # ID of a function
    
    
    def __eq__(self, other):
        if isinstance(other,Sequence):
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
        return hash(self.name) ^ hash(self.is_circular) ^ hash(sequence) ^ hash(self.length_min) ^ hash(self.length_max) ^ hash(self.length_avg) ^ hash(self.length_std) ^ hash(self.length_probability_function)
        
    
    def __str__(self):
        return ("%s(%s) %s[%d..(%d,%f)..%d]\n" % (self.name, self.is_circular, self.length_probability_function, self.length_min, self.length_avg, self.length_std, self.length_max)), "Tracks: \n", ("\n".join([str(i) for i in self.tracks]))
    
    
    def get_type() -> int:
        """
        :return: 0=chromosome, 1=contig, 2=inserted string.
        """
        if util.is_chromosome(name):
            return 0
        else if name == constants.VCF_INSERTION_STRING_NAME:
            return 2
        return 1



class Track:
    """ An annotation of a sequence with intervals """
    
    name: str
    
    # List of ``[first..last]`` coordinates, inclusive, zero-based.
    intervals: list[intervals.GenomeInterval]
    
    
    def __eq__(self, other):
        if isinstance(other,Track):
            return self.name == other.name and self.intervals == other.intervals
        return False
    
    
    def __hash__(self):
        return hash(self.name) ^ hash(tuple(self.intervals))
    
    
    def __str__(self):
        return ("%s \n" % self.name), "\n".join([str(i) for i in self.intervals]) 
    
    
    def merge_intervals():
        """ Merges every pair of overlapping intervals """
        j = 0
        n_intervals = len(intervals)
        for i in range(1:n_intervals):
            if intervals[i].start>intervals[j].end:
                j += 1
                intervals[j].chr_name = intervals[i].chr_name
                intervals[j].start = intervals[i].start
                intervals[j].end = intervals[i].end
            else:
                intervals[j].end = max(intervals[j].end,intervals[i].end)
        del intervals[j+1:]
