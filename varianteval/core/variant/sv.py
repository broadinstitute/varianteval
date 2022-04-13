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
import varianteval.core.constants as constants



class Event:
    """A human-understandable variant.
    
    This is either a walk in the variation graph (if the variant is supported by
    evidence of new context), or a set of intervals of the reference with 
    anomalous copy number (if the variant is supported just by coverage
    evidence).
    """
    
    # One of the known SV types defined in $constants.py$.
    event_type: int
    
    # Sequence of new adjacencies that descibes a walk in the variation graph
    # (might be None). Consecutive adjacencies are assumed to be connected by
    # intervals of the reference, traversed in the orientation determined by the
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




class Adjacency:
    """An edge of the variation graph that does not belong to the reference. 
    
    A \emph{new adjacency} is the observation that two nucleotides are close to 
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




class Reference_Interval:
    """An interval of the reference, with anomalous copy number.
    
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




class Breakpoint:
    """An uncertain position in a sequence.
    
    Every breakpoint must be involved in an adjacency or in a reference 
    interval.
    """
    
    sequence: Sequence
    
    # If we cut $sequence$ at position $p$, we need to know whether we are
    # using $[1..p]$ or $[p..n]$ in the adjacency. These are defined in
    # $constants.py$.
    side: int
    
    # To model uncertainty, a breakpoint position is a probability distribution
    # over an interval (the interval has size one iff there is no uncertainty).
    #
    # Remark: some callers might not return an interval, but just a probability
    # distribution. In this case we set $position_first = position_last = -1$.
    #
    # Remark: if $sequence$ is circular, then it can happen that 
    # $position_first>position_last$.
    position_first: int  # Inclusive, zero-based.
    position_last: int  # Inclusive, zero-based.
    position_avg: int
    position_std: float
    position_probability_function: int  # ID of a function
    
    # One count for each individual. The count tells how many haplotypes in the
    # individual (0,1,2,...) use an adjacency that involves the breakpoint. Use
    # -1 for unknown.
    #
    # Remark: we give a genotype to a breakpoint, as well, since the same
    # breakpoint might be involved in multiple, distinct events in the 
    # haplotypes of the same individual.
    genotype: list[int]
    
    
    def relax_by_radius(radius: int):
        """Ensures that $[position_avg-radius..position_avg+radius]$ belongs to 
        the uncertainty interval.
        """
        position_first = max(0,min(position_avg-radius,position_min))
        position_last = min(max(position_avg+radius,position_max),sequence.length_max)
        if position_probability_function == constants.PROBABILITY_FUNCTION.DELTA:
            position_probability_function == constants.PROBABILITY_FUNCTION.UNIFORM
            position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
        elif position_probability_function == constants.PROBABILITY_FUNCTION.UNIFORM:
            position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
        else:
            pass  # We keep the STD of the original function
    
    
    @staticmethod
    def build_precise_beakpoint(sequence: Sequence, position: int, side: int) -> Breakpoint:
        """Builds a breakpoint with no uncertainty"""
        
        out = Breakpoint()
        out.sequence = sequence
        out.side = side
        out.position_first = position
        out.position_last = position
        out.position_avg = position
        out.position_std = 0
        out.position_probability_function = constants.PROBABILITY_FUNCTION.DELTA
        return out


    @staticmethod
    def build_uniform_breakpoint(sequence: Sequence, position_first: int, position_last: int, side: int) -> Breakpoint:
        """Builds a breakpoint with uniform uncertainty"""
        
        out = Breakpoint()
        out.sequence = sequence
        out.side = side
        out.position_first = position_first
        out.position_last = position_last
        out.position_avg = (position_last+position_first)//2;
        out.position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
        out.position_probability_function = constants.PROBABILITY_FUNCTION.UNIFORM
        return out




class Sequence:
    """A distinct chromosome, contig, or insertion sequence.
    
    Remark: a sequence might have no breakpoint associated with it.
    """
    
    # General properties
    name: str
    is_circular: bool
    sequence: str  # Nucleotides. Might be NULL if not needed.
    tracks: list[Track]
    
    # The length of the sequence may be uncertain (e.g. if it is an insertion).
    length_min: int
    length_max: int
    length_avg: int
    length_std: float
    length_probability_function: int  # ID of a function
    
    
    @staticmethod
    def build_precise_sequence(name: str, is_circular: bool, length: int, sequence: str) -> Sequence:
        """Builds a sequence with no length uncertainty.
        
        Args:
           sequence: can be None.
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

       
    @staticmethod
    def build_uniform_sequence(name: str, is_circular: bool, length_min: int, length_max: int, sequence: str) -> Sequence:
        """Builds a sequence with uniform length uncertainty
        
        Args:
           sequence: can be None.
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




class Track:
    """An annotation of the intervals of a sequence, like in genome browsers."""
    
    name: str
    
    # List of $(first,last)$ coordinates, zero-based.
    intervals: list[tuple[int,int]]
    
    
    def merge_intervals(min_overlap: int):
        """Merges every pair of intervals that share at least $min_overlap$ 
        bases
        """




def relax_breakpoints_by_track(breakpoints: list[Breakpoint], track: Track):
    """If the uncertainty interval A of a breakpoint overlaps with an interval B
    in $track$, then A is reset to A \union B.
    
    Args:
       breakpoints: assumed to be sorted by $position_first$;
       track: procedure $merge_intervals()$ is assumed to have already been
         executed.
    """
    if (breakpoints is None) or len(breakpoints) == 0 or (track is None) or len(track) == 0: return
    i = 0; j = 0
    n_breakpoints = len(breakpoints)
    n_tracks = len(tracks)
    while i<n_breakpoints and j<n_tracks:
        if breakpoints[i].position_last < track[j][0]:
            i += 1
            continue
        if breakpoints[i].position_first > track[j][1]:
            j += 1
            continue
        breakpoints[i].position_first = min(breakpoints[i].position_first,track[j][0])
        breakpoints[i].position_last = max(breakpoints[i].position_last,track[j][1])
        i += 1
    