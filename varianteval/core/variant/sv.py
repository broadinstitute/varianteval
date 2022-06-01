"""
The main data structures for representing variants.

.. image:: ../../../docs/sv.png

Remark: phasing is not currently modeled and it is left to future work.
"""

import copy
import math
import pickle
from typing import Any, Callable
import pysam

import varianteval.core.constants
import varianteval.core.vcf as vcf
import varianteval.core.utils as utils
import varianteval.core.variant.intervals as intervals


# All distinct objects of a given type. The same object might be pointed to by
# several other objects. Identical insertion strings are collapsed into a single
# object.
events: dict[int, list[Event]] = {}
reference_intervals: dict[int, list[Reference_Interval]] = {}
adjacencies: dict[int, list[Adjacency]] = {}
breakpoints: dict[int, list[Breakpoint]] = {}
sequences: dict[int, list[Sequence]] = {}




# ---------------------- DATA STRUCTURES LOADING PROCEDURES --------------------

def load_vcf(path: str, record_filter: Callable[[VariantRecord],bool], contig_lengths: dict[str,int]):
    """
    Adds to the global dictionaries all the objects related to a given VCF file.
    
    Remark: only contigs involved in some event are loaded in ``sequences``. The
    basepairs of a contig are not loaded.
    
    Warning: genotype information is not loaded, and objects are not tagged with
    a sample ID yet. This if left to the future.
    
    Warning: the procedure assumes just one ALT per record. This is done just
    for simplicity and should be generalized in the future.
    
    :param record_filter: a function that decides which VCF records to keep
    (e.g. only records with PASS annotation, or only long SVs);
    :param contig_lengths: the keys are assumed to contain every contig name
    used in ``path`` (lowercase).
    """
    tmp_sequence = Sequence()
    tmp_breakpoint = Breakpoint()
    tmp_adjacency = Adjacency()
    tmp_interval = Reference_Interval()
    tmp_event = Event()
    vcf_file = VariantFile(path,'r')
    for record in vcf_file.fetch():
        if not record_filter(record): continue
        sv_type = vcf.get_sv_type(record)
        if sv_type in [SV_TYPE.DEL, SV_TYPE.DEL_ME]: 
            load_event_del(record,contig_lengths,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type in [SV_TYPE.INS, SV_TYPE.INS_ME, SV_TYPE.INS_NOVEL]: 
            load_event_ins(record,contig_lengths,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type in [SV_TYPE.DUP, SV_TYPE.DUP_TANDEM]:
            load_event_dup(record,contig_lengths,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.DUP_INT:
            # Warning: don't know how to load sparse dup records yet.
            pass
        elif sv_type == SV_TYPE.CNV: 
            load_event_cnv(record,contig_lengths,tmp_event,tmp_interval,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.INV: 
            load_event_inv(record,contig_lengths,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.INV_DUP: 
            load_event_invdup(record,contig_lengths,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.DEL_INV:
            # Warning: a DEL/INV is assumed to be a sequence ``aC'e`` where
            # ``C'`` is the reverse-complement of ``C`` and the reference is
            # ``aBCDe``. For now we assume that the length of B,D is unknown,
            # so we cannot create adjacencies.
            pass
        elif sv_type == SV_TYPE.BND: 
            load_event_bnd(record,False,contig_lengths,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.TRA: 
            load_event_bnd(record,True,contig_lengths,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.CPX:
            # Warning: don't know how to handle generic complex calls yet.
            pass
        elif sv_type == SV_TYPE.UNK:
            # Warning: don't know how to handle unknown calls yet.
            pass
    vcf_file.close()


def load_event_del(record: VariantRecord, contig_lengths: dict[str,int], tmp_event: Event, tmp_adjacency: Adjacency, tmp_breakpoint: Breakpoint, tmp_sequence: Sequence):
    """
    Converts ``record`` into a deletion event and stores it in global dictionary
    ``events``.
    
    Remark: we assume that a DEL record means that the new adjacency between the
    end and the beginning of the deleted unit is observed by the caller. This
    should be implemented better by using field SVCLAIM from the VCF 4.4 spec.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.DEL
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence = retrieve_instance(tmp_sequence,sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    tmp_adjacency.breakpoint1 = retrieve_instance(tmp_breakpoint,breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    tmp_adjacency.breakpoint2 = retrieve(tmp_breakpoint,breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,adjacencies))
    retrieve_instance(tmp_event,events)


def load_event_ins(record: VariantRecord, contig_lengths: dict[str,int], tmp_event: Event, tmp_adjacency: Adjacency, tmp_breakpoint: Breakpoint, tmp_sequence: Sequence):
    """
    Converts ``record`` into an insertion event and stores it in global 
    dictionary ``events``.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.INS
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence1 = retrieve_instance(tmp_sequence,sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence1
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    breakpoint1 = retrieve_instance(tmp_breakpoint,breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence1
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint4 = retrieve_instance(tmp_breakpoint,breakpoints)
    
    load_sequence(record,False,contig_lengths,tmp_sequence)    
    sequence2 = retrieve_instance(tmp_sequence,sequences)
    set_precise_beakpoint(tmp_breakpoint,sequence2,0,BREAKPOINT_SIDE.LEFT)
    breakpoint2 = retrieve_instance(tmp_breakpoint,breakpoints)
    set_precise_beakpoint(tmp_breakpoint,sequence2,sequence2.length-1,BREAKPOINT_SIDE.RIGHT)
    breakpoint3 = retrieve_instance(tmp_breakpoint,breakpoints)
    
    tmp_adjacency.breakpoint1 = breakpoint1
    tmp_adjacency.breakpoint2 = breakpoint2
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,adjacencies))
    tmp_adjacency.breakpoint1 = breakpoint3
    tmp_adjacency.breakpoint2 = breakpoint4
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,adjacencies))
    retrieve_instance(tmp_event,events)


def load_event_dup(record: VariantRecord, contig_lengths: dict[str,int], tmp_event: Event, tmp_adjacency: Adjacency, tmp_breakpoint: Breakpoint, tmp_sequence: Sequence):
    """
    Converts ``record`` into a tandem duplication event and stores it in global
    dictionary ``events``. 
    
    Remark: we assume that a DUP record means that the new adjacency between the
    end and the beginning of the tandem unit is observed by the caller. This
    should be implemented better by using field SVCLAIM from the VCF 4.4 spec.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.INS
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence = retrieve_instance(tmp_sequence,sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    tmp_adjacency.breakpoint1 = retrieve_instance(tmp_breakpoint,breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    tmp_adjacency.breakpoint2 = retrieve_instance(tmp_breakpoint,breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,adjacencies))
    retrieve_instance(tmp_event,events)


def load_event_cnv(record: VariantRecord, contig_lengths: dict[str,int], tmp_event: Event, tmp_interval: Reference_Interval, tmp_breakpoint: Breakpoint, tmp_sequence: Sequence):
    """
    Converts ``record`` into a copy-number variation event and stores it in
    global dictionary ``events``. We assume that a CNV record means that the 
    evidence for a change in copy number comes from coverage only, i.e. that no
    new adjacency is observed by the caller. This follows the VCF 4.4 
    spec. of SVCLAIM.
    
    Warning: the procedure does not set the ``copy_number`` field of the
    reference interval yet.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.CNV
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence = retrieve_instance(tmp_sequence,sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    tmp_interval.breakpoint1 = retrieve_instance(tmp_breakpoint,breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    tmp_breakpoint.copy_number = 0.0
    tmp_interval.breakpoint2 = retrieve_instance(tmp_breakpoint,breakpoints)
    tmp_event.reference_intervals.append(retrieve_instance(tmp_interval,intervals))
    retrieve_instance(tmp_event,events)


def load_event_inv(record: VariantRecord, contig_lengths: dict[str,int], tmp_event: Event, tmp_adjacency: Adjacency, tmp_breakpoint: Breakpoint, tmp_sequence: Sequence):
    """
    Converts ``record`` into an inversion event and stores it in global
    dictionary ``events``.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.INV
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence = retrieve_instance(tmp_sequence,sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint2 = retrieve_instance(tmp_breakpoint,breakpoints)
    tmp_breakpoint.move_left(1)
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    breakpoint1 = retrieve_instance(tmp_breakpoint,breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    breakpoint3 = retrieve_instance(tmp_breakpoint,breakpoints)
    tmp_breakpoint.move_right(1)
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint4 = retrieve_instance(tmp_breakpoint,breakpoints)

    tmp_adjacency.breakpoint1 = retrieve_instance(breakpoint1,breakpoints)
    tmp_adjacency.breakpoint2 = retrieve_instance(breakpoint3,breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,adjacencies))
    tmp_adjacency.breakpoint1 = retrieve_instance(breakpoint2,breakpoints)
    tmp_adjacency.breakpoint2 = retrieve_instance(breakpoint4,breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,adjacencies))
    retrieve_instance(tmp_event,events)


def load_event_invdup(record: VariantRecord, contig_lengths: dict[str,int], tmp_event: Event, tmp_adjacency: Adjacency, tmp_breakpoint: Breakpoint, tmp_sequence: Sequence):
    """
    Converts ``record`` into an inverted duplication event and stores it in
    global dictionary ``events``. An INVDUP is assumed to be a sequence 
    ``xVV'y``, where ``V'`` is the reverse-complement of ``V`` and the reference
    is ``xVy``. This is not equivalent to Sniffles 1's definition (Supplementary 
    Figure 2.4), in which it seems to be ``xVy'V'z`` where the reference is 
    ``xVyz``. So we are interpreting Sniffles 1's INVDUP as if it were Sniffles
    1's DUP+INV (see again Supplementary Figure 2.4).
    
    We assume that an INVDUP record means that the new adjacencies between the 
    end and the end of the tandem unit, and between the beginning and the right
    neighborhood of the tandem unit, are observed by the caller.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.INVDUP
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence = retrieve_instance(tmp_sequence,sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint2 = retrieve_instance(tmp_breakpoint,breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    breakpoint3 = retrieve_instance(tmp_breakpoint,breakpoints)
    tmp_breakpoint.move_right(1)
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint4 = retrieve_instance(tmp_breakpoint,breakpoints)

    tmp_adjacency.breakpoint1 = retrieve_instance(breakpoint3,breakpoints)
    tmp_adjacency.breakpoint2 = retrieve_instance(breakpoint3,breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,adjacencies))
    tmp_adjacency.breakpoint1 = retrieve_instance(breakpoint2,breakpoints)
    tmp_adjacency.breakpoint2 = retrieve_instance(breakpoint4,breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,adjacencies))
    retrieve_instance(tmp_event,events)


def load_event_bnd(record: VariantRecord, is_tra: bool, contig_lengths: dict[str,int], tmp_event: Event, tmp_adjacency: Adjacency, tmp_breakpoint: Breakpoint, tmp_sequence: Sequence):
    """
    Converts ``record`` into an event that contains a single new adjacency.
    
    Warning: instead, we should collect in the same event all the BND records
    with the same VCF EVENT tag... This is left to the future...
    
    :param is_tra: the breakend record is encoded using the TRA (true) or the
    BND (false) convention.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.INVDUP
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence1 = retrieve_instance(tmp_sequence,sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence1
    if is_tra:
        tmp_tuple = vcf.ct2side(record)
        if tmp_tuple[0] == 1:
            tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
        else:
            tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    else:
        alt_field = record.alts[0]
        tmp_breakpoint.side = vcf.get_breakend_chromosome_side(alt_field,True)
    tmp_adjacency.breakpoint1 = retrieve_instance(tmp_breakpoint,breakpoints)
    load_sequence(record,False,contig_lengths,tmp_sequence)    
    sequence2 = retrieve_instance(tmp_sequence,sequences)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence2
    if is_tra:
        if tmp_tuple[1] == 1:
            tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
        else:
            tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    else:
        alt_field = record.alts[0]
        tmp_breakpoint.side = vcf.get_breakend_chromosome_side(alt_field,False)
    tmp_adjacency.breakpoint2 = retrieve_instance(tmp_breakpoint,breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,adjacencies))


def retrieve_instance(query: Any, dictionary: Any) -> Any:
    """
    If ``dictionary`` contains an object that is equivalent to ``query``, the 
    procedure returns that object. Otherwise, it creates a shallow copy of 
    ``query`` and adds it to ``dictionary``.
    
    Remark: this procedure is designed to work on different types.
    
    :return: the object in ``dictionary`` that is equivalent to ``query``.
    """
    
    if isinstance(query,Event) or isinstance(query,Reference_Interval) or isinstance(query,Adjacency): query.canonize()
    out = None
    key = hash(query)
    if key in dictionary.keys():
        i = dictionary[key].index(query)
        if i >= 0:
            out = dictionary[key][i]
        else:
            out = copy.copy(query)
            dictionary[key].append(out)
    else:
        out = copy.copy(query)
        events[key] = [out]
    return out


def load_breakpoint(record: VariantRecord, is_first_pos: bool, breakpoint: Breakpoint):
    """
    Loads in ``breakpoint`` the information about the first or the second
    position in ``record``.
    
    Remark: the procedure does not set the following fields of the breakpoint,
    which are left to the caller to complete: ``sequence, side, genotype``.
    """
    SIGMA_MULTIPLE = 3  # Arbitrary
    
    breakpoint.clear()
    if is_first_pos:
        # Remark: VCF's pos is one-based, but it actually refers to the position
        # before the variant.
        breakpoint.position_avg = int(record.pos)
        interval = vcf.get_confidence_interval(record,True,SIGMA_MULTIPLE)
        breakpoint.position_first = breakpoint.position_avg + interval[0]
        breakpoint.position_last = breakpoint.position_avg + interval[1]
        breakpoint.position_std = vcf.get_confidence_interval_std(record,True)
        if breakpoint.position_std != 0.0:
            breakpoint.position_probability_function = PROBABILITY_FUNCTION.GAUSSIAN
        elif breakpoint.position_first == breakpoint.position_last:
            breakpoint.position_probability_function = PROBABILITY_FUNCTION.DELTA
        else:
            breakpoint.position_probability_function = PROBABILITY_FUNCTION.UNIFORM
    else:
        sv_type = vcf.get_sv_type(record)
        # Remark: VCF's pos is one-based, but it actually refers to the position
        # before the variant.
        position = int(record.pos)
        if sv_type in [SV_TYPE.DEL, SV_TYPE.DEL_ME]:
            length = int(record.info.get(VCF_SVLEN_STR))
            if length < 0: length = -length
            breakpoint.position_avg = position + length - 1
            interval = vcf.get_confidence_interval(record,False,SIGMA_MULTIPLE)
            breakpoint.position_first = breakpoint.position_avg + interval[0]
            breakpoint.position_last = breakpoint.position_avg + interval[1]
            breakpoint.position_std = vcf.get_confidence_interval_std(record,False)
            if breakpoint.position_std != 0.0:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.GAUSSIAN
            elif breakpoint.position_first == breakpoint.position_last:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.DELTA
            else:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.UNIFORM
        elif sv_type in [SV_TYPE.INS, SV_TYPE.INS_ME, SV_TYPE.INS_NOVEL]:
            breakpoint.position_avg = int(record.pos) + 1
            interval = vcf.get_confidence_interval(record,True,SIGMA_MULTIPLE)
            breakpoint.position_first = breakpoint.position_avg + interval[0]
            breakpoint.position_last = breakpoint.position_avg + interval[1]
            breakpoint.position_std = vcf.get_confidence_interval_std(record,True)
            if breakpoint.position_std != 0.0:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.GAUSSIAN
            elif breakpoint.position_first == breakpoint.position_last:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.DELTA
            else:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.UNIFORM
        elif sv_type in [SV_TYPE.DUP, SV_TYPE.DUP_TANDEM, SV_TYPE.INV, SV_TYPE.CNV, SV_TYPE.INV_DUP, SV_TYPE.DEL_INV]:
            # Remark: VCF's end position is one-based.
            breakpoint.position_avg = int(record.info.get(VCF_END_STR)) - 1
            interval = vcf.get_confidence_interval(record,False,SIGMA_MULTIPLE)
            breakpoint.position_first = breakpoint.position_avg + interval[0]
            breakpoint.position_last = breakpoint.position_avg + interval[1]
            breakpoint.position_std = vcf.get_confidence_interval_std(record,False)
            if breakpoint.position_std != 0.0:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.GAUSSIAN
            elif breakpoint.position_first == breakpoint.position_last:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.DELTA
            else:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.UNIFORM
        elif sv_type == SV_TYPE.BND:
            # Setting the other breakpoint to precise for now. Later it should
            # be merged with the description of the same breakpoint in the
            # symmetrical VCF record.
            position = vcf.get_breakend_chromosome_position(record.alts[0])
            set_precise_beakpoint(breakpoint,None,position,-1)
        elif sv_type == SV_TYPE.TRA:
            # Setting the other breakpoint to precise. This is arbitrary.
            # Remark: VCF's end position is one-based.
            position = int(record.info.get(VCF_END_STR)) - 1
            set_precise_beakpoint(breakpoint,None,position,-1)


def load_sequence(record: VariantRecord, is_first_sequence: bool, contig_lengths: dict[str,int], sequence: Sequence):
    """
    Loads in global dictionary ``sequence`` the information about the sequence 
    of the first or second position in ``record`` (the first position is the one
    encoded by columns CHROM,POS).
    
    :param contig_lengths: keys are assumed to contain every contig name used in
    ``record`` (lowercase).
    """
    IDENTITY_THRESHOLD = 1  # We tolerate off-by-one errors in SVLEN
    
    sequence.clear()
    sv_type = vcf.get_sv_type(record)
    if is_first_sequence:
        contig_name = record.contig.lower()
        is_circular = utils.is_mitochondrion(contig_name)
        set_precise_sequence(sequence,contig_name,is_circular,contig_lengths[contig_name],None)
    elif sv_type in [SV_TYPE.INS, SV_TYPE.INS_ME, SV_TYPE.INS_NOVEL]:
        length_prime = record.info.SVLEN
        alt_field = record.alts[0]
        if len(alt_field) == 0 or VCF_INSERTION_ALT in alt_field:
            basepairs = "n" * length_prime
        else:
            if alt_field[0] == record.ref:
                basepairs = alt_field[1:]
            else:     
                basepairs = alt_field
            length = len(basepairs)
            if length_prime < length-IDENTITY_THRESHOLD or length_prime > length+IDENTITY_THRESHOLD:
                print("ERROR: the length of the inserted string is different from the SVLEN field: %d != %d" % (length,length_prime))
                exit()
        set_precise_sequence(sequence,VCF_INSERTION_STRING_NAME,False,length,basepairs)
    elif sv_type == SV_TYPE.BND:
        alt_field = record.alts[0]
        contig_name = vcf.get_breakend_chromosome(alt_field)
        is_circular = utils.is_mitochondrion(contig_name)
        set_precise_sequence(sequence,contig_name,False,contig_lengths[contig_name],None)


def serialize(output_file: str):
    """
    All core data structures. ``pickle`` ensures that the same object is stored
    only once when it is pointed to by several other objects.
    """
    with open(output_file, "wb") as f:
        pickle.dump(events, f)
        pickle.dump(reference_intervals, f)
        pickle.dump(adjacencies, f)
        pickle.dump(breakpoints, f)
        pickle.dump(sequences, f)


def deserialize(input_file: str):
    """
    All core data structures
    """
    with open(input_file, "rb") as f:
        events = pickle.load(f)
        reference_intervals = pickle.load(f)
        adjacencies = pickle.load(f)
        breakpoints = pickle.load(f)
        sequences = pickle.load(f)


def new_precise_beakpoint(sequence: Sequence, position: int, side: int) -> Breakpoint:
    """ Builds a breakpoint with no uncertainty """
    
    out = Breakpoint()
    set_precise_beakpoint(out,sequence,position,side)
    return out


def set_precise_beakpoint(breakpoint: Breakpoint, sequence: Sequence, position: int, side: int):
    breakpoint.sequence = sequence
    breakpoint.side = side
    breakpoint.position_first = position
    breakpoint.position_last = position
    breakpoint.position_avg = position
    breakpoint.position_std = 0
    breakpoint.position_probability_function = PROBABILITY_FUNCTION.DELTA


def new_uniform_breakpoint(sequence: Sequence, position_first: int, position_last: int, side: int) -> Breakpoint:
    """ Builds a breakpoint with uniform uncertainty """
    
    out = Breakpoint()
    set_uniform_breakpoint(out,sequence,position_first,position_last,side)
    return out


def set_uniform_breakpoint(breakpoint: Breakpoint, sequence: Sequence, position_first: int, position_last: int, side: int):
    breakpoint.sequence = sequence
    breakpoint.side = side
    breakpoint.position_first = position_first
    breakpoint.position_last = position_last
    breakpoint.position_avg = (position_last+position_first)//2;
    breakpoint.position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
    breakpoint.position_probability_function = PROBABILITY_FUNCTION.UNIFORM


def new_precise_sequence(name: str, is_circular: bool, length: int, sequence: str) -> Sequence:
    """
    Builds a sequence with no length uncertainty.
    
    :param sequence: can be None.
    """
    
    out = Sequence()
    set_precise_sequence(out,name,is_circular,length,sequence)
    return out


def set_precise_sequence(out: Sequence, name: str, is_circular: bool, length: int, sequence: str):
    out.name = name
    out.is_circular = is_circular
    out.sequence = sequence
    out.length_min = length
    out.length_max = length
    out.length_avg = length
    out.length_std = 0
    out.length_probability_function = PROBABILITY_FUNCTION.DELTA

   
def new_uniform_sequence(name: str, is_circular: bool, length_min: int, length_max: int, sequence: str) -> Sequence:
    """
    Builds a sequence with uniform length uncertainty.
    
    :param sequence: can be None.
    """
    
    out = Sequence()
    set_uniform_sequence(out,name,is_circular,length_min,length_max,sequence)
    return out


def set_uniform_sequence(Sequence: out, name: str, is_circular: bool, length_min: int, length_max: int, sequence: str):
    out.name = name
    out.is_circular = is_circular
    out.sequence = sequence
    out.length_first = length_min
    out.length_max = length_max
    out.length_avg = (length_max+length_min)//2;
    out.length_std = math.sqrt(float((length_max-length_min+1)**2)/12)
    out.length_probability_function = PROBABILITY_FUNCTION.UNIFORM


def check_consistency():
    """
    TODO: makes sure that all the data structures are in a valid state...
    """
    pass




# -------------------- BREAKPOINT SIMPLIFICATION PROCEDURES --------------------


def sort_and_compact_breakpoints():
    """
    TODO: Merges similar breakpoints...
    """
    pass


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
    i = 0
    j = 0
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
        breakpoints[i].update_position(True)
        i += 1




# ---------------------- TRACK SIMPLIFICATION PROCEDURES -----------------------


def merge_track_intervals(track: Track):
    """ Merges every pair of overlapping intervals in ``track``. """
    j = 0
    n_intervals = len(track.intervals)
    for i in range(1:n_intervals):
        if track.intervals[i].start>track.intervals[j].end:
            j += 1
            track.intervals[j].chr_name = track.intervals[i].chr_name
            track.intervals[j].start = track.intervals[i].start
            track.intervals[j].end = track.intervals[i].end
        else:
            track.intervals[j].end = max(track.intervals[j].end,track.intervals[i].end)
    del track.intervals[j+1:]




# -------------------------- MAIN DATA STRUCTURES ------------------------------

class Event:
    """
    A variant of arbitrary complexity. This is either a walk in the variation 
    graph (if the variant is supported by evidence of new context), or a set of 
    intervals of the reference with anomalous copy number (if the variant is 
    supported just by evidence of anomalous coverage), or both.
    
    Remark: this is similar to the notion of event in the VCF specification.
    """
    
    def __init__(self):
       # One of the known SV types defined in ``constants.py``.
       self.event_type: int = SV_TYPE.UNK
    
       # Sequence of new adjacencies that descibes a walk in the variation
       # graph, in arbitrary orientation. The order of this list is important.
       # Consecutive adjacencies are assumed to be connected by intervals of
       # the reference, taken in the orientation determined by the adjacencies.
       # The sequence can be None.
       #
       # Remark: an event might contain just one new adjacency (e.g. in a
       # deletion or in a telomere-telomere fusion).
       #
       # Remark: several adjacencies in the same event might use the same
       # breakpoint.
       self.adjacencies: list[Adjacency] = []
    
       # List of intervals of the reference with anomalous copy number but with
       # no evidence of new context (might be None).
       self.reference_intervals: list[Reference_Interval] = []
    
       # One count for every individual. The count specifies how many haplotypes
       # in the individual contain the event (0,1,2,...). Use -1 for unknown.
       # 
       # Remark: the reference is always assumed not to contain the event.
       self.genotype: list[int] = []
    
    
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
        
    
    def clear():
        event_type = SV_TYPE.UNK
        adjacencies = []
        reference_intervals = []
        genotype = []
        
    
    def canonize():
        for a in adjacencies: a.canonize()
        reversed_list = copy.copy(adjacencies)
        reversed_list.reverse()
        if reversed_list < adjacencies: adjacencies = reversed_list
        reference_intervals.sort()
        

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
    
    def __init__(self, breakpoint1: Breakpoint, breakpoint2: Breakpoint):
        # One of the two breakpoints might be None, to represent e.g. a
        # chromosomal break that is not rejoined.
        self.breakpoint1: breakpoint1
        self.breakpoint2: breakpoint2
    
        # One count for each individual. The count tells how many haplotypes in
        # the individual (0,1,2,...) use the adjacency at least once. Use -1 for
        # unknown.
        #
        # Remark: we give a genotype to an adjacency, as well, since the same
        # adjacency might be involved in multiple, distinct events in the 
        # haplotypes of the same individual.
        self.genotype: list[int] = []
    
    
    def __eq__(self, other):
        if isinstance(other,Adjacency):
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
        return str(self.breakpoint1), " -- ", str(self.breakpoint2), "\n", "Genotype: ", (" ".join([str(i) for i in self.genotype]))
    
    
    def canonize():
        if breakpoint2 < breakpoint1:
            tmp_breakpoint = breakpoint1
            breakpoint1 = breakpoint2
            breakpoint2 = tmp_breakpoint


class Reference_Interval:
    """
    An interval of the reference, with anomalous copy number.
    
    This class models the case in which the caller has evidence for a clear 
    change in coverage, but it doesn't have evidence for new context.
    """
    
    def __init__(self, breakpoint1: Breakpoint, breakpoint2: Breakpoint, copy_number: float):
        # Cannot be None.
        self.breakpoint1: Breakpoint = breakpoint1
        self.breakpoint2: Breakpoint = breakpoint2
        self.copy_number: float = copy_number
    
        # One count for each individual. The count tells how many haplotypes in
        # the individual (0,1,2,...) contain this interval with this specific
        # copy number. Use -1 for unknown.
        #
        # Remark: we give a genotype to an interval, as well, just for
        # generality. I.e. we allow the same interval to be involved in
        # multiple, distinct events in the haplotypes of the same individual.
        self.genotype: list[int] = []
    
    
    def __eq__(self, other):
        if isinstance(other,Reference_Interval):
            return self.breakpoint1 == other.breakpoint1 and \
                   self.breakpoint2 == other.breakpoint2 and \
                   self.copy_number == other.copy_number
        return False
    
    
    def __hash__(self):
        return hash(self.breakpoint1) ^ hash(self.breakpoint2) ^ hash(copy_number)
    
    
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
        return "From: ", str(self.breakpoint1), "\n To: ", str(self.breakpoint2), "\n", ("Copy number: %f \n" % self.copy_number), "Genotype: ", (" ".join([str(i) for i in self.genotype]))
    
    
    def canonize():
        if breakpoint2 < breakpoint1:
            tmp_breakpoint = breakpoint1
            breakpoint1 = breakpoint2
            breakpoint2 = tmp_breakpoint


class Breakpoint:
    """
    Assume that every position of the reference is a strip of paper with a left
    side and a right side. A breakpoint is an uncertain position with certain 
    side.
    """
    
    def __init__(self): 
        self.sequence: Sequence = None
    
        # If we cut a sequence of length ``n`` at position ``p`` (zero-based),
        # we need to know whether we are using ``[0..p]`` (true) or ``[p..n-1]``
        # (false) in an adjacency. These are defined in ``constants.py``.
        self.side: int = BREAKPOINT_SIDE.UNKNOWN
    
        # To model uncertainty, a breakpoint position is a probability
        # distribution over an interval (the interval has size one iff there is
        # no uncertainty).
        #
        # Remark: some callers might not return an interval, but just a
        # probability distribution. In this case we set ``position_first =
        # position_last = -1``.
        #
        # Remark: if ``sequence`` is circular, then it can happen that 
        # ``position_first > position_last``.
        self.position_first: int = -1  # Inclusive, zero-based.
        self.position_last: int = -1  # Inclusive, zero-based.
        self.position_avg: int = -1
        self.position_std: float = 0.0
        self.position_probability_function: int = PROBABILITY_FUNCTION.UNKNOWN  # Defined in ``constants.py``.
    
        # One count for each individual. The count tells how many haplotypes in
        # the individual (0,1,2,...) use an adjacency that involves the
        # breakpoint. Use -1 for unknown.
        #
        # Remark: we give a genotype to a breakpoint, as well, since the same
        # breakpoint might be involved in multiple, distinct events in the 
        # haplotypes of the same individual.
        self.genotype: list[int] = []
        
    
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
            left = "<"
            right = "]"
        else:
            left: "["
            right: ">"
        return ("%s %s%s%d..(%d,%f)..%d%s" % (self.sequence.name, self.position_probability_function, left, self.position_first, self.position_avg, self.position_std, self.position_last, right)), "Genotype: ", (" ".join([str(i) for i in self.genotype]))
    
    
    def clear():
        sequence = None
        side = BREAKPOINT_SIDE.UNKNOWN
        position_first = -1
        position_last = -1
        position_avg = -1
        position_std = -1
        position_probability_function = PROBABILITY_FUNCTION.UNKNOWN
    
    
    def is_similar(other: Breakpoint, identity_threshold: int, jaccard_threshold: float) -> bool:
        """
        :return: TRUE iff the two breakpoints are intervals of size one and are
        at most ``identity_threshold`` bps apart, or if at least one breakpoint
        is an interval of size greater than one and the two intervals have
        Jaccard similarity at least ``jaccard_threshold``.
        """
        if self.sequence != other.sequence or self.side != other.side:
            return False
        if position_first == position_last:
            if other.position_first == other.position_last:
                return abs(position_avg, other.position_avg) <= identity_threshold
            else:
                jaccard = float(min(other.position_last,position_last)-max(other.position_first,position_first)+1) / (max(other.position_last,position_last)-min(other.position_first,position_first)+1)
                return jaccard >= jaccard_threshold
        else:
            jaccard = float(min(other.position_last,position_last)-max(other.position_first,position_first)+1) / (max(other.position_last,position_last)-min(other.position_first,position_first)+1)
            return jaccard >= jaccard_threshold
        
    
    def update_position(update_avg: bool):
        """ 
        Given ``position_{first,last}``, resets all other ``position_*`` fields.
        
        :param update_avg: if FALSE, ``position_avg`` is not updated.
        """
        if position_first == position_last:
            position_probability_function == PROBABILITY_FUNCTION.DELTA
            if update_avg:
                position_avg == position_first
            position_std = 0
        else:
            if position_probability_function in [PROBABILITY_FUNCTION.DELTA, PROBABILITY_FUNCTION.UNIFORM]:
                position_probability_function == PROBABILITY_FUNCTION.UNIFORM
                if update_avg:
                    position_avg = (position_first+position_last)//2
                position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
            else:
                pass  # We keep the STD and AVG of the original function
    

    def move_right(bps: int):
        position_avg = min(position_avg+bps, sequence.length_max-1)
        position_first = min(position_first+bps, sequence.length_max-1)
        position_last = min(position_last+bps, sequence.length_max-1)
        update_position(False)
        
        
    def move_left(bps: int):
        position_avg = max(position_avg-bps, 0)
        position_first = max(position_first-bps, 0)
        position_last = max(position_last-bps, 0)
        update_position(False)

    
    def relax_by_radius(radius: int):
        """
        Ensures that ``[position_avg-radius..position_avg+radius]`` belongs to 
        the uncertainty interval.
        """
        position_first = max(0,min(position_avg-radius,position_first))
        position_last = min(max(position_avg+radius,position_last),sequence.length_max)
        if position_probability_function == PROBABILITY_FUNCTION.DELTA:
            position_probability_function == PROBABILITY_FUNCTION.UNIFORM
            position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
        elif position_probability_function == PROBABILITY_FUNCTION.UNIFORM:
            position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
        else:
            pass  # We keep the STD of the original function


class Sequence:
    """
    A distinct chromosome, contig, or inserted sequence.
    Only sequences with some breakpoint should be loaded in memory.
    """
    
    def __init__(self):
        # General properties
        self.name: str = None  # Can be None if it is not a known sequence.
        self.is_circular: bool = False
        self.sequence: str = None  # The actual nucleotides. Can be None if not needed.
        self.tracks: list[Track] = []
    
        # The length of the sequence may be uncertain (e.g. if it is an
        # insertion).
        self.length_min: int = 0
        self.length_max: int = 0
        self.length_avg: int = 0
        self.length_std: float = 0.0
        self.length_probability_function: int = PROBABILITY_FUNCTION.UNKNOWN  # ID of a function
    
    
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
        return ("%s(%s) %s[%d..(%d,%f)..%d]\n" % (self.name, self.is_circular, self.length_probability_function, self.length_min, self.length_avg, self.length_std, self.length_max)), "Tracks: \n", ("\n".join([str(i) for i in self.tracks]))
    
    
    def clear():
        name = None
        is_circular = false
        sequence = None
        tracks = None
        length_min = 0
        length_max = 0
        length_avg = 0
        length_std = 0
        length_probability_function = PROBABILITY_FUNCTION.UNKNOWN
    
    
    def get_type() -> int:
        """
        :return: 0=chromosome, 1=contig, 2=inserted string.
        """
        if util.is_chromosome(name): return 0
        elif name == VCF_INSERTION_STRING_NAME: return 2
        return 1


class Track:
    """ Interval annotation of a sequence """
    
    def __init__(self):
        self.name: str = None
    
        # List of ``[first..last]`` coordinates, inclusive, zero-based.
        self.intervals: list[intervals.GenomeInterval] = []
    
    
    def __eq__(self, other):
        if isinstance(other,Track):
            return self.name == other.name and self.intervals == other.intervals
        return False
    
    
    def __hash__(self):
        return hash(self.name) ^ hash(tuple(self.intervals))
        
    
    def __str__(self):
        return ("%s \n" % self.name), "\n".join([str(i) for i in self.intervals])
