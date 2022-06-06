from pysam import VariantFile, VariantRecord
from typing import Any, Callable
import copy

from varianteval.core.constants import *
import varianteval.core.vcf as corevcf
import varianteval.core.utils as utils
import varianteval.core.variant.sv as sv


# --- VCF field accessors and iteration utils for the pysam VCF library ---
def get_sv_type(pysam_vcf_record):
    return SV_TYPE_ENCODING[pysam_vcf_record.info[VCF_SVTYPE_STR]] if VCF_SVTYPE_STR in pysam_vcf_record.info else SV_TYPE.UNK


def get_gt(pysam_vcf_record, sample_id=0):
    sample_name = pysam_vcf_record.samples[sample_id].name
    if VCF_GT_STR in pysam_vcf_record.samples[sample_name]:
        gt = pysam_vcf_record.samples[sample_name][VCF_GT_STR]
    else:
        gt = (None, None)
    return VCF2GT_ENCODING[gt]


def get_sv_len(pysam_vcf_record):
    if VCF_SVLEN_STR in pysam_vcf_record.info:
        return abs(int(pysam_vcf_record.info[VCF_SVLEN_STR] if not isinstance(pysam_vcf_record.info[VCF_SVLEN_STR], tuple)
                       else pysam_vcf_record.info[VCF_SVLEN_STR][0]))
    else:
        return abs(pysam_vcf_record.stop - pysam_vcf_record.pos)


def record2sv(record):
    """Converts the VCF record into the interval SV event object"""
    sv_type = get_sv_type(record)
    sv_len = get_sv_len(record)
    # TODO: construct the internal representation
    return SV(sv_type)


def filter_record(record, filter_funcs=None):
    """Applies a set of filters to the record"""
    if filter_funcs is None:
        filter_funcs = []
    pass


def vcf_iter(vcf_fname, filter_funcs=None):
    vcf_file = VariantFile(vcf_fname)
    for rec in vcf_file.fetch():
        if not filter_record(rec, filter_funcs):
            # TODO: handle multi-line events
            yield record2sv(rec)
    # consolidate multi-line calls


def vcf2callset(vcf_fname, filter_funcs=None):
    """Loads a VCF into an SV callset"""
    vcf_file = VariantFile(vcf_fname)
    callset = SVCallset()
    for rec in vcf_file.fetch():
        if not filter_record(rec, filter_funcs):
            # TODO: handle multi-line events
            callset.add(record2sv(rec))
    # TODO: consolidate multi-line calls
    return callset


def svs2vcf(vcf_fname, sv_callset):
    """Outputs a VCF given an SV callset"""
    pass


def preprocess_vcf(vcf_fname):
    pass


def get_breakend_chromosome(string: str) -> str:
    """
    :param str: the ALT field of a BND record of a VCF file.
    """
    q: int = string.find(":")
    p: int = string.find("[")
    if p == -1: p = string.find("]")
    return string[p+1:q]


def get_breakend_chromosome_position(string: str) -> int:
    """
    :param str: the ALT field of a BND record of a VCF file.
    """
    q: int = string.find(":")
    p: int = string.rfind("[")
    if p < 0: p = string.rfind("]")
    return int(string[q+1:p]) - 1  # VCF positions are one-based


def get_breakend_chromosome_side(string: str, position_id: bool) -> int:
    """
    :param string: the ALT field of a BND record of a VCF file;
    :param position_id: TRUE=first position in the VCF record (encoded by the 
    CHROM and POS fields); FALSE=second position in the VCF record (encoded by 
    the ALT field).
    :return: value defined in ``constants.BREAKPOINT_SIDE``.
    """
    if position_id:
        if string[0] == "[" or string[0] == "]":
            return BREAKPOINT_SIDE.LEFT
        else:
            return BREAKPOINT_SIDE.RIGHT
    else:
        if string.find("[") >= 0:
            return BREAKPOINT_SIDE.LEFT
        else:
            return BREAKPOINT_SIDE.RIGHT


def get_confidence_interval(record: VariantRecord, is_first_pos: bool, sigma_multiple: int) -> tuple[int, int]:
    """
    Returns the range of uncertainty around a position ``X`` (implicit from the 
    context).
    
    :param is_first_pos: TRUE=first position in a VCF record; FALSE=last 
    position;
    :param sigma_multiple: if the record does not contain coordinates but just 
    the STD, the procedure assumes an interval that extends this multiple of the
    STD on both sides (2 or 3 captures most of a Gaussian);
    :return: First element: quantity to be added to ``X`` to get the first 
    position of the uncertainty interval (typically negative or zero). Second 
    element: quantity to be added to ``X`` to get the last position of the 
    uncertainty interval (typically positive or zero). If ``record`` does not 
    contain uncertainty information on one side of the interval, the
    corresponding element of the tuple is set to zero. If ``record``does not
    contain any uncertainty information on ``X``, both elements of the tuple are
    set to zero and ``X`` is assumed to be certain.
    """
    quantum: int
    value: float	

    if is_first_pos:
        if VCF_CIPOS_STR in record.info.keys():
            return tuple(map(int, record.info.get(VCF_CIPOS_STR).split(VCF_CI_SEPARATOR)))
            # The sign is already ok
    else:
        if VCF_CIEND_STR in record.info.keys():
            return tuple(map(int, record.info.get(VCF_CIEND_STR).split(VCF_CI_SEPARATOR)))
            # The sign is already ok
        if VCF_CILEN_STR in record.info.keys():
            return tuple(map(int, record.info.get(VCF_CILEN_STR).split(VCF_CI_SEPARATOR)))
            # The sign is already ok
    for key in (VCF_STD_START if is_first_pos else VCF_STD_END):
        if key in record.info.keys():
            value = float(record.info.get(key))
            quantum = int(value*SIGMA_MULTIPLE)
            return (-quantum, quantum)
    return (0, 0)


def get_confidence_interval_std(record: VariantRecord, is_first_pos: bool) -> float:
    """
    Returns the STD of the range of uncertainty around a position ``X`` 
    (implicit from the context).
    
    :param is_first_pos: TRUE=first position in a VCF record; FALSE=last 
    position.
    """

    for key in (VCF_STD_START if is_first_pos else VCF_STD_END):
        if key in record.info.keys():
            return float(record.info.get(key))                
    return 0.0


def ct2side(record: VariantRecord) -> tuple[int,int]:
    """
    Reads info about the CT field of a TRA record.
    
    :return: First element: left side. Second element: right side.
    """

    value = record.info.get(VCF_CT_STR)
    if value == VCF_CT_325_STR:
        return (1, 0)
    elif value == VCF_CT_523_STR:
        return (0, 1)
    elif value == VCF_CT_525_STR:
        return (0, 0)
    elif value == VCF_CT_323_STR:
        return (1, 1)
    else:
        raise Exception("Invalid value of the CT field")




# --------------------- LOADING THE MAIN DATA STRUCTURES -----------------------

def load_vcf(path: str, record_filter: Callable[[VariantRecord],bool], contig_lengths: dict[str,int], callset: sv.Callset):
    """
    Adds to ``callset`` all the objects related to a given VCF file.
    
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
    tmp_sequence = sv.Sequence()
    tmp_breakpoint = sv.Breakpoint()
    tmp_adjacency = sv.Adjacency()
    tmp_interval = sv.Reference_Interval()
    tmp_event = sv.Event()
    vcf_file = VariantFile(path,'r')
    for record in vcf_file.fetch():
        if not record_filter(record): continue
        sv_type = corevcf.get_sv_type(record)
        if sv_type in [SV_TYPE.DEL, SV_TYPE.DEL_ME]: 
            load_event_del(record,contig_lengths,callset,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type in [SV_TYPE.INS, SV_TYPE.INS_ME, SV_TYPE.INS_NOVEL]: 
            load_event_ins(record,contig_lengths,callset,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type in [SV_TYPE.DUP, SV_TYPE.DUP_TANDEM]:
            load_event_dup(record,contig_lengths,callset,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.DUP_INT:
            # Warning: don't know how to load sparse dup records yet.
            pass
        elif sv_type == SV_TYPE.CNV: 
            load_event_cnv(record,contig_lengths,callset,tmp_event,tmp_interval,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.INV: 
            load_event_inv(record,contig_lengths,callset,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.INV_DUP: 
            load_event_invdup(record,contig_lengths,callset,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.DEL_INV:
            # Warning: a DEL/INV is assumed to be a sequence ``aC'e`` where
            # ``C'`` is the reverse-complement of ``C`` and the reference is
            # ``aBCDe``. For now we assume that the length of B,D is unknown,
            # so we cannot create adjacencies.
            pass
        elif sv_type == SV_TYPE.BND: 
            load_event_bnd(record,False,contig_lengths,callset,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.TRA: 
            load_event_bnd(record,True,contig_lengths,callset,tmp_event,tmp_adjacency,tmp_breakpoint,tmp_sequence)
        elif sv_type == SV_TYPE.CPX:
            # Warning: don't know how to handle generic complex calls yet.
            pass
        elif sv_type == SV_TYPE.UNK:
            # Warning: don't know how to handle unknown calls yet.
            pass
    vcf_file.close()


def load_event_del(record: VariantRecord, contig_lengths: dict[str,int], callset: sv.Callset, tmp_event: sv.Event, tmp_adjacency: sv.Adjacency, tmp_breakpoint: sv.Breakpoint, tmp_sequence: sv.Sequence):
    """
    Converts ``record`` into a deletion event and stores it in dictionary
    ``callset.events``.
    
    Remark: we assume that a DEL record means that the new adjacency between the
    end and the beginning of the deleted unit is observed by the caller. This
    should be implemented better by using field SVCLAIM from the VCF 4.4 spec.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.DEL
    
    load_sequence(record,True,contig_lengths,callset,tmp_sequence)
    sequence = retrieve_instance(tmp_sequence,callset.sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    tmp_adjacency.breakpoint1 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    tmp_adjacency.breakpoint2 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,callset.adjacencies))
    retrieve_instance(tmp_event,callset.events)


def load_event_ins(record: VariantRecord, contig_lengths: dict[str,int], callset: sv.Callset, tmp_event: sv.Event, tmp_adjacency: sv.Adjacency, tmp_breakpoint: sv.Breakpoint, tmp_sequence: sv.Sequence):
    """
    Converts ``record`` into an insertion event and stores it in dictionary 
    ``callset.events``.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.INS
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence1 = retrieve_instance(tmp_sequence,callset.sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence1
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    breakpoint1 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence1
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint4 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    
    load_sequence(record,False,contig_lengths,tmp_sequence)    
    sequence2 = retrieve_instance(tmp_sequence,callset.sequences)
    set_precise_beakpoint(tmp_breakpoint,sequence2,0,BREAKPOINT_SIDE.LEFT)
    breakpoint2 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    set_precise_beakpoint(tmp_breakpoint,sequence2,sequence2.length-1,BREAKPOINT_SIDE.RIGHT)
    breakpoint3 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    
    tmp_adjacency.breakpoint1 = breakpoint1
    tmp_adjacency.breakpoint2 = breakpoint2
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,callset.adjacencies))
    tmp_adjacency.breakpoint1 = breakpoint3
    tmp_adjacency.breakpoint2 = breakpoint4
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,callset.adjacencies))
    retrieve_instance(tmp_event,callset.events)


def load_event_dup(record: VariantRecord, contig_lengths: dict[str,int], callset: sv.Callset, tmp_event: sv.Event, tmp_adjacency: sv.Adjacency, tmp_breakpoint: sv.Breakpoint, tmp_sequence: sv.Sequence):
    """
    Converts ``record`` into a tandem duplication event and stores it in
    dictionary ``callset.events``. 
    
    Remark: we assume that a DUP record means that the new adjacency between the
    end and the beginning of the tandem unit is observed by the caller. This
    should be implemented better by using field SVCLAIM from the VCF 4.4 spec.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.INS
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence = retrieve_instance(tmp_sequence,callset.sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    tmp_adjacency.breakpoint1 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    tmp_adjacency.breakpoint2 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,callset.adjacencies))
    retrieve_instance(tmp_event,callset.events)


def load_event_cnv(record: VariantRecord, contig_lengths: dict[str,int], callset: sv.Callset, tmp_event: sv.Event, tmp_interval: sv.Reference_Interval, tmp_breakpoint: sv.Breakpoint, tmp_sequence: sv.Sequence):
    """
    Converts ``record`` into a copy-number variation event and stores it in
    dictionary ``callset.events``. We assume that a CNV record means that the 
    evidence for a change in copy number comes from coverage only, i.e. that no
    new adjacency is observed by the caller. This follows the VCF 4.4 
    spec. of SVCLAIM.
    
    Warning: the procedure does not set the ``copy_number`` field of the
    reference interval yet.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.CNV
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence = retrieve_instance(tmp_sequence,callset.sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    tmp_interval.breakpoint1 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    tmp_breakpoint.copy_number = 0.0
    tmp_interval.breakpoint2 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    tmp_event.reference_intervals.append(retrieve_instance(tmp_interval,callset.intervals))
    retrieve_instance(tmp_event,callset.events)


def load_event_inv(record: VariantRecord, contig_lengths: dict[str,int], callset: sv.Callset, tmp_event: sv.Event, tmp_adjacency: sv.Adjacency, tmp_breakpoint: sv.Breakpoint, tmp_sequence: sv.Sequence):
    """
    Converts ``record`` into an inversion event and stores it in dictionary 
    ``callset.events``.
    """
    tmp_event.clear()
    tmp_event.event_type = SV_TYPE.INV
    
    load_sequence(record,True,contig_lengths,tmp_sequence)
    sequence = retrieve_instance(tmp_sequence,callset.sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint2 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    tmp_breakpoint.move_left(1)
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    breakpoint1 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    breakpoint3 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    tmp_breakpoint.move_right(1)
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint4 = retrieve_instance(tmp_breakpoint,callset.breakpoints)

    tmp_adjacency.breakpoint1 = retrieve_instance(breakpoint1,callset.breakpoints)
    tmp_adjacency.breakpoint2 = retrieve_instance(breakpoint3,callset.breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,callset.adjacencies))
    tmp_adjacency.breakpoint1 = retrieve_instance(breakpoint2,callset.breakpoints)
    tmp_adjacency.breakpoint2 = retrieve_instance(breakpoint4,callset.breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,callset.adjacencies))
    retrieve_instance(tmp_event,callset.events)


def load_event_invdup(record: VariantRecord, contig_lengths: dict[str,int], callset sv.Callset, tmp_event: sv.Event, tmp_adjacency: sv.Adjacency, tmp_breakpoint: sv.Breakpoint, tmp_sequence: sv.Sequence):
    """
    Converts ``record`` into an inverted duplication event and stores it in
    dictionary ``callset.events``. An INVDUP is assumed to be a sequence 
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
    sequence = retrieve_instance(tmp_sequence,callset.sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint2 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence
    tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
    breakpoint3 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    tmp_breakpoint.move_right(1)
    tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    breakpoint4 = retrieve_instance(tmp_breakpoint,callset.breakpoints)

    tmp_adjacency.breakpoint1 = retrieve_instance(breakpoint3,callset.breakpoints)
    tmp_adjacency.breakpoint2 = retrieve_instance(breakpoint3,callset.breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,callset.adjacencies))
    tmp_adjacency.breakpoint1 = retrieve_instance(breakpoint2,callset.breakpoints)
    tmp_adjacency.breakpoint2 = retrieve_instance(breakpoint4,callset.breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,callset.adjacencies))
    retrieve_instance(tmp_event,callset.events)


def load_event_bnd(record: VariantRecord, is_tra: bool, contig_lengths: dict[str,int], callset: sv.Callset, tmp_event: sv.Event, tmp_adjacency: sv.Adjacency, tmp_breakpoint: sv.Breakpoint, tmp_sequence: sv.Sequence):
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
    sequence1 = retrieve_instance(tmp_sequence,callset.sequences)
    load_breakpoint(record,True,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence1
    if is_tra:
        tmp_tuple = corevcf.ct2side(record)
        if tmp_tuple[0] == 1:
            tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
        else:
            tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    else:
        alt_field = record.alts[0]
        tmp_breakpoint.side = corevcf.get_breakend_chromosome_side(alt_field,True)
    tmp_adjacency.breakpoint1 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    load_sequence(record,False,contig_lengths,tmp_sequence)    
    sequence2 = retrieve_instance(tmp_sequence,callset.sequences)
    load_breakpoint(record,False,tmp_breakpoint)
    tmp_breakpoint.sequence = sequence2
    if is_tra:
        if tmp_tuple[1] == 1:
            tmp_breakpoint.side = BREAKPOINT_SIDE.RIGHT
        else:
            tmp_breakpoint.side = BREAKPOINT_SIDE.LEFT
    else:
        alt_field = record.alts[0]
        tmp_breakpoint.side = corevcf.get_breakend_chromosome_side(alt_field,False)
    tmp_adjacency.breakpoint2 = retrieve_instance(tmp_breakpoint,callset.breakpoints)
    tmp_event.adjacencies.append(retrieve_instance(tmp_adjacency,callset.adjacencies))


def retrieve_instance(query: Any, dictionary: Any) -> Any:
    """
    If ``dictionary`` contains an object that is equivalent to ``query``, the 
    procedure returns that object. Otherwise, it creates a shallow copy of 
    ``query`` and adds it to ``dictionary``.
    
    Remark: this procedure is designed to work on different types.
    
    :return: the object in ``dictionary`` that is equivalent to ``query``.
    """
    
    if isinstance(query,sv.Event) or isinstance(query,sv.Reference_Interval) or isinstance(query,sv.Adjacency): query.canonize()
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


def load_breakpoint(record: VariantRecord, is_first_pos: bool, breakpoint: sv.Breakpoint):
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
        interval = corevcf.get_confidence_interval(record,True,SIGMA_MULTIPLE)
        breakpoint.position_first = breakpoint.position_avg + interval[0]
        breakpoint.position_last = breakpoint.position_avg + interval[1]
        breakpoint.position_std = corevcf.get_confidence_interval_std(record,True)
        if breakpoint.position_std != 0.0:
            breakpoint.position_probability_function = PROBABILITY_FUNCTION.GAUSSIAN
        elif breakpoint.position_first == breakpoint.position_last:
            breakpoint.position_probability_function = PROBABILITY_FUNCTION.DELTA
        else:
            breakpoint.position_probability_function = PROBABILITY_FUNCTION.UNIFORM
    else:
        sv_type = corevcf.get_sv_type(record)
        # Remark: VCF's pos is one-based, but it actually refers to the position
        # before the variant.
        position = int(record.pos)
        if sv_type in [SV_TYPE.DEL, SV_TYPE.DEL_ME]:
            length = int(record.info.get(VCF_SVLEN_STR))
            if length < 0: length = -length
            breakpoint.position_avg = position + length - 1
            interval = corevcf.get_confidence_interval(record,False,SIGMA_MULTIPLE)
            breakpoint.position_first = breakpoint.position_avg + interval[0]
            breakpoint.position_last = breakpoint.position_avg + interval[1]
            breakpoint.position_std = corevcf.get_confidence_interval_std(record,False)
            if breakpoint.position_std != 0.0:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.GAUSSIAN
            elif breakpoint.position_first == breakpoint.position_last:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.DELTA
            else:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.UNIFORM
        elif sv_type in [SV_TYPE.INS, SV_TYPE.INS_ME, SV_TYPE.INS_NOVEL]:
            breakpoint.position_avg = int(record.pos) + 1
            interval = corevcf.get_confidence_interval(record,True,SIGMA_MULTIPLE)
            breakpoint.position_first = breakpoint.position_avg + interval[0]
            breakpoint.position_last = breakpoint.position_avg + interval[1]
            breakpoint.position_std = corevcf.get_confidence_interval_std(record,True)
            if breakpoint.position_std != 0.0:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.GAUSSIAN
            elif breakpoint.position_first == breakpoint.position_last:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.DELTA
            else:
                breakpoint.position_probability_function = PROBABILITY_FUNCTION.UNIFORM
        elif sv_type in [SV_TYPE.DUP, SV_TYPE.DUP_TANDEM, SV_TYPE.INV, SV_TYPE.CNV, SV_TYPE.INV_DUP, SV_TYPE.DEL_INV]:
            # Remark: VCF's end position is one-based.
            breakpoint.position_avg = int(record.info.get(VCF_END_STR)) - 1
            interval = corevcf.get_confidence_interval(record,False,SIGMA_MULTIPLE)
            breakpoint.position_first = breakpoint.position_avg + interval[0]
            breakpoint.position_last = breakpoint.position_avg + interval[1]
            breakpoint.position_std = corevcf.get_confidence_interval_std(record,False)
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
            position = corevcf.get_breakend_chromosome_position(record.alts[0])
            set_precise_beakpoint(breakpoint,None,position,-1)
        elif sv_type == SV_TYPE.TRA:
            # Setting the other breakpoint to precise. This is arbitrary.
            # Remark: VCF's end position is one-based.
            position = int(record.info.get(VCF_END_STR)) - 1
            set_precise_beakpoint(breakpoint,None,position,-1)


def load_sequence(record: VariantRecord, is_first_sequence: bool, contig_lengths: dict[str,int], sequence: sv.Sequence):
    """
    Loads in ``sequence`` the information about the sequence of the first or 
    second position in ``record`` (the first position is the one encoded by
    columns CHROM,POS).
    
    :param contig_lengths: keys are assumed to contain every contig name used in
    ``record`` (lowercase).
    """
    IDENTITY_THRESHOLD = 1  # We tolerate off-by-one errors in SVLEN
    
    sequence.clear()
    sv_type = corevcf.get_sv_type(record)
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
        contig_name = corevcf.get_breakend_chromosome(alt_field)
        is_circular = utils.is_mitochondrion(contig_name)
        set_precise_sequence(sequence,contig_name,False,contig_lengths[contig_name],None)


def new_precise_beakpoint(sequence: sv.Sequence, position: int, side: int) -> sv.Breakpoint:
    """ Builds a breakpoint with no uncertainty """
    
    out = sv.Breakpoint()
    set_precise_beakpoint(out,sequence,position,side)
    return out


def set_precise_beakpoint(breakpoint: sv.Breakpoint, sequence: sv.Sequence, position: int, side: int):
    breakpoint.sequence = sequence
    breakpoint.side = side
    breakpoint.position_first = position
    breakpoint.position_last = position
    breakpoint.position_avg = position
    breakpoint.position_std = 0
    breakpoint.position_probability_function = PROBABILITY_FUNCTION.DELTA


def new_uniform_breakpoint(sequence: sv.Sequence, position_first: int, position_last: int, side: int) -> sv.Breakpoint:
    """ Builds a breakpoint with uniform uncertainty """
    
    out = sv.Breakpoint()
    set_uniform_breakpoint(out,sequence,position_first,position_last,side)
    return out


def set_uniform_breakpoint(breakpoint: sv.Breakpoint, sequence: sv.Sequence, position_first: int, position_last: int, side: int):
    breakpoint.sequence = sequence
    breakpoint.side = side
    breakpoint.position_first = position_first
    breakpoint.position_last = position_last
    breakpoint.position_avg = (position_last+position_first)//2;
    breakpoint.position_std = math.sqrt(float((position_last-position_first+1)**2)/12)
    breakpoint.position_probability_function = PROBABILITY_FUNCTION.UNIFORM


def new_precise_sequence(name: str, is_circular: bool, length: int, characters: str) -> sv.Sequence:
    """
    Builds a sequence with no length uncertainty.
    
    :param sequence: can be None.
    """
    
    out = sv.Sequence()
    set_precise_sequence(out,name,is_circular,length,characters)
    return out


def set_precise_sequence(out: sv.Sequence, name: str, is_circular: bool, length: int, characters: str):
    out.name = name
    out.is_circular = is_circular
    out.sequence = characters
    out.length_min = length
    out.length_max = length
    out.length_avg = length
    out.length_std = 0
    out.length_probability_function = PROBABILITY_FUNCTION.DELTA

   
def new_uniform_sequence(name: str, is_circular: bool, length_min: int, length_max: int, characters: str) -> sv.Sequence:
    """
    Builds a sequence with uniform length uncertainty.
    
    :param sequence: can be None.
    """
    
    out = sv.Sequence()
    set_uniform_sequence(out,name,is_circular,length_min,length_max,characters)
    return out


def set_uniform_sequence(out: sv.Sequence, name: str, is_circular: bool, length_min: int, length_max: int, characters: str):
    out.name = name
    out.is_circular = is_circular
    out.sequence = characters
    out.length_first = length_min
    out.length_max = length_max
    out.length_avg = (length_max+length_min)//2;
    out.length_std = math.sqrt(float((length_max-length_min+1)**2)/12)
    out.length_probability_function = PROBABILITY_FUNCTION.UNIFORM
