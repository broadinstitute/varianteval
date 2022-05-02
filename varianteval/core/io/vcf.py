from pysam import VariantFile
from pysam import VariantRecord

import varianteval.core.constants


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
    q = string.find(":")
    p = string.find("[")
    if p == -1: p = string.find("]")
    return string[p+1:q]


def get_breakend_chromosome_position(string: str) -> int:
    """
    :param str: the ALT field of a BND record of a VCF file.
    """
	q = string.find(":")
	p = string.rfind("[")
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
    else
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
    uncertainty interval (typically positive or zero).
    """
	quantum: int
    value: float	
    
	if is_first_pos:
		if VCF_CIPOS_STR in record.info.keys():
            return tuple(map(int, record.info.get(VCF_CIPOS_STR).split(VCF_CI_SEPARATOR)))
            # The sign is already ok
        for string in VCF_STD_START:
    		if string in record.info.keys():
                value = float(record.info.get(string))
                quantum = int(value*SIGMA_MULTIPLE)
                return (-quantum, quantum)
	}
	else:
		if VCF_CIEND_STR in record.info.keys():
            return tuple(map(int, record.info.get(VCF_CIEND_STR).split(VCF_CI_SEPARATOR)))
            # The sign is already ok
		if VCF_CILEN_STR in record.info.keys():
            return tuple(map(int, record.info.get(VCF_CILEN_STR).split(VCF_CI_SEPARATOR)))
            # The sign is already ok
        for string in VCF_STD_END:
    		if string in record.info.keys():
                value = float(record.info.get(string))
                quantum = int(value*SIGMA_MULTIPLE)
                return (-quantum, quantum)
	return (0, 0)
}


def get_confidence_interval_std(record: VariantRecord, is_first_pos: bool) -> float:
    """
    Returns the STD of the range of uncertainty around a position ``X`` 
    (implicit from the context).
    
    :param is_first_pos: TRUE=first position in a VCF record; FALSE=last 
    position.
    """
	if is_first_pos:
        for string in VCF_STD_START:
    		if string in record.info.keys():
                return float(record.info.get(string))
	else:
        for string in VCF_STD_END:
    		if string in record.info.keys():
                return float(record.info.get(string))
	return 0.0
}


def ct2side(record: VariantRecord) -> tuple[int,int]:
    """
    Reads info about the CT field of a TRA record.
    
    :return: First element: left side. Second element: right side.
    """
    out: tuple[int,int]
    
    value = record.info.get(VCF_CT_STR)
	if value == VCF_CT_325_STR:
        out[0] = 1
        out[1] = 0
	elif value == VCF_CT_523_STR:
        out[0] = 0
        out[1] = 1
	elif value == VCF_CT_525_STR:
        out[0] = 0
        out[1] = 0
	elif value == VCF_CT_323_STR:
        out[0] = 1
        out[1] = 1
    return out