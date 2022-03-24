from pysam import VariantFile
from varianteval.core.constants import *
from varianteval.core.variant.sv import *

# --- VCF field accessors and iteration utils for the pysam VCF library ---
def get_sv_type(pysam_vcf_record):
    return SV_TYPE_ENCODING[pysam_vcf_record.info['SVTYPE']] if 'SVTYPE' in pysam_vcf_record.info else SV_TYPE.UNK

def get_gt(pysam_vcf_record, sample_id=0):
    sample_name = pysam_vcf_record.samples[sample_id].name
    if 'GT' in pysam_vcf_record.samples[sample_name]:
        gt = pysam_vcf_record.samples[sample_name]['GT']
    else:
        gt = (None, None)
    return VCF2GT_ENCODING[gt]

def get_sv_len(pysam_vcf_record):
    if 'SVLEN' in pysam_vcf_record.info:
        return abs(int(pysam_vcf_record.info['SVLEN'] if not isinstance(pysam_vcf_record.info['SVLEN'], tuple)
                       else pysam_vcf_record.info['SVLEN'][0]))
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
