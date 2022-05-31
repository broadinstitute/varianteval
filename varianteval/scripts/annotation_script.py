"""
Script to add annotation to input callset with respect to specified
comparison callsets
"""

import sys
import re
import numpy as np
from pysam import VariantFile
import argparse
from collections import defaultdict
from intervaltree import IntervalTree


def overlapping_intervals(intrvl_A, intrvl_B, overlap_threshold=0.5, bp_margin=5000):
    '''
    Helper method to determine whether two records' start-end
    intervals overlap within a given breakpoint margin, and whether
    the intervals' overlap is >=0.5*(svlen) for interval A and B
    '''
    bp_check = (np.abs(intrvl_A[0] - intrvl_B[0]) < bp_margin \
                and np.abs(intrvl_A[1] - intrvl_B[1]) < bp_margin)
    intrvl_overlap = (min(intrvl_A[1], intrvl_B[1])) - (max(intrvl_A[0], intrvl_B[0]))
    overlap_check = (intrvl_overlap >= overlap_threshold * (intrvl_A[1] - intrvl_A[0])) and \
                    (intrvl_overlap >= overlap_threshold * (intrvl_B[1] - intrvl_B[0]))

    return bp_check and overlap_check


def build_callset_dict(callset_vcf):
    record_dict = defaultdict(IntervalTree)
    callset = VariantFile(callset_vcf)
    for rec in callset.fetch():
        record_dict[rec.chrom].addi(rec.start, rec.stop, rec.info['SVTYPE'])
    return record_dict


def annotate_source_vcf(source_vcf, annotated_vcf, comparison_vcfs, comparison_labels, overlap_threshold=0.5,
                        bp_margin=5000, type_constraint='agnostic'):
    """
    Inputs
        source_vcf: vcf that be augmented with annotations
        annotated_vcf: output vcf path
        comparison_vcfs: list of vcfs to be checked for overlap with records in source_vcf
        comparison_labels: list of labels to be used in annotation (corresponding to comparison_vcfs list)
        overlap_threshold: threshold for labeled two records as overlapping
    Output
        Annotated vcf placed in the specified directory with records' INFO fields amended with labels indicating
        overlap with the input comparison vcfs
    """
    assert len(comparison_vcfs) == len(comparison_labels), \
        'comparison_vcfs and comparison_labels do not have matching length'

    vcf_in_file = VariantFile(source_vcf)
    header = vcf_in_file.header
    # OVERLAP may be redundant with the vcf-specific labels, but also gives a compact field for just showing
    # the callset overlaps without the extra information about the conflicting record
    if 'OVERLAP' not in header.info:
        header.info.add('OVERLAP', number=1, type='String', description="Callset also containing this record (the "
                                                                        "source of the overlap)")
    if 'OVERLAPINFO' not in header.info:
        header.info.add('OVERLAPINFO', number=1, type='String',
                        description="Annotation info about overlapping comparison records")
    vcf_annt_file = VariantFile(annotated_vcf, 'w', header=header)

    # dictionary hierarchy for the labels and comparison vcfs:
    # {label: {chrom: [records]}}, for [records] stored in an intervaltree with data giving svtype of the event
    comparison_dicts = {}
    for i in range(len(comparison_vcfs)):
        comparison_dicts[comparison_labels[i]] = build_callset_dict(comparison_vcfs[i])

    for rec in vcf_in_file.fetch():
        for comp_label in comparison_labels:
            # comparison_dicts[comp_label][rec.chrom] an interval tree giving the typed SV intervals from the
            # comparison vcf corresponding to comp_label whose chromosome matches the current source record's
            for comp_start, comp_stop, comp_svtype in comparison_dicts[comp_label][rec.chrom].overlap(rec.start,
                                                                                                      rec.stop):
                # logic for checking SV type under the specified type constraint
                type_constraint_met = (type_constraint == 'agnostic') or \
                                      (type_constraint == 'same' and rec.info['SVTYPE'] == comp_svtype) or \
                                      (type_constraint == 'different' and rec.info['SVTYPE'] != comp_svtype)
                if type_constraint_met and overlapping_intervals((rec.start, rec.stop), (comp_start, comp_stop),
                                                                 overlap_threshold=overlap_threshold,
                                                                 bp_margin=bp_margin):
                    # comp_start, comp_stop = comp_intrvl
                    comp_overlap = (min(rec.stop, comp_stop)) - (max(rec.start, comp_start))
                    comp_dist_from_start = np.abs(rec.start - comp_start)
                    annot_str = "_%s-%d-%d-%d-%d-%d-%s" % (comp_label, comp_start, comp_stop,
                                                           comp_stop - comp_start, comp_overlap,
                                                           comp_dist_from_start, comp_svtype)
                    if 'OVERLAPINFO' in rec.info:
                        rec.info['OVERLAPINFO'] += annot_str
                    else:
                        rec.info['OVERLAPINFO'] = annot_str
                    if 'OVERLAP' not in rec.info:
                        rec.info['OVERLAP'] = comp_label
                    elif comp_label not in rec.info['OVERLAP']:
                        rec.info['OVERLAP'] += ('_' + comp_label)
        vcf_annt_file.write(rec)
    vcf_in_file.close()
    vcf_annt_file.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='VCF annotation')
    parser.add_argument('--output_vcf', type=str, default='./annotated_callset.vcf', help='output path for annotated '
                                                                                          'file')
    parser.add_argument('--source_vcf', type=str, help='input vcf that will have annotations add to (a copy of) it')
    parser.add_argument('--comparison_panel', type=str, help='Tab-separated panel file giving label/vcf pairs for the'
                                                             'comparison vcfs in each line of the text file')
    parser.add_argument('--overlap_thr', type=float, default=0.8,
                        help='Overlap threshold to be applied when looking for overlapping calls of different type')
    parser.add_argument('--overlap_bp_margin', type=int, default=5000,
                        help='breakpoint margin to be applied in determining overlap of two intervals')
    parser.add_argument('--allowed_overlap_type', choices=['agnostic', 'same', 'different'], default='agnostic',
                        help='SV type constraints on what events are considered overlapping: if \'agnostic\' then '
                             'overlaps will be allowed across type; if \'same\' then overlaps will require events'
                             'of the same type; if \'different\' then overlaps will require events of different type')
    args = parser.parse_args()

    # constructing the list of labels and vcf paths from the comparison panel file
    comparison_vcfs, comparison_labels = [], []
    with open(args.comparison_panel) as f:
        for line in f:
            comparison_vcfs.append(line.split()[1])
            comparison_labels.append(line.split()[0])

    annotate_source_vcf(args.source_vcf, args.output_vcf, comparison_vcfs, comparison_labels,
                        overlap_threshold=float(args.overlap_thr), bp_margin=args.overlap_bp_margin,
                        type_constraint=args.allowed_overlap_type)
