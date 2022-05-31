'''
Script to write IGV batch scripts with instructions to screenshot each record from an input vcf with optional
filtering/iteration over SURVIVOR SUPP_VEC fields
'''
import numpy as np
import sys
from pysam import VariantFile
import os
import argparse


def generate_script_for_supp_vec(input_vcf, bam_path, output_batch_script_path, igv_screenshot_dir, min_svlen,
                                 max_svlen, trg_supp_vec=None, genome='hg38', colorby_ins_size=True,
                                 groupby_pair_orientation=True, viewaspairs=True):
    """
    method to generate a batch script for an input vcf
    """
    out = open(output_batch_script_path, 'w')
    out.write('new\n')
    out.write(f'genome {genome}\n')
    out.write('preference SAM.MAX_VISIBLE_RANGE 1000\n')
    out.write('preference SAM.SHOW_MISMATCHES FALSE\n')
    out.write(f'load {bam_path}\n')

    write_first_time = True
    input_vcf_file = VariantFile(input_vcf)
    for rec in input_vcf_file.fetch():
        svlen = rec.info['SVLEN']
        if isinstance(svlen, tuple):
            svlen = int(np.abs(rec.info['SVLEN'][0]))
        if min_svlen <= svlen <= max_svlen:
            # setting the margin on either side of the interval to svlen/10
            screenshot_margin = svlen // 10
            # if a target supp vec is given, then skip this record if it doesn't match
            if trg_supp_vec:
                if rec.info['SUPP_VEC'] != trg_supp_vec:
                    continue
            start_pos = str(rec.start - screenshot_margin)
            end_pos = str(rec.stop + screenshot_margin)
            svtype = rec.info['SVTYPE']
            out.write('goto ' + rec.chrom + ':' + start_pos + '-' + end_pos + '\n')
            if write_first_time:
                if colorby_ins_size:
                    out.write('colorby INSERT_SIZE\n')
                if groupby_pair_orientation:
                    out.write('group PAIR_ORIENTATION\n')
                if viewaspairs:
                    out.write('viewaspairs\n')
                out.write('maxPanelHeight 1000\n')
                out.write('snapshotDirectory ' + igv_screenshot_dir + '\n')
                write_first_time = False
            out.write('collapse\n')
            out.write(f'snapshot {rec.chrom}_{start_pos}_{end_pos}_{svtype}.png\n')
    out.close()


def input_vcf_sweep(input_vcf, bam_path, output_dir, min_svlen, max_svlen, genome='hg38',
                    colorby_ins_size=True, groupby_pair_orientation=True,
                    viewaspairs=True):
    """
    sweep through input vcf creating subdirectories for each support vector
    and populate the subdirectory with the corresponding batch script
    """
    input_vcf_file = VariantFile(input_vcf)
    support_vectors = set()
    for rec in input_vcf_file.fetch():
        support_vectors.add(rec.info['SUPP_VEC'])
    for supp_vec in support_vectors:
        supp_vec_subdir = output_dir + '/' + supp_vec
        os.makedirs(supp_vec_subdir)
        generate_script_for_supp_vec(input_vcf, bam_path, supp_vec_subdir + '/IGV_batch_script.txt', supp_vec_subdir,
                                     min_svlen, max_svlen, trg_supp_vec=supp_vec, genome=genome,
                                     colorby_ins_size=colorby_ins_size,
                                     groupby_pair_orientation=groupby_pair_orientation,
                                     viewaspairs=viewaspairs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Batch script generation')
    parser.add_argument('--input_vcf', help='Input vcf file')
    parser.add_argument('--bam_path', help='Path to BAM to be given to IGV')
    parser.add_argument('--output_dir', help='Output directory where support vector sub-directories will be generated')
    parser.add_argument('--trg_supp_vec', help='Target support vector to define subset of records that '
                                               'will be included in batch script (NA to be used in '
                                               'cases when the records dont have a SUPP_VEC field)')
    parser.add_argument('--igv_screenshot_dir', help='Screenshot directory to be specified in batch script - useful'
                                                     'when iterating over multiple sets of screenshots')
    parser.add_argument('--min_svlen', type=int, default=0, help='Minimum sv length to include in batch script')
    parser.add_argument('--max_svlen', type=int, default=250000, help='Maximum sv length to include in batch script')
    parser.add_argument('--ref_genome', type=str, default='hg38', help='Reference genome to load into IGV')
    parser.add_argument('--mode', choices=['non_survivor', 'survivor_single', 'survivor_multiple'],
                        default='non_survivor',
                        help='Different modes for generating a batch script: non_survivor - input VCF doesnt include'
                             'SUPP_VEC field and all records will be added to the batch script; survivor_single - input'
                             'VCF includes SUPP_VEC field and only the records with SUPP_VEC==args.trg_supp_vec will'
                             'be added to the batch script; survivor_multiple - input VCF includes SUPP_VEC field and'
                             'the screenshots from the records with the different SUPP_VEC values will be placed in'
                             'directories named with the corresponding SUPP_VEC binary string')
    # flags for different IGV visualization settings
    parser.add_argument('--colorby_insert_size', type=bool, default=True, help='Adds IGV colorby INSERT SIZE')
    parser.add_argument('--groupby_pair_orientation', type=bool, default=True, help='Adds groupby PAIR ORIENTATION')
    parser.add_argument('--viewaspairs', type=bool, default=True, help='Adds viewaspairs')
    args = parser.parse_args()

    if args.igv_screenshot_dir is None:
        igv_screenshot_dir = args.output_dir
    else:
        igv_screenshot_dir = args.igv_screenshot_dir

    if args.mode == 'non_survivor':
        generate_script_for_supp_vec(args.input_vcf, args.bam_path, args.output_dir + '/IGV_batch_script.txt',
                                     igv_screenshot_dir, args.min_svlen, args.max_svlen, genome=args.ref_genome,
                                     colorby_ins_size=args.colorby_insert_size,
                                     groupby_pair_orientation=args.groupby_pair_orientation,
                                     viewaspairs=args.viewaspairs)
    elif args.mode == 'survivor_single':
        assert args.trg_supp_vec is not None, '--trg_supp_vec must be specified when using survivor_single mode'
        generate_script_for_supp_vec(args.input_vcf, args.bam_path, args.output_dir + '/IGV_batch_script.txt',
                                     igv_screenshot_dir, args.min_svlen, args.max_svlen,
                                     trg_supp_vec=args.trg_supp_vec,
                                     genome=args.ref_genome, colorby_ins_size=args.colorby_insert_size,
                                     groupby_pair_orientation=args.groupby_pair_orientation,
                                     viewaspairs=args.viewaspairs)
    elif args.mode == 'survivor_multiple':
        input_vcf_sweep(args.input_vcf, args.bam_path, args.output_dir, args.min_svlen, args.max_svlen,
                        genome=args.ref_genome, colorby_ins_size=args.colorby_insert_size,
                        groupby_pair_orientation=args.groupby_pair_orientation,
                        viewaspairs=args.viewaspairs)
