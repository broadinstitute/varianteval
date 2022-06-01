"""
Script to write IGV batch scripts with instructions to screenshot each record from an input vcf
"""
import numpy as np
from pysam import VariantFile
import argparse


def generate_script(input_vcf, bam_path, output_batch_script_path, igv_screenshot_dir, min_svlen,
                    max_svlen, trg_supp_vec=None, genome='hg38', colorby_ins_size=True,
                    groupby_pair_orientation=True, viewaspairs=True, screenshot_margin=2000):
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Batch script generation')
    parser.add_argument('--input_vcf', help='Input vcf file')
    parser.add_argument('--bam_path', help='Path to BAM to be given to IGV')
    parser.add_argument('--output_dir', help='Output directory where support vector sub-directories will be generated')
    parser.add_argument('--igv_screenshot_dir', help='Screenshot directory to be specified in batch script - useful'
                                                     'when iterating over multiple sets of screenshots')
    parser.add_argument('--min_svlen', type=int, default=0, help='Minimum sv length to include in batch script')
    parser.add_argument('--max_svlen', type=int, default=250000, help='Maximum sv length to include in batch script')
    parser.add_argument('--ref_genome', type=str, default='hg38', help='Reference genome to load into IGV')
    # flags for different IGV visualization settings
    parser.add_argument('--colorby_insert_size', type=bool, default=True, help='Adds IGV colorby INSERT SIZE')
    parser.add_argument('--groupby_pair_orientation', type=bool, default=True, help='Adds groupby PAIR ORIENTATION')
    parser.add_argument('--viewaspairs', type=bool, default=True, help='Adds viewaspairs')
    parser.add_argument('--screenshot_margin', type=int, default=2000, help='Margin to be added to either side of event'
                                                                            'interval for IGV screenshot (set in bp)')
    args = parser.parse_args()

    if args.igv_screenshot_dir is None:
        igv_screenshot_dir = args.output_dir
    else:
        igv_screenshot_dir = args.igv_screenshot_dir

    generate_script(args.input_vcf, args.bam_path, args.output_dir + '/IGV_batch_script.txt',
                    igv_screenshot_dir, args.min_svlen, args.max_svlen, genome=args.ref_genome,
                    colorby_ins_size=args.colorby_insert_size,
                    groupby_pair_orientation=args.groupby_pair_orientation,
                    viewaspairs=args.viewaspairs, screenshot_margin=args.screenshot_margin)