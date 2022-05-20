"""
Script to parse an annotated VCF (with annotations indicating rows that are found in comparison databases) and generate
binary overlap matrices, and upsetplots from those. Assumes vcfs that are given are generated from SURVIVOR (i.e.,
include SUPP_VEC value in the info field)
"""

import sys
import argparse
from pysam import VariantFile
import numpy as np
import pandas as pd
from upsetplot import from_memberships
from upsetplot import plot as upsetplot
from upsetplot import UpSet
import matplotlib.pyplot as plt
from matplotlib import cm
from venn import venn, pseudovenn


def annotated_overlap_mtx(annotated_vcf, ovlp_column=True):
    '''
    method to generate overlap matrix from each row's SUPP_VEC value and OVERLAP value
    '''
    outfile_path = annotated_vcf[:-4] + '.mat.txt'
    vcf_in_file = VariantFile(annotated_vcf)
    out = open(outfile_path, 'w')

    for rec in vcf_in_file.fetch():
        # get the space-separated support vector
        ovlp_line = ' '.join(rec.info['SUPP_VEC'])
        if ovlp_column:
            if 'OVERLAP' not in rec.info:
                ovlp_line += ' NONE\n'
            else:
                ovlp_line += ' ' + rec.info['OVERLAP'] + '\n'
        else:
            ovlp_line += '\n'
        out.write(ovlp_line)

    out.close()


def upset_from_annotations(overlap_df, fig_title, stacked=True):
    overlap_df = overlap_df.set_index(getattr(overlap_df, overlap_df.columns[0]) == 1)
    for col in overlap_df.columns[1:-1]:
        overlap_df.set_index(getattr(overlap_df, col) == 1, append=True, inplace=True)

    if not stacked:
        overlap_df.set_index(getattr(overlap_df, overlap_df.columns[-1]) == 1, append=True, inplace=True)

    upset = UpSet(overlap_df, show_counts=True, sort_by='cardinality', sort_categories_by=None,
                  intersection_plot_elements=10 * int(not stacked))  # min_subset_size=4
    if stacked:
        upset.add_stacked_bars(by='OVLP', colors=cm.Pastel1, elements=10)
    upset.plot()

    plt.savefig(fig_title)


def high_order_venn(overlap_df, fig_title, venn_type='full'):
    """
    Function to create venn diagrams for 4, 5, or 6 sets -- the venn package used here taken input of sets of elements
    and finds the overlap itself, rather than the matplotlib venn package which explicitly takes each intersection size
    venn_type='full','pseudo' -- flag to generate either a full venn diagram or a pseudo venn diagram
    """
    # transformation of binary overlap matrix into dict of sets – each row of the overlap matrix represents an SV so
    # the sets for the respective tools comprise the indices at which that tool column had a 1
    sets = {col.split()[0]: set([i for i in range(len(overlap_df[col])) if overlap_df[col].iloc[i] == 1]) for col in
            overlap_df.columns}
    if venn_type == 'full':
        venn(sets, legend_loc='upper left')
    else:
        pseudovenn(sets, legend_loc='upper left')
    plt.savefig(fig_title)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='VCF comparison')
    parser.add_argument('--annot_vcf', help='Comparison-annotated vcf file')
    parser.add_argument('--output_path', help='Output path for plot')
    parser.add_argument('--tools', help='Tool names', nargs='+')
    parser.add_argument('--plot_type', choices=['stacked_upset', 'upset', 'venn', 'pseudovenn'], default='stacked_upset',
                        help='Flag to indicate whether to generate an upset or a venn diagram')
    args = parser.parse_args()

    col_names = args.tools

    annotated_vcf = args.annot_vcf
    annotated_overlap_mtx(annotated_vcf, ovlp_column=(args.plot_type == 'stacked_upset'))
    if args.plot_type == 'stacked_upset':
        col_names.append('OVLP')
    overlap_df = pd.read_csv(annotated_vcf[:-4] + '.mat.txt', sep=' ', names=col_names)
    if args.plot_type in ['upset', 'stacked_upset']:
        upset_from_annotations(overlap_df, args.output_path, stacked=(args.plot_type == 'stacked_upset'))
    elif args.plot_type == 'venn':
        high_order_venn(overlap_df, args.output_path[:-4] + '_venn.png')
    else:
        high_order_venn(overlap_df, args.output_path[:-4] + '_pseudovenn.png', venn_type='pseudo')
