from cyvcf2 import VCF
import pyranges as pr


def read(file):
    for variant in VCF(file):
        print(f'{variant.CHROM} {variant.start}')