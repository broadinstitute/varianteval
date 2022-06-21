from enum import Enum


# IO constants
BAM_TYPE = Enum("BAM_TYPE", "SHORT LONG LINKED")
BED_TYPE = Enum("BED_TYPE", "BED BEDPE")


# SV-related constants
GT = Enum("GT", "HOM_REF HET HOM_ALT UNK")
SV_TYPE = Enum("SV_TYPE", "DEL DEL_ME DEL_INV INS INS_ME INS_NOVEL DUP DUP_TANDEM DUP_INT INV INV_DUP CNV BND TRA TRA_BALANCED TRA_UNBALANCED CHROMOTHRIPSIS CHROMOPLEXY BFB DOUBLEMINUTE CPX UNK")
SV_TYPE_ENCODING = { "DEL": SV_TYPE.DEL,
                     "DEL:ME": SV_TYPE.DEL_ME,
                     "DEL/INV": SV_TYPE.DEL_INV,
                     "INS": SV_TYPE.INS,
                     "INS:ME": SV_TYPE.INS_ME,
                     "INS:NOVEL": SV_TYPE.INS_NOVEL,
                     "DUP": SV_TYPE.DUP,
                     "DUP:TANDEM": SV_TYPE.DUP_TANDEM,
                     "DUP:INT": SV_TYPE.DUP_INT,
                     "INV": SV_TYPE.INV,
                     "INVDUP": SV_TYPE.INV_DUP,
                     "CNV": SV_TYPE.CNV,
                     "BND": SV_TYPE.BND,
                     "TRA": SV_TYPE.TRA,
                     "TRA:BALANCED": SV_TYPE.TRA_BALANCED,
                     "TRA:UNBALANCED": SV_TYPE.TRA_UNBALANCED,
                     "CHROMOTHRIPSIS": SV_TYPE.CHROMOTHRIPSIS,
                     "CHROMOPLEXY": SV_TYPE.CHROMOPLEXY,
                     "BFB": SV_TYPE.BFB,
                     "DOUBLEMINUTE": SV_TYPE.DOUBLEMINUTE,
                     "CPX": SV_TYPE.CPX,
                     "UNK": SV_TYPE.UNK 
                   }
VCF2GT_ENCODING = {(0, 0): GT.HOM_REF, (0, 1): GT.HET, (1, 1): GT.HOM_ALT, (1, 0): GT.HET, (None, None): GT.UNK}
BED2GT_ENCODING = {"0/0": GT.HOM_REF, "0/1": GT.HET, "1/1": GT.HOM_ALT, "1/0": GT.HET, "./.": GT.UNK}
GT2BED_ENCODING = {GT.HOM_REF: "0/0", GT.HOM_ALT: "1/1", GT.HET: "0/1", GT.UNK: "./."}
GT2VCF_ENCODING = {GT.HOM_REF: (0, 0), GT.HOM_ALT: (1, 1), GT.HET: (0, 1), GT.UNK: (None, None)}

STRAND = Enum("STRAND", "FWD REV")
STRAND_ENCODING = {"+": STRAND.FWD, "-": STRAND.REV}


# Genome annotation constants
REPEAT_TYPE = Enum("REPEAT_TYPE", "SD SR SINE LINE MS VNTR ALL")


# Tools supported by the framework
SV_CALLER_TYPE = Enum("SV_CALLER_TYPE", "GENERIC TRUTH_SET SNIFFLES PBSV CUTESV MANTA DELLY LUMPY")
LONG_READ_CALLERS = {SV_CALLER_TYPE.SNIFFLES, SV_CALLER_TYPE.PBSV, SV_CALLER_TYPE.CUTESV}
SHORT_READ_CALLERS = {SV_CALLER_TYPE.MANTA, SV_CALLER_TYPE.DELLY, SV_CALLER_TYPE.LUMPY}


# Representation of variants
BREAKPOINT_SIDE = Enum("BREAKPOINT_SIDE", "LEFT RIGHT UNKNOWN")
PROBABILITY_FUNCTION = Enum("PROBABILITY_FUNCTION", "DELTA UNIFORM GAUSSIAN UNKNOWN")


# Basic VCF constants
VCF_FILTER_PASS = "PASS"
VCF_PRECISE_STR = "PRECISE"
VCF_IMPRECISE_STR = "IMPRECISE"
VCF_END_STR = "END"
VCF_SVTYPE_STR = "SVTYPE"
VCF_SVLEN_STR = "SVLEN"
VCF_CHR2_STR = "CHR2"
VCF_CT_STR = "CT"
VCF_CT_325_STR = "'3to5'"
VCF_CT_523_STR = "'5to3'"
VCF_CT_525_STR = "'5to5'"
VCF_CT_323_STR = "'3to3'"
VCF_INSERTION_STRING_NAME = "ins"
VCF_INSERTION_ALT = "<INS>"
VCF_GT_STR = "GT"

# Confidence intervals of positions.
#
# Remark: some callers report a standard deviation instead of a confidence interval. Sniffles reports additional
# interval information in its BEDPE output. Some callers use CILEN to express a "confidence interval around
# inserted/deleted material between breakends": we interpret CILEN exactly like CIEND, and we ignore it for insertions
# (since representing variable-length insertion strings complicates our code).
#
VCF_CI_SEPARATOR = ","
VCF_CIPOS_STR = "CIPOS"
VCF_CIEND_STR = "CIEND"
VCF_STD_START = ["STD_quant_start", "STD_POS1"]
VCF_STD_END = ["STD_quant_stop", "STD_POS2"]
VCF_CILEN_STR = "CILEN"
