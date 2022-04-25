from enum import Enum, unique

# IO constants
BAM_TYPE = Enum('BAM_TYPE', 'SHORT LONG LINKED')
BED_TYPE = Enum('BED_TYPE', 'BED BEDPE')

# SV-related constants
@unique
class GT(IntEnum):
    HOM_REF = auto()
    HET = auto()
    HOM_ALT = auto()
    UNK = auto()


SV_TYPE = Enum('SV_TYPE', 'DEL DUP INV INS TRA BND CPX UNK')
SV_TYPE_ENCODING = {"DEL": SV_TYPE.DEL,
                    "DUP": SV_TYPE.DUP,
                    "INV": SV_TYPE.INV,
                    "INS": SV_TYPE.INS,
                    "TRA": SV_TYPE.TRA,
                    "BND": SV_TYPE.BND,
                    "CPX": SV_TYPE.CPX,
                    "UNK": SV_TYPE.UNK}
VCF2GT_ENCODING = {(0, 0): GT.HOM_REF, (0, 1): GT.HET, (1, 1): GT.HOM_ALT, (1, 0): GT.HET, (None, None): GT.UNK}
BED2GT_ENCODING = {"0/0": GT.HOM_REF, "0/1": GT.HET, "1/1": GT.HOM_ALT, "1/0": GT.HET, "./.": GT.UNK}
GT2BED_ENCODING = {GT.HOM_REF: "0/0", GT.HOM_ALT: "1/1", GT.HET: "0/1", GT.UNK: "./."}
GT2VCF_ENCODING = {GT.HOM_REF: (0, 0), GT.HOM_ALT: (1, 1), GT.HET: (0, 1), GT.UNK: (None, None)}

STRAND = Enum('STRAND', 'FWD REV')
STRAND_ENCODING = {"+": STRAND.FWD, "-": STRAND.REV}

# Genome annotation constants
REPEAT_TYPE = Enum('REPEAT_TYPE', 'SD SR SINE LINE MS VNTR ALL')

# Tools supported by the framework
SV_CALLER_TYPE = Enum('SV_CALLER_TYPE', 'GENERIC TRUTH_SET SNIFFLES PBSV CUTESV MANTA DELLY LUMPY')
LONG_READ_CALLERS = {SV_CALLER_TYPE.SNIFFLES, SV_CALLER_TYPE.PBSV, SV_CALLER_TYPE.CUTESV}
SHORT_READ_CALLERS = {SV_CALLER_TYPE.MANTA, SV_CALLER_TYPE.DELLY, SV_CALLER_TYPE.LUMPY}

# Representation of variants
BREAKPOINT_SIDE = IntEnum('BREAKPOINT_SIDE', 'LEFT RIGHT'):
PROBABILITY_FUNCTION = IntEnum('PROBABILITY_FUNCTION', 'DELTA UNIFORM GAUSSIAN')

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

# Confidence intervals of positions.
#
# Remark: some callers report a standard deviation instead of a confidence
# interval. Sniffles reports additional interval information in its BEDPE
# output. Some callers use CILEN to express a "confidence interval around
# inserted/deleted material between breakends": we interpret CILEN exactly like
# CIEND, and we ignore it for insertions (since representing variable-length
# insertion strings complicates our code).
#
VCF_CI_SEPARATOR = ",";
VCF_CIPOS_STR = "CIPOS";
VCF_CIEND_STR = "CIEND";
VCF_STD_START1_STR = "STD_quant_start";
VCF_STD_START2_STR = "STD_POS1";
VCF_STD_END1_STR = "STD_quant_stop";
VCF_STD_END2_STR = "STD_POS2";
VCF_CILEN_STR = "CILEN";