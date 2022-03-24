from enum import Enum

# IO constants
BAM_TYPE = Enum('BAM_TYPE', 'SHORT LONG LINKED')
BED_TYPE = Enum("BED_TYPE", 'BED BEDPE')

# SV-related constants
class GT(Enum):
    HOM_REF = 0
    HET = 1
    HOM_ALT = 2
    UNK = 3


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
