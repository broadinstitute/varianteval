import unittest
from varianteval.core.io.vcf import vcf_iter, vcf2callset
import varianteval.core.constants as constants

import os.path

HERE = os.path.dirname(__file__)


def test_vcf_iter():
    f = os.path.join(HERE, "test_data/GRCh37/eval/sv/HG002_hs37d5_pacbio-integrated.chr22.vcf.gz")
    for sv in vcf_iter(f):
        print(sv.sv_type)
        assert sv.sv_type in [constants.SV_TYPE.DEL, constants.SV_TYPE.INS]

def test_vcf2callset():
    f = os.path.join(HERE, "test_data/GRCh37/eval/sv/HG002_hs37d5_pacbio-integrated.chr22.vcf.gz")
    callset = vcf2callset(f)
    assert len(callset) == 644
