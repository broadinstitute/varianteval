from varianteval import reader

import os.path

HERE = os.path.dirname(__file__)


def test_read_vcf():
    f = os.path.join(HERE, "test_data/small/HG002_GRCh37_1_22_v4.2.1_benchmark.chr22.vcf.gz")

    reader.read(f)
    assert 1 == 1