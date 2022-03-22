from varianteval import record

import os.path

HERE = os.path.dirname(__file__)


def test_add():
    assert record.add(2, 3) == 5