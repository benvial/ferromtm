from ferromtm.models import compute
import numpy.testing as npt


def test_run():
    res = compute.run(3)
    npt.assert_almost_equal(res, 9, decimal=3)
