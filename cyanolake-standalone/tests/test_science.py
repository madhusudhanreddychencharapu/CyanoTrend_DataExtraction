import numpy as np
import pytest

from cyanolake.l2gen import l2_products
from cyanolake.science import additional_olci_indices, ci_family, ndci_index


def test_detection_gates_and_negative_input():
    a = ci_family([0.01], [0.03], [0.02], [0.04])
    assert a["ci_detection"][0]
    assert a["SS_681"][0] < 0 and a["SS_665"][0] > 0
    assert a["CI_cyano"][0] == pytest.approx(0.01363636, rel=1e-5)
    negative = ci_family([-0.01], [0.03], [0.02], [0.04])
    assert negative["finite"][0] and negative["ci_detection"][0]
    absent = ci_family([0.01], [0.03], [0.05], [0.04])
    assert not absent["ci_candidate"][0] and np.isnan(absent["CI_cyano"][0])
    wrong_shape = ci_family([0.04], [0.01], [0.005], [0.04])
    assert wrong_shape["ci_candidate"][0] and not wrong_shape["ci_detection"][0]


def test_strict_detection_threshold_and_nan():
    result = ci_family([0.01], [0.03], [0.02], [0.04])
    equal = ci_family([0.01], [0.03], [0.02], [0.04], detection_limit=float(result["CI"][0]))
    assert not equal["ci_detection"][0]
    bad = ci_family([np.nan], [0.03], [0.02], [0.04])
    assert not bad["finite"][0] and np.isnan(bad["CI"][0])


def test_ndci_zero_denominator():
    values, valid = ndci_index(np.array([0.02, 0.1, np.nan]), np.array([0.04, -0.1, 0.2]))
    assert values[0] == pytest.approx(1 / 3)
    assert valid.tolist() == [True, False, False]
    assert np.isnan(values[1:]).all()


def test_band_mapping_and_added_indices():
    products, required = l2_products("full_swath_rhos", False)
    assert "rhos_884" in products and "rhos_885" not in products
    assert {
        "rhos_490",
        "rhos_560",
        "rhos_620",
        "rhos_665",
        "rhos_681",
        "rhos_709",
        "rhos_754",
        "rhos_865",
        "rhos_884",
    } <= required
    arrays = {
        w: np.array([v])
        for w, v in {
            490: 0.01,
            560: 0.02,
            620: 0.01,
            665: 0.03,
            681: 0.02,
            709: 0.04,
            754: 0.05,
            865: 0.02,
            885: 0.01,
        }.items()
    }
    waves = {w: float(884 if w == 885 else w) for w in arrays}
    values = additional_olci_indices(arrays, waves)
    assert values["FAI"][0] == pytest.approx(0.02445, rel=1e-5)
    assert values["MPH_peak_nm"][0] == 754
