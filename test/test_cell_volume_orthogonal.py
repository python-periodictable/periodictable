"""Exact orthogonal invariants and neighboring-angle controls."""
from decimal import Decimal
from itertools import permutations
from math import inf, isclose, nextafter

import numpy as np
import pytest

from periodictable import formula
from periodictable.util import cell_volume


@pytest.mark.parametrize("lengths", list(permutations((2, 3, 4))))
@pytest.mark.parametrize("angles", [(), (90,), (90, 90), (90, 90, 90),
                                   (None, 90, None), (None, None, 90)])
def test_exact_orthogonal_volume(lengths, angles):
    assert cell_volume(*lengths, *angles) == 24
    assert formula("Si").volume(*lengths, *angles) == 24*1e-24


@pytest.mark.parametrize("args, keywords, expected", [
    ((2,), {}, 8), ((2, 3), {}, 12),
    ((2,), {"alpha": 90}, 8),
    ((2,), {"gamma": 90}, 8),
    ((2,), {"b": 3, "gamma": 90}, 12),
    ((2,), {"c": 4, "beta": 90}, 16),
])
def test_exact_scalar_defaults(args, keywords, expected):
    assert cell_volume(*args, **keywords) == expected


@pytest.mark.parametrize("angle", [nextafter(90, 0), nextafter(90, inf)])
@pytest.mark.parametrize("index", range(3))
def test_adjacent_right_angle(angle, index):
    angles = [90, 90, 90]
    angles[index] = angle
    # The correction is quadratic in a one-ULP deviation from pi/2.
    assert isclose(cell_volume(2, 3, 4, *angles), 24, rel_tol=5e-16)


@pytest.mark.parametrize("scalar", [float, np.float32, np.float64, np.int64])
def test_scalar_types(scalar):
    actual = cell_volume(*(scalar(x) for x in (2, 3, 4, 90, 90, 90)))
    assert actual == 24
    # Do not return an integer merely because all lengths are integers.
    assert isinstance(actual, (float, np.floating))


@pytest.mark.parametrize("angles", [
    (0, 90, 90), (30, 60, 90), (120, 120, 120), (90, 90, 180),
])
def test_flat_cells(angles):
    assert cell_volume(2, 3, 4, *angles) == 0


@pytest.mark.parametrize("angles", [(10, 10, 90)])
def test_invalid_cells(angles):
    with pytest.raises(ValueError):
        cell_volume(2, 3, 4, *angles)
    with pytest.raises(ValueError):
        formula("Si").volume(2, 3, 4, *angles)


@pytest.mark.parametrize("args", [(), ("2", 3, 4), (Decimal(2), 3, 4),
                                  (2, 3, 4, "90")])
def test_incompatible_parameters(args):
    with pytest.raises(TypeError):
        cell_volume(*args)


@pytest.mark.parametrize("lengths", [(1e308, 1e-308, 1),
                                     (1e154, 1e154, 1e-308),
                                     (1e-100, 1e-100, 1e-100),
                                     (1e100, 1e100, 1e100)])
def test_extreme_orthogonal_length_product(lengths):
    # Preserve the existing left-to-right product, not a general range redesign.
    expected = lengths[0]*lengths[1]*lengths[2]*1.0
    assert cell_volume(*lengths) == expected
    assert cell_volume(*lengths, 90, 90, 90) == expected


# cos(60) = 1/2 gives Gram determinant 1/2 and V = 12 sqrt(2).
# Mixed cases put each unsigned type in every operand position.
REAL_ANGLE_TYPES = [np.uint8, np.uint16, np.uint32, np.uint64,
                    np.int8, np.int16, np.int32, np.int64]
REAL_ANGLE_CASES = [(dtype, dtype, dtype) for dtype in REAL_ANGLE_TYPES]
REAL_ANGLE_CASES += [(dtype, float, np.int64) for dtype in REAL_ANGLE_TYPES[:4]]
REAL_ANGLE_CASES += [(float, dtype, np.int64) for dtype in REAL_ANGLE_TYPES[:4]]
REAL_ANGLE_CASES += [(float, np.int64, dtype) for dtype in REAL_ANGLE_TYPES[:4]]


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("dtypes", REAL_ANGLE_CASES)
def test_real_angle_sixty_analytic(dtypes):
    angles = [dtype(60) for dtype in dtypes]
    # Decimal square root is independent of the production trig identity.
    expected = float(Decimal(12)*Decimal(2).sqrt())
    assert cell_volume(2, 3, 4, *angles) == pytest.approx(expected, rel=2e-15)
    assert formula("Si").volume(2, 3, 4, *angles) == pytest.approx(
        expected*1e-24, rel=2e-15)


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("dtypes", REAL_ANGLE_CASES)
def test_real_angle_ninety_no_warnings(dtypes):
    angles = [dtype(90) for dtype in dtypes]
    assert cell_volume(2, 3, 4, *angles) == 24
    assert cell_volume(2, alpha=angles[0]) == 8
    assert formula("Si").volume(2, 3, 4, *angles) == 24*1e-24
