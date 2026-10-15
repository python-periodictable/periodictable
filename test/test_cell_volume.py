"""Check lattice volumes against independent Cartesian geometry."""
from itertools import permutations
from math import atan2, degrees, dist, inf, nextafter, prod, sqrt

import pytest

from periodictable import formula
from periodictable.util import cell_volume


def _dot(u, v):
    return sum(x*y for x, y in zip(u, v))


def _cross(u, v):
    return (u[1]*v[2] - u[2]*v[1],
            u[2]*v[0] - u[0]*v[2],
            u[0]*v[1] - u[1]*v[0])


def _parameters(vectors):
    lengths = [dist(v, (0, 0, 0)) for v in vectors]
    # atan2 avoids losing the small angles before cell_volume is called.
    angles = [degrees(atan2(dist(_cross(u, v), (0, 0, 0)), _dot(u, v)))
              for u, v in ((vectors[1], vectors[2]),
                           (vectors[0], vectors[2]),
                           (vectors[0], vectors[1]))]
    return (*lengths, *angles)


@pytest.mark.parametrize("shear", [0, 10, 100, 1000, 10000, 100000])
@pytest.mark.parametrize("order", list(permutations(range(3))))
def test_sheared_primitive_cell(shear, order):
    # Integer shears change the basis, but not the lattice or primitive volume.
    basis = [(1, 0, 0), (shear, 1, 0), (shear, 0, 1)]
    vectors = [basis[k] for k in order]
    expected = abs(_dot(vectors[0], _cross(vectors[1], vectors[2])))
    assert expected == 1
    parameters = _parameters(vectors)
    assert cell_volume(*parameters) == pytest.approx(expected, rel=1e-12)

    # Formula.volume converts Ang^3 to cm^3; derived density must not depend
    # on which primitive basis was used to describe the same crystal.
    crystal = formula("Si")
    volume = crystal.volume(*parameters)
    assert volume == pytest.approx(expected*1e-24, rel=1e-12, abs=0)
    density = crystal.molecular_mass/volume
    assert density == pytest.approx(crystal.molecular_mass/1e-24, rel=1e-12)


@pytest.mark.parametrize("vectors", [
    [(2, 0, 0), (0, 3, 0), (0, 0, 4)],
    [(2, 1, 0), (-1, 3, 1), (1, -2, 4)],
    [(1, 0, 0), (-10000, 1, 0), (0, 0, 2)],
])
def test_cartesian_volume(vectors):
    expected = abs(_dot(vectors[0], _cross(vectors[1], vectors[2])))
    assert cell_volume(*_parameters(vectors)) == pytest.approx(expected, rel=1e-10)


def test_lattice_parameter_defaults():
    assert cell_volume(2) == 8
    assert cell_volume(2, 3) == 12
    assert cell_volume(2, 3, 4) == prod((2, 3, 4))
    # Supplying alpha alone describes equal rhombohedral angles.
    assert cell_volume(2, alpha=60) == pytest.approx(8/2**0.5)
    assert cell_volume(2, alpha=60, beta=90) == pytest.approx(
        cell_volume(2, 2, 2, 60, 90, 60))
    with pytest.raises(TypeError, match="missing lattice parameters"):
        cell_volume()
    with pytest.raises(ValueError):
        cell_volume(1, alpha=10, beta=10, gamma=90)


@pytest.mark.parametrize("angles", [
    (1, 122, 123), (1, 38, 39), (30, 60, 90),
    (120, 120, 120), (90, 90, 180),
])
@pytest.mark.parametrize("order", list(permutations(range(3))))
def test_exact_planar_cell(angles, order):
    # Exact degree relationships describe coplanar vectors: volume is zero.
    assert cell_volume(1, 1, 1, *(angles[k] for k in order)) == 0


@pytest.mark.parametrize("angles", [
    (1, 122, nextafter(123, inf)),
    (1, 38, nextafter(39, inf)),
    (30, 60, nextafter(90, inf)),
    (120, 120, nextafter(120, inf)),
])
@pytest.mark.parametrize("order", list(permutations(range(3))))
def test_just_outside_planar_cell(angles, order):
    # Do not hide invalid geometry by clamping small negative determinants.
    with pytest.raises(ValueError):
        cell_volume(1, 1, 1, *(angles[k] for k in order))


@pytest.mark.parametrize("parameters, expected", [
    ((1e308, 1e-308, 1), 1),
    ((1e154, 1e154, 1e-308), 1),
    ((1e308, 1e-308, 1, 60, 60, 60), 1/sqrt(2)),
])
def test_representable_volume_without_early_overflow(parameters, expected):
    assert cell_volume(*parameters) == pytest.approx(expected, rel=1e-15)
