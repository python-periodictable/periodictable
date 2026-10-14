from math import exp, log, nextafter, inf
from decimal import Decimal, localcontext
from itertools import permutations

import numpy as np
import pytest

import periodictable as pt
from periodictable.activation import ActivationEnvironment, Sample


@pytest.mark.parametrize('offset', [-1, 0, 1])
@pytest.mark.parametrize('exposure', [0.0, 1.0, 10000.0])
def test_coincident_capture_rates(offset, exposure):
    isotope = pt.Co[59]
    product = next(row for row in isotope.neutron_activation
                   if row.reaction == 'act' and row.daughter == 'Co-60')
    decay_rate = log(2) / product.Thalf_hrs
    flux = decay_rate / (3600e-24 * (product.thermalXS - product.thermalXS_parent))
    if offset:
        flux = nextafter(flux, 0.0 if offset < 0 else inf)

    sample = Sample('Co', mass=1)
    sample.calculate_activation(ActivationEnvironment(fluence=flux),
                                exposure=exposure, rest_times=[0.0, 10.0])
    activity = next(values for row, values in sample.activity.items()
                    if row.reaction == 'act' and row.daughter == 'Co-60')

    # Coincident depletion/removal rates have a finite two-state limit.
    capture_rate = flux * product.thermalXS * 3600e-24
    production = flux * product.thermalXS * 1e-24 / isotope.isotope * 1.6276e19
    initial = production * decay_rate * exposure * exp(-capture_rate * exposure)
    expected = [initial, initial * exp(-decay_rate * 10)]
    np.testing.assert_allclose(activity, expected, rtol=1e-12, atol=0)


@pytest.mark.parametrize('time', [0.0, 1.0, 1000.0])
def test_repeated_integrals(time):
    from periodictable.activation import _two_rate_integral, _three_rate_integral

    rate = 0.2
    np.testing.assert_allclose(_two_rate_integral(time, rate, rate),
                               time * exp(-rate*time), rtol=1e-14, atol=0)
    np.testing.assert_allclose(_three_rate_integral(time, rate, rate, rate),
                               time**2 * exp(-rate*time) / 2, rtol=1e-14, atol=0)


@pytest.mark.parametrize('rates', [
    (1., 1.+1e-12, 1.+2e-12),
    (0., 0.004, 0.009),
    (0., 0.004, 0.010),
    (0., 0.004, 0.011),
    (0.2, 1.2, 3.),
])
def test_three_distinct_rates(rates):
    from periodictable.activation import _three_rate_integral

    # Evaluate the original divided difference at independent high precision.
    with localcontext() as context:
        context.prec = 100
        a, b, c = map(Decimal.from_float, rates)
        expected = float((-a).exp()/((b-a)*(c-a))
                         + (-b).exp()/((a-b)*(c-b))
                         + (-c).exp()/((a-c)*(b-c)))
    for ordering in permutations(rates):
        np.testing.assert_allclose(_three_rate_integral(1., *ordering),
                                   expected, rtol=1e-12, atol=0)
