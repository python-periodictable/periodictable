import pytest

from periodictable import formula, mix_by_volume, mix_by_weight, neutron_sld


@pytest.mark.parametrize("mixer", [mix_by_weight, mix_by_volume])
@pytest.mark.parametrize("keyword", ["density", "natural_density"])
@pytest.mark.parametrize("density", [0., 0.7])
def test_explicit_mixture_density(mixer, keyword, density):
    water = formula("H2O", density=1.)
    heavy_water = formula("D2O", natural_density=1.)
    mixture = mixer(water, 1., heavy_water, 1., **{keyword: density})

    assert getattr(mixture, keyword) == pytest.approx(density)
    if density == 0.:
        assert neutron_sld(mixture) == (0., 0., 0.)
