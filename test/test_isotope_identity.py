import pytest

from periodictable import C, D, Fe, H, formula


@pytest.mark.parametrize("text, atom", [
    ("H2O", H), ("D2O", D), ("C6H12O6", C), ("Fe[57]2O3", Fe[57]),
])
@pytest.mark.parametrize("portion", [0., 0.25, 0.5, 1.])
def test_replace_same_atom_preserves_material(text, atom, portion):
    original = formula(text, density=2.)
    result = original.replace(atom, atom, portion=portion)

    assert result is not original
    assert result.atoms == original.atoms
    assert result.mass == original.mass
    assert result.density == original.density
    result.density = 3.
    assert original.density == 2.
