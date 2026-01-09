# This program is in the public domain
# Author: Paul Kienzle
"""
Helper functions
"""
from math import sqrt

def parse_uncertainty(s: str) -> tuple[float, float]|tuple[None, None]:
    """
    Given a floating point value plus uncertainty return the pair (val, unc).

    Format is val, val(unc), [nominal] or [low,high].

    The val(unc) form is like 23.0035(12), but also 23(1), 23.0(1.0), or 
    maybe even 23(1.0). This parser does not handle exponential notation
    such as 1.032(4)E10

    The nominal form has zero uncertainty, as does a bare value.

    The [low,high] form is assumed to be a rectangular distribution of 1-sigma
    equivalent width (high-low)/sqrt(12).

    An empty string is returned as None,None rather than 0,inf.
    """
    if s == "": # missing
        # TODO: maybe 0 +/- inf ?
        return None, None

    # Parse [nominal] or [low,high]
    if s.startswith('['):
        s = s[1:-1]
        parts = s.split(',')
        if len(parts) > 1:
            low, high = float(parts[0]), float(parts[1])
            # Use equivalent 1-sigma width for a rectangular distribution
            return (high+low)/2, (high-low)/sqrt(12)
        else:
            return float(parts[0]), 0

    # Parse value(unc) with perhaps '#' at the end
    parts = s.split('(')
    if len(parts) > 1:
        # Split the value and uncertainty.
        value, unc = parts[0], parts[1].split(')')[0]
        # Count digits after the decimal for value and produce
        # 0.00...0{unc} with the right number of zeros.
        # e.g., 23.0035(12) but not 23(1) or 23.0(1.0) or 23(1.0)
        if '.' not in unc and '.' in value:
            zeros = len(value.split('.')[1]) - len(unc)
            unc = "0." + ("0"*zeros) + unc
        return float(value), float(unc)

    # Plain value with no uncertainty
    return float(s), 0

def cell_volume(a=None, b=None, c=None, alpha=None, beta=None, gamma=None) -> float:
    r"""
    Compute cell volume from lattice parameters.

    :Parameters:
        *a*, *b*, *c* : float | |Ang|
            Lattice spacings.  *a* is required.
            *b* and *c* default to *a*.
        *alpha*, *beta*, *gamma* : float | |deg|
            Lattice angles.  *alpha* defaults to 90\ |deg|.
            *beta* and *gamma* default to *alpha*.

    :Returns:
        *V* : float | |Ang^3|
            Cell volume

    :Raises:
        *TypeError* : missing or invalid parameters

    The following formula works for all lattice types:

    .. math::

        V = a b c \sqrt{1 - \cos^2 \alpha - \cos^2 \beta - \cos^2 \gamma
                          + 2 \cos \alpha \cos \beta \cos \gamma}
    """
    from math import cos, radians, sqrt
    if a is None:
        raise TypeError('missing lattice parameters')
    if b is None:
        b = a
    if c is None:
        c = a
    calpha = cos(radians(alpha)) if alpha is not None else 0
    cbeta = cos(radians(beta)) if beta is not None else calpha
    cgamma = cos(radians(gamma)) if gamma is not None else calpha
    V = a*b*c*sqrt(1 - calpha**2 - cbeta**2 - cgamma**2 + 2*calpha*cbeta*cgamma)
    return V
