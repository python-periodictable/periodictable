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
            unc = f"0.{'0' * zeros}{unc}"
        return float(value), float(unc)

    # Plain value with no uncertainty
    return float(s), 0

def from_subscript(value: str) -> str:
    """
    Convert unicode subscript characters to normal characters. This allows us to parse,
    for example, H₂O as H2O.
    """
    codepoints = {
        '\u2080': '0', '\u2081': '1', '\u2082': '2', '\u2083': '3',
        '\u2084': '4', '\u2085': '5', '\u2086': '6', '\u2087': '7',
        '\u2088': '8', '\u2089': '9', '\u208a': '+', '\u208b': '-',
        '\u208c': '=', '\u208d': '(', '\u208e': ')',

        '\u2090': 'a', '\u2091': 'e', '\u2092': 'o', '\u2093': 'x',
        '\u2095': 'h', '\u2096': 'k', '\u2097': 'l',
        '\u2098': 'm', '\u2099': 'n', '\u209a': 'p', '\u209b': 's',
        '\u209c': 't',
    }
    return ''.join(codepoints.get(char, char) for char in str(value))

def from_superscript(value: str) -> str:
    """
    Convert unicode superscript characters to normal characters. This allows us to parse,
    for example, Ca²⁺ as Ca{2+}.
    """
    codepoints = {
        '\u2070': '0', '\u00B9': '1', '\u00B2': '2', '\u00B3': '3',
        '\u2074': '4', '\u2075': '5', '\u2076': '6', '\u2077': '7',
        '\u2078': '8', '\u2079': '9', '\u207a': '+', '\u207b': '-',
        '\u207c': '=', '\u207d': '(', '\u207e': ')',

        '\u2071': 'i', '\u207f': 'n',
    }
    return ''.join(codepoints.get(char, char) for char in str(value))

def unicode_subscript(value: str) -> str:
    # Unicode subscript codepoints. Note that decimal point looks okay as subscript
    codepoints = {
        '0': '\u2080', '1': '\u2081', '2': '\u2082', '3': '\u2083',
        '4': '\u2084', '5': '\u2085', '6': '\u2086', '7': '\u2087',
        '8': '\u2088', '9': '\u2089', '+': '\u208a', '-': '\u208b',
        '=': '\u208c', '(': '\u208d', ')': '\u208e',

        'a': '\u2090', 'e': '\u2091', 'o': '\u2092', 'x': '\u2093',
        'h': '\u2095', 'k': '\u2096', 'l': '\u2097',
        'm': '\u2098', 'n': '\u2099', 'p': '\u209a', 's': '\u209b',
        't': '\u209c',

        '\u2013': '\u208b', # en-dash is same as dash
        '\u2014': '\u208b', # em-dash is same as dash
    }
    return ''.join(codepoints.get(char, char) for char in str(value))

def unicode_superscript(value: str) -> str:
    # Unicode subscript codepoints. Note that decimal point looks okay as subscript
    codepoints = {
        #'.': '\u00B0',  # degree symbol looks too much like zero
        #'.': ' \u02D9',  # dot above modifier looks okay in a floating string, but risky
        #'.': ' \u0307',  # space with dot above?
        #'.': '\u22C5', # math dot operator
        '.': '\u1427',  # Canadian aboriginal extended block dot (looks good on mac)
        '2': '\u00B2', '3': '\u00B3',
        '1': '\u00B9',
        '0': '\u2070', 'i': '\u2071',
        '4': '\u2074', '5': '\u2075', '6': '\u2076', '7': '\u2077',
        '8': '\u2078', '9': '\u2079', '+': '\u207a', '-': '\u207b',
        '=': '\u207c', '(': '\u207d', ')': '\u207e', 'n': '\u207f',

        '\u2013': '\u207b', # en-dash is same as dash
        '\u2014': '\u207b', # em-dash is same as dash
    }
    return ''.join(codepoints.get(char, char) for char in str(value))


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
