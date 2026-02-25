##############################################################################
# Copyright (c) 2009 Trustees of the Columbia University
# in the City of New York.  All rights reserved.
#
# DANSE Diffraction group, Simon J. L. Billinge
#
# File coded by:    Pavol Juhas
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
#  * Neither the name of Columbia University nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
# TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
##############################################################################


"""
Cromer-Mann formula for calculating x-ray scattering factors.
"""

# module version
__id__ = "$Id: cromermann.py 1051 2010-01-30 01:01:43Z juhas $"

import os

import numpy
from numpy.typing import ArrayLike

from . import core


def getCMformula(symbol: str) -> "CromerMannFormula":
    """
    Obtain Cromer-Mann formula and coefficients for a specified element.

    *symbol* : string
        symbol of an element

    Return instance of CromerMannFormula.
    """
    if not _cmformulas:
        _update_cmformulas()
    return _cmformulas[symbol]


def fxrayatq(symbol: str, Q: ArrayLike, charge: int|None=None) -> ArrayLike:
    """
    Return x-ray scattering factors of an element at a given Q.

    *symbol* : string
         symbol of an element or ion, e.g., "Ca", "Ca2+"
    *Q* : float or [float] | |1/Ang|
         Q value
    *charge* : int
         ion charge, overrides any valence suffixes such as "-", "+", "3+".

    Return float or numpy array.
    """
    stol = numpy.array(Q) / (4 * numpy.pi)
    rv = fxrayatstol(symbol, stol, charge)
    return rv


def fxrayatstol(symbol: str, stol: ArrayLike, charge: int|None=None) -> ArrayLike:
    """
    Calculate x-ray scattering factors at specified sin(theta)/lambda

    *symbol* : string
        symbol of an element or ion, e.g., "Ca", "Ca2+"
    *stol* : float or [float] | |1/Ang|
        sin(theta)/lambda
    *charge* : int
        ion charge, overrides any valence suffixes such as "-", "+", "3+".

    Return float or numpy.array.
    """
    # resolve lookup symbol smbl, by default symbol
    smbl = symbol
    # build standard element or ion symbol
    if charge is not None:
        smbl = symbol.rstrip('012345678+-')
        if charge:
            smbl += ("%+i" % charge)[::-1]
    # convert Na+ or Cl- to Na1+, Cl1-
    elif symbol[-1:] in '+-' and not symbol[-2:-1].isdigit():
        smbl = (symbol[:-1] + "1" + symbol[-1:])
    # smbl is resolved here
    cmf = getCMformula(smbl)
    rv = cmf.atstol(stol)
    return rv


class CromerMannFormula:
    """
    Cromer-Mann formula for x-ray scattering factors.
    Coefficient storage and evaluation.

    Class data:

    *stollimit* : float | |1/Ang|
        maximum sin(theta)/lambda for which the formula works

    Attributes:

    *symbol* : string
        symbol of an element
    *a* : [float]
        a-coefficients
    *b* : [float]
        b-coefficients
    *c* : float
        c-coefficient
    """

    # obtained from tables/f0_WaasKirf.dat and the associated reference
    # D. Waasmaier, A. Kirfel, Acta Cryst. (1995). A51, 416-413
    # http://dx.doi.org/10.1107/S0108767394013292
    stollimit: float = 6
    a: numpy.ndarray
    b: numpy.ndarray
    c: float
    symbol: str

    def __init__(self, symbol, a, b, c):
        """
        Create a new instance of CromerMannFormula for specified element.

        No return value
        """
        self.symbol = symbol
        self.a = numpy.asarray(a, dtype=float)
        self.b = numpy.asarray(b, dtype=float)
        self.c = float(c)

    def atstol(self, stol: ArrayLike) -> ArrayLike:
        """
        Calculate x-ray scattering factors at specified sin(theta)/lambda

        *stol* : float or [float] | |1/Ang|
            sin(theta)/lambda

        Return float or numpy.array.
        """
        stolflat = numpy.array(stol).flatten()
        n = len(stolflat)
        stol2row = numpy.reshape(stolflat ** 2, (1, n))
        bcol = self.b.reshape((len(self.a), 1))
        bstol2 = numpy.dot(bcol, stol2row)
        adiag = numpy.diag(self.a)
        rvrows = numpy.dot(adiag, numpy.exp(-bstol2))
        rvflat = rvrows.sum(axis=0) + self.c
        rvflat[stolflat > self.stollimit] = numpy.nan
        # when stol is scalar, addition of zero converts the rv array to float
        rv = rvflat.reshape(numpy.shape(stol)) + 0.0
        return rv

# class CromerMannFormula


def _update_cmformulas() -> None:
    """
    Update the static dictionary of CromerMannFormula instances.
    """
    path = os.path.join(core.get_data_path('.'), 'f0_WaasKirf.dat')
    # Data file contains:
    #   #{C,D,F,UD,UO,UT,UIDL,UF0TYPE}              header lines
    #
    #   #S <Z> <El>(#[+-])?                         data blocks
    #   #N 11
    #   #L a1 a2 a3 a4 a5 c b1 b2 b3 b4 b5
    #    <a1 a2 a3 a4 a5> <c> <b1 b2 b3 b4 b5>
    with open(path) as fp:
        for line in fp:
            if line.startswith("#S"):
                symbol = line.split()[2]
            elif line.startswith(" "):
                if symbol is None:
                    raise RuntimeError(f"Should not be here; {path} has been corrupted.")
                w = line.split()
                a = list(map(float, w[0:5]))
                b = list(map(float, w[6:11]))
                c = float(w[5])
                cmf = CromerMannFormula(symbol, a, b, c)
                _cmformulas[cmf.symbol] = cmf
                symbol = None

_cmformulas: dict[str, CromerMannFormula] = {}

# End of file
