# -*- coding: iso-8859-15 -*-
# This program is public domain
# Author: Paul Kienzle
# Based on spreadsheet by Les Slaback (1998).
r"""
Calculate expected neutron activation from time spent in beam line.

Accounts for burnup and 2n, g production.

This is only a gross estimate.  Many effects are not taken into account, such as
self-shielding in the sample and secondary activation from the decay products.

**Introduction to neutron activation terminology**

See a text!! These are just a few notes to refresh the memory of those who
already know the topic.

Reactions: The most common reaction is the simple addition of a neutron,
thereby increasing its atomic number by one but leaving it as the same
element. Not all activation products are radioactive. Those that are not are
not included in this database. Beware that there are production chains that
depend upon these intermediate products. This database does not address those
more complicated processes.

There are exceptions to the simple addition process such as the n,alpha
reaction of Li-6. These then result in a different element as the activation
product. These are identified in the reaction column. Radioactive products
also have the potential of undergoing a neutron reaction. This is
accounted for in this database for selected products (generally those that
result in significant half-life products). These reactions and products
should only be important at very high fluences, i.e., >1E16 n/cm2.

Neutron energy: The majority of the reactions are initiated by thermal neutrons.
As a practical matter the number of thermal neutrons is usually measured with a
cadmium filter. This excludes the neutrons above the cadmium absorption
threshold, called the cadmium cutoff. For most materials the resulting measured
'thermal' neutron fluence is adequate to determine the thermal neutron
activation. A few materials have large cross-sections for the neutrons above the
cadmium cutoff. Hence a 'cadmium ratio' is needed to predict the number of
neutrons present above this cutoff, as seen by the specific reaction of
interest.

For the same neutron spectrum two different elements will have different
cadmium ratios. Similarly, for the same reaction but different neutron
environments the ratio will also vary.  Copper is a good material on which
to predict the thermal-to-epithermal fluence ratio since its two cross-sections
are about equal.  If you use the ratio of another reaction, correct the
activity ratio of the nuclide by the cross-section ratio in order to derive
the extimated fluence ratio.  Take care.  By definition, the cadmium ratio
is an activity ratio, not a fluence ratio.  [Note: This ratio is divided
into the specfied thermal neutron fluence rate to get the epithermal rate.]

Fast neutron reactions: There are also those reactions that depend on
higher energy neutrons. These can be of a great variety: n,gamma; n,p;
n,d; n,alpha; n,triton; etc. Any they are highly variable in how high the
neutron energy must be to initiate the reaction. In general the cross-
sections are relatively small, e.g., 10's of millibarns and less as compared
to 1000's of barns in some cases for thermal reactions. Plus the fast
neutron fluences are much smaller so that in general the amount of the
activation product is much less than for thermal production processes.

Barn: As they say- just how big is the broad side of that barn? For neutrons
it is 1 E-24 cm2, e.g. 0.000000000000000000000001 square cm. Millibarns
are of course for the _______ (leprechauns?).

Calulation:
(number of atoms)*(cross-section)*(no. of neutrons)*(decay correction)*C.F.

Number of atoms: mass*avagadro's number*isotopic abundance/mole.weight
  where the mass is that of the target element, not the whole sample.

Decay correction: Since the radioactive product decays while it is being
  made one must correct for this using (1-exp(-0.693*t/hlflf))
  where t is the exposure time and hlflf is the product halflife

C.F: includes unit conversion factors, e.g., convert to microcuries


Example::

    >>> from periodictable import activation
    >>> env = activation.ActivationEnvironment(fluence=1e5, Cd_ratio=70, fast_ratio=50, location="BT-2")
    >>> sample = activation.Sample("Co30Fe70", 10)
    >>> sample.calculate_activation(env, exposure=10, rest_times=[0, 1, 24, 360])
    >>> sample.show_table()
                                           ----------------- activity (uCi) ------------------
    isotope  product   reaction  half-life        0 hrs        1 hrs       24 hrs      360 hrs
    -------- --------- -------- ---------- ------------ ------------ ------------ ------------
    Co-59    Co-60          act   5.2714 y    0.0004957    0.0004957    0.0004955     0.000493
    Co-59    Co-60m+        act   10.467 m        1.664      0.03129          ---          ---
    -------- --------- -------- ---------- ------------ ------------ ------------ ------------
                                     total        1.664      0.03181    0.0005079    0.0005045
    -------- --------- -------- ---------- ------------ ------------ ------------ ------------

    >>> print("%.3f"%sample.decay_time(0.001)) # number of hours to reach 1 nCi
    2.046

The default rest times used above show the sample activity at the end of neutron
activation and after 1 hour, 1 day, and 15 days.

The neutron activation table, *activation.dat*,\ [#Shleien1998]_ contains
details about the individual isotopes, with interaction cross sections taken
from from IAEA-273\ [#IAEA1987]_, and half-lives from NUBASE2020\ [#Kondev2021]_

Activation can be run from the command line using::

    $ python -m periodictable.activation FORMULA

where FORMULA is the chemical formula for the material.

.. [#Shleien1998] Shleien, B., Slaback, L.A., Birky, B.K., 1998.
   Handbook of health physics and radiological health.
   Williams & Wilkins, Baltimore.

.. [#IAEA1987] IAEA (1987)
   Handbook on Nuclear Activation Data.
   TR 273 (International Atomic Energy Agency, Vienna, Austria, 1987).
   http://cds.cern.ch/record/111089/files/IAEA-TR-273.pdf

.. [#Kondev2021] F.G. Kondev, M. Wang, W.J. Huang, S. Naimi, and G. Audi, 2021.
   Chin. Phys. C45, 030001. DOI:10.1088/1674-1137/abddae
   Data retrieved from https://www.anl.gov/phy/reference/nubase-2020-nubase4mas20

Code original developed for spreadsheet by Les Slaback of NIST.
"""

# Comments on ACT_2N_X.xls from ACT_CALC.TXT
#
# Fast neutron cross-sections added 2/1/94
#
# One database error detected and corrected at this time.
#
# Corrections made 3/4/94
#
# 1. Halflife for entry 21 entered.       (previous entry blank)
# 2. Cross-section for entry 311 entered  (previous entry blank)
# 3. 2 beta mode entries - added parent lambda's
# 4. Changed range of look up function in 2 columns to include last row.
# 5. With 9E15 y halflife the calculation always returned 0 microcuries!
#     The 9E15 exceeds a math precision resulting in zero being assigned to the
#     the function 1-e(-x).  For any value of x<1E-10 the approximation x+x^2/2
#     is more accurate, and for values of X<1E-16 this approximation returns a
#     non-zero value while 1-e(-x)=0.
#     Appropriate changes have been made, but in fact only two entries are
#     affected by this, one with a 12 Ty halflife and one with a 9000Ty t1/2.
#     The equations for "b" mode production were not changed because they are
#     more complex, and none of the halflives currently in the database related
#     to this mode are a problem. Take care if you add more with very long halflives.
#     [PAK: These have since been updated to use expm1()]
# 6. The cross-section for the reaction Sr-88(n,gamma)Sr-89 was reduced from
#    .058 to .0058 b.  The entry of .058 in IAEA273 is in error, based on a
#    number of other references (including IAEA156).
# 7. The unit for the halflife of Pm-151 coorected from 'm' to 'h'
#
# Changes made in April 1994
#
# 1. Burnup cross-sections added to the database and nuclides produced by
#     two neutron additions (2n,gamma) have been added.  The activation equations
#     have been changed to account for burnup.  This does not become significant
#     until exposures exceed 1e17 n/cm2 (e.g., 1e10 n/cm2/sec*1000 hrs), even for
#     those with very large cross-sections.  Note that 'burnup' can be viewed as
#     loss of the intended n,gamma product or as a 2n,gamma production mode.  Both
#     effects are included in the database and equations.
# 2.ACT_CALC.WQ1 does not have the burnup equations or related cross-section data.
#     ACT_2N.WQ1 has the added equations and data, but requires manual entry of the
#     cross-section database index numbers.
#     ACT_2N_X.WQ1 allows direct entry of the chemical element to retrieve ALL
#     entries from the database for that element
#     Results from both have been compared to assure that the new spreadsheet
#     is correct.  Also the Au-197(2n,g)Au-199 reaction has been checked against
#     a textbook example (Friedlander,Kennedy,Miller:Nuclear Chemistry).
# 3.An educational note related to this addition:
# Computing the equation [exp(-x)-exp(-y)] in that syntax is better than using
# [exp(-x)*(1-exp(x-y))].  The latter format blows up when large
# values of X and Y are encountered due to computational limitations for these
# functions in a PC.  But this problem was only encountered in excess of
# 1E22 n/cm2.
# 4. 61 (2n,g) reactions have been added.  You will note many more
# radionuclides with secondary capture cross-sections.  Many of these go to
# stable nuclides (particularly the larger x-section ones - logically enough).
# The isomer-ground state pairs that have 2n modes to the same resulting
# nuclide are treated one of several ways.
#     - if the obvious dominate pathway is only through the ground state then
#     only production via that mode is included.
#     - if the combination of parent halflives and production cross-sections
#     make it unclear as to which is the dominant production mode then both
#     are calculated.  Beware, these are not necessarily additive values.
#     Sometimes the ground state also reflects production via the isomer.
#     See the specific notes for any particular reaction.
# 5.None of the entries for production via 'b' (beta decay ingrowth from a
# neutron induced parent) account for burnup.  Those with significant burnup
# cross-sections have a specific note reminding of this.  This is an issue
# only at high fluence rates.
# 6.Note that the database entries have different meanings for 'b' and '2n'
# production modes.  For instance, the first set of cross-sections (proceeding
# horizontally for a particular database entry) is not the cross-section of
# the 2n reaction, but that to produce the parent.  The second set of cross-
# sections which are the burnup cross-sections for n,g; n,p; n,alpha; etc.
# reactions are the production cross-sections for the 2n reactions.
# 7.If you want to determine how much burnup is occurring do one of the
# following:
#     - do the calculation separately in ACT_CALC and ACT_2N.
#     - do the calculation at a low fluence, e.g., 1E7, prorate to the fluence
#     of interest, and compare to the result calculated directly.  Make
#     certain the same exposure time is used.  You cannot prorate this
#     parameter.
#
# Additions/changes made July 1997
#
# The following n,p and n,alpha reactions have both a thermal cross-section (as
# per IAEA273) and a fast cross-section.  In all cases the database entry was for
# the fast cross-section but indicated it was a thermal reaction.  That has been
# corrected (and verified against IAEA156) and a second entry made for the thermal
# reaction.  Despite the entry in IAEA 273 there is some question as to whether
# the thermal induced reaction is energetically possible.  Dick Lindstrom's
# calculation shows that the coulomb barrier should prevent the thermal reaction
# from being possible.
#
# The entries changed, and added, are for the following:
#
#     35Cl (n,p)
#     35Cl (n,alpha)
#     33S (n,p)
#     39K (n,alpha)
#     40Ca (n,p)
#     58Ni (n,alpha)


from math import exp, log, expm1
import os
from collections.abc import Callable, Sequence
from typing import cast

from .formulas import formula as build_formula, Formula, FormulaInput
from .core import Element, Isotope, isisotope, get_data_path

LN2 = log(2)

def table_abundance(iso: Isotope) -> float:
    """
    Isotopic abundance in % from the periodic table package.
    """
    return iso.abundance if iso.abundance else 0.

def IAEA1987_isotopic_abundance(iso: Isotope) -> float:
    """
    Isotopic abundance in % from the IAEA, as provided in the activation.dat table.

    Note: this will return an abundance of 0 if there is no neutron activation for
    the isotope even though for isotopes such as H[1], the natural abundance may in
    fact be rather large.

    IAEA 273: Handbook on Nuclear Activation Data, 1987.
    """
    activation = getattr(iso, "neutron_activation", None)
    if activation is not None:
        return activation[0].abundance
    return 0.

class Sample:
    """
    Sample properties.

    *formula* : chemical formula

        Chemical formula.  Any format accepted by
        :func:`.formulas.formula` can be used, including
        formula string.

    *mass* : float | g

        Sample mass.

    *name* : string

        Name of the sample (defaults to formula).
    """
    formula: Formula
    mass: float
    name: str
    activity: dict["ActivationResult", list[float]]
    environment: "ActivationEnvironment"
    exposure: float
    rest_times: tuple[float, ...]

    def __init__(self, formula: FormulaInput, mass: float, name: str|None=None):
        self.formula = build_formula(formula)
        self.mass = mass               # cell F19
        self.name = name if name else str(self.formula) # cell F20
        self.activity = {}

        # The following are set in calculation_activation
        #self.environment = None
        self.exposure = 0.
        self.rest_times = ()

    def calculate_activation(
            self,
            environment: "ActivationEnvironment",
            exposure: float=1,
            rest_times: tuple[float, ...]=(0, 1, 24, 360),
            abundance: Callable[[Isotope], float]=table_abundance,
            ):
        """
        Calculate sample activation (uCi) after exposure to a neutron flux.

        *environment* is the exposure environment.

        *exposure* is the exposure time in hours (default is 1 h).

        *rest_times* are deactivation times in hours (default [0, 1, 24, 360]).

        *abundance* is a function that returns the relative abundance of an
        isotope.  By default it uses :func:`table_abundance` with natural
        abundance defined in :mod:`periodictable.mass`, but there is the
        alternative :func:`IAEA1987_isotopic_abundance` in the activation data
        table.
        """
        self.activity = {}
        self.environment = environment
        self.exposure = exposure
        self.rest_times = rest_times
        for el, frac in self.formula.mass_fraction.items():
            if isisotope(el):
                A = activity(cast(Isotope, el), self.mass*frac, environment, exposure, rest_times)
                self._accumulate(A)
            else:
                for iso_num in el.isotopes:
                    iso: Isotope = cast(Element, el)[iso_num]
                    iso_mass = self.mass*frac*abundance(iso)*0.01
                    if iso_mass:
                        A = activity(iso, iso_mass, environment, exposure, rest_times)
                        self._accumulate(A)

    def decay_time(self, target: float, tol: float=1e-10):
        """
        After determining the activation, compute the number of hours required to achieve
        a total activation level after decay.
        """
        if not self.rest_times or not self.activity:
            return 0

        # Find the smallest rest time (probably 0 hr)
        k, t_k = min(enumerate(self.rest_times), key=lambda x: x[1])
        # Find the activity at that time, and the decay rate
        data = [
            (Ia[k], LN2/a.Thalf_hrs)
            for a, Ia in self.activity.items()
            # TODO: not sure why Ia is zero, but it messes up the initial value guess if it is there
            if Ia[k] > 0.0
            ]
        # Need an initial guess near the answer otherwise find_root gets confused.
        # Small but significant activation with an extremely long half-life will
        # dominate at long times, but at short times they will not affect the
        # derivative. Choosing a time that satisfies the longest half-life seems
        # to work well enough.
        guess = max(-log(target/Ia)/La + t_k for Ia, La in data)
        # With times far from zero the time resolution in the exponential is
        # poor. Adjust the start time to the initial guess, rescaling intensities
        # to the activity at that time.
        adj = [(Ia*exp(-La*(guess-t_k)), La) for Ia, La in data]
        #print(adj)
        # Build f(t) = total activity at time T minus target activity and its
        # derivative df/dt. f(t) will be zero when activity is at target
        f = lambda t: sum(Ia*exp(-La*t) for Ia, La in adj) - target
        df = lambda t: sum(-La*Ia*exp(-La*t) for Ia, La in adj)
        #print("data", data, [])
        t, ft = find_root(0, f, df, tol=tol)
        percent_error = 100*abs(ft)/target
        if percent_error > 0.1:
            #return 1e100*365*24 # Return 1e100 rather than raising an error
            msg = (
                "Failed to compute decay time correctly (%.1g error). Please"
                " report material, mass, flux and exposure.") % percent_error
            raise RuntimeError(msg)
        # Return time at least zero hours after removal from the beam. Correct
        # for time adjustment we used to stablize the fit.
        return max(t+guess, 0.0)

    def _accumulate(self, activity: dict["ActivationResult", list[float]]):
        for el, activity_el in activity.items():
            el_total = self.activity.get(el, [0]*len(self.rest_times))
            self.activity[el] = [T+v for T, v in zip(el_total, activity_el)]

    def show_table(self, cutoff: float=0.0001, format: str="%.4g"):
        """
        Tabulate the daughter products.

        *cutoff=1* : float | uCi

              The minimum activation value to show.

        *format="%.1f"* : string

              The number format to use for the activation.
        """
        # TODO: need format="auto" which picks an appropriate precision based on
        # cutoff and/or activation level.

        # Track individual rows with more than 1 uCi of activation, and total activation
        # Replace any activation below the cutoff with '---'
        rows = []
        total = [0.]*len(self.rest_times)
        for el, activity_el in sorted_activity(self.activity):
            total = [t+a for t, a in zip(total, activity_el)]
            if all(a < cutoff for a in activity_el):
                continue
            activity_str = [format%a if a >= cutoff else "---" for a in activity_el]
            rows.append([el.isotope, el.daughter, el.reaction, el.Thalf_str]+activity_str)
        footer = ["", "", "", "total"] + [format%t if t >= cutoff else "---" for t in total]

        # If no significant total activation then don't print the table
        if all(t < cutoff for t in total):
            print("No significant activation")
            return

        # Print the table header, with an overbar covering the various rest times
        # Print a dashed separator above and below each column
        header = ["isotope", "product", "reaction", "half-life"] \
                 + ["%g hrs"%vi for vi in self.rest_times]
        separator = ["-"*8, "-"*9, "-"*8, "-"*10] + ["-"*12]*len(self.rest_times)
        cformat = "%-8s %-9s %8s %10s " + " ".join(["%12s"]*len(self.rest_times))

        width = sum(len(c)+1 for c in separator[4:]) - 1
        if width < 16:
            width = 16
        overbar = "-"*(width//2-8) + " activity (uCi) " + "-"*((width+1)//2-8)
        offset = sum(len(c)+1 for c in separator[:4]) - 1
        print(" "*(offset+1)+overbar)
        print(cformat%tuple(header))
        print(cformat%tuple(separator))

        # Print the significant table rows, or indicate that there were no
        # significant rows if the total is significant but none of the
        # individual isotopes
        if rows:
            for r in rows:
                print(cformat%tuple(r))
        else:
            print("No significant isotope activation")
        print(cformat%tuple(separator))

        # If there is more than one row, or if there is enough marginally
        # significant activation that the total is greater then the one row
        # print the total in the footer
        if len(rows) != 1 or any(c != t for c, t in zip(rows[0][4:], footer[4:])):
            print(cformat%tuple(footer))
            print(cformat%tuple(separator))

def find_root(
        x: float,
        f: Callable[[float], float],
        df: Callable[[float], float],
        max: int=20,
        tol: float=1e-10,
        ):
    r"""
    Find zero of a function.

    Returns when $|f(x)| < tol$ or when max iterations have been reached,
    so check that $|f(x)|$ is small enough for your purposes.

    Returns x, f(x).
    """
    fx = f(x)
    for _ in range(max):
        #print(f"step {_}: {x=} {fx=} df/dx={df(x)} dx={fx/df(x)}")
        if abs(fx) < tol:
            break
        x -= fx / df(x)
        fx = f(x)
    return x, fx


def sorted_activity(
        activity: dict["ActivationResult", list[float]],
        ) -> list[tuple["ActivationResult", list[float]]]:
    """Interator over activity pairs sorted by isotope then daughter product."""
    return sorted(activity.items(), key=lambda x: (x[0].isotope, x[0].daughter))


class ActivationEnvironment:
    """
    Neutron activation environment.

    The activation environment provides details of the neutron flux at the
    sample position.

    *fluence* : float | n/cm^2/s

        Thermal neutron fluence on sample.  For COLD neutrons enter equivalent
        thermal neutron fluence.

        **Warning**: For very high fluences, e.g., >E16 to E17 n/cm2, the
        equations give erroneous results because of the precision limitations.
        If there is doubt simply do the calculation at a lower flux and
        proportion the result. This will not work for the cascade reactions,
        i.e., two neutron additions.

    *Cd_ratio* : float

        Neutron cadmium ratio.  Use 0 to suppress epithermal contribution.

        This is to account for those nuclides that have a significant
        contribution to the activation due to epithermal neutrons. This is
        tabulated in the cross-section database as the 'resonance cross-
        section'. Values can range from 4 to more than 100.

        ....... Use 0 for your initial calculation ..........

        If you do a specific measurement for the nuclide and spectrum
        of interest you simply apply a correction to the thermal based
        calculation, i.e., reduce the fluence by the appropriate factor.

        This computation can be based on a Cd ratio of a material that has
        no significant resonance cross-section, or that has been corrected
        so that it reflects just the thermal to epithermal fluence ratio.
        The computation simply adds a portion of this resonance cross-section
        to the thermal cross-section based on the presumption that the
        Cd ratio reflects the fluence ratio. Copper is a good material on
        which to equate the Cd ratio to the thermal-epithermal fluence ratio.

        For other materials, correct for the thermal:epithermal cross-section
        ratio.

        At the NBSR this ranges from 12 in mid-core, to 200 at RT-4,
        to >2000 at a filtered cold neutron guide position.

    *fast_ratio* : float

        Thermal/fast ratio needed for fast reactions. Use 0 to suppress fast
        contribution.

        This is very reaction dependent. You in essence need to know the
        answer before you do the calculation! That is, this ratio depends
        upon the shape of the cross-section curve as well as the spectrum
        shape above the energy threshold of the reaction. But at least you
        can do 'what-if' calculations with worse case assumptions (in the
        absence of specific ratios).

        The fast cross-sections in this database are weighted for a fast
        maxwellian spectrum so the fast/thermal ratio should be just a
        fluence correction (i.e., a fluence ratio), not an energy correction.

        Fast neutron cross-sections from IAEA273 are manually spectrum weighted.
        Those from IAEA156 are fission spectrum averaged as tabulated.

        Use 50 (for NBSR calculations) as a starting point if you simply
        exploring for possible products.

    """
    fluence: float
    Cd_ratio: float
    fast_ratio: float
    location: str

    def __init__(self, fluence=1e5, Cd_ratio=0., fast_ratio=0., location=""):
        self.fluence = fluence     # cell F13
        self.Cd_ratio = Cd_ratio   # cell F15
        self.fast_ratio = fast_ratio # cell F17
        self.location = location   # cell F21

    # Cell Q1
    @property
    def epithermal_reduction_factor(self):
        """
        Used as a multiplier times the resonance cross section to add to the
        thermal cross section for all thermal induced reactions.
        """
        return 1./self.Cd_ratio if self.Cd_ratio >= 1 else 0

COLUMN_NAMES = [
    "_symbol",     # 0 AF
    "_index",      # 1 AG
    "Z",           # 2 AH
    "symbol",      # 3 AI
    "A",           # 4 AJ
    "isotope",     # 5 AK
    "abundance",   # 6 AL
    "daughter",    # 7 AM
    "_Thalf",      # 8 AN
    "_Thalf_unit", # 9 AO
    "isomer",      # 10 AP
    "percentIT",   # 11 AQ
    "reaction",    # 12 AR
    "fast",        # 13 AS
    "thermalXS",   # 14 AT
    "gT",          # 15 AU
    "resonance",   # 16 AV
    "Thalf_hrs",   # 17 AW
    "Thalf_str",   # 18 AX
    "Thalf_parent", # 19 AY
    "thermalXS_parent",  # 20 AZ
    "resonance_parent",  # 21 BA
    "comments",          # 22 BB
]
INT_COLUMNS = [1, 2, 4]
BOOL_COLUMNS = [13]
FLOAT_COLUMNS = [6, 11, 14, 15, 16, 17, 19, 20, 21]
UNITS_TO_HOURS = {'y': 8760, 'd': 24, 'h': 1, 'm': 1/60, 's': 1/3600}

def activity(
        isotope: Isotope,
        mass: float,
        env: ActivationEnvironment,
        exposure: float,
        rest_times: Sequence[float],
        ) -> dict["ActivationResult", list[float]]:
    """
    Compute isotope specific daughter products after the given exposure time and
    rest period.

    Activations are listed in *isotope.neutron_activation*. Most of the
    activations (n,g n,p n,a n,2n) define a single step process, where a neutron
    is absorbed yielding the daughter and some prompt radiation. The daughter
    itself will decay during exposure, yielding a balance between production and
    emission. Any exposure beyound about eight halflives will not increase
    activity for that product.

    Activity for daughter products may undergo further neutron capture, reducing
    the activity of the daughter product but introducing a grand daughter with
    its own activity.

    The data tables for activation are only precise to about three significant
    figures. Any changes to the calculations below this threshold, e.g., due to
    slightly different mass or abundance, are therefore of little concern.

    Differences in formulas compare to the NCNR activation spreadsheet:

    * Column M: Use ln(2) rather than 0.693
    * Column N: Use ln(2) rather than 0.693
    * Column O: Rewrite to use expm1(x) = exp(x) - 1::
        = L (1 - exp(-M Y43)/(1-(M/N) ) + exp(-N Y43)/((N/M)-1) )
        = L (1 - N exp(-M Y43)/(N-M) + M exp(-N Y43)/(N-M) )
        = L/(N-M) ( (N-M) - N exp(-M Y43) + M exp(-N Y43) )
        = L/(N-M) ( N(1 - exp(-M Y43)) + M (exp(-N Y43) - 1) )
        = L/(N-M) ( -N expm1(-M Y43) + M expm1(-N) )
        = L/(N-M) (M expm1(-N Y43) - N expm1(-M Y43))
    * Column X: Rewrite to use expm1(x) = exp(x) - 1::
        = W ((abs(U)<1e-10 and abs(V)<1e-10) ? (V-U + (V-U)(V+U)/2) : (exp(-U)-exp(-V)))
        = W (exp(-U) - exp(-V))
        = W exp(-V) (exp(V-U) - 1) = W exp(-U) (1 - exp(U-V))
        = W exp(-V) expm1(V-U) = -W exp(-U) expm1(U-V)
        = (U > V) ? (W exp(-V) expm1(V-U)) : (-W exp(-U) expm1(U-V))

    Differences in the data tables:

    * AW1462 (W-186 => W-188 2n) t1/2 in hrs is not converting days to hours
    * AK1495 (Au-198 => Au-199 2n) target should be Au-197
    * AN1428 (Tm-169 => Tm-171 2n) t1/2 updated to Tm-171 rather than Tm-172
    * AN1420 (Er-162 => Ho-163 b) t1/2 updated to 4570 y from 10 y
    * AT1508 (Pb-208 => Pb-209 act) Thermal (b) x 1000 to convert from mbarns to barns
    """
    # TODO: is the table missing 1-H => 3-H ?
    # 0nly activations which produce radioactive daughter products are
    # included. Because 1-H => 2-H (act) is not in the table, is this why
    # there is no entry for 1-H => 3-H (2n).

    result: dict["ActivationResult", list[float]] = {}
    if not hasattr(isotope, 'neutron_activation'):
        return result

    for ai in isotope.neutron_activation:
        # Ignore fast neutron interactions if not using fast ratio
        if ai.fast and env.fast_ratio == 0:
            continue
        # Column D: elemental % mass content of sample
        #    mass fraction and abundance already included in mass calculation, so not needed
        # Column E: target nuclide and comment
        #    str(isotope), ai.comment
        # Column F: Nuclide Produced
        #    ai.daughter
        # Column G: Half-life
        #    ai.Thalf_str
        # Column H: initial effective cross-section [barn]
        #    env.epithermal_reduction_factor:$Q$1 = 1/env.Cd_ratio:$F$15
        initialXS = ai.thermalXS + env.epithermal_reduction_factor*ai.resonance
        # Column I: reaction
        #    ai.reaction
        # Column J: fast?
        #    ai.fast
        # Column K: effective reaction flux [n/cm^2/s]
        #    env.fluence:$F$13  env.fast_ratio:$F$17
        flux = env.fluence/env.fast_ratio if ai.fast else env.fluence
        # Column L: root part of activation calculation [uCi]
        #    mass:$F$19
        #    Decay correction portion done in column M
        #    The given mass is sample mass * sample fraction * isotope abundance
        #    The constant converts from Bq to uCi via avogadro's number with
        #        N_A[atoms] / 3.7e10[Bq/Ci] * 1e6 [uCi/Ci] ~ 1.627605611e19
        Bq_to_uCi = 1.6276e19
        #Bq_to_uCi = constants.avogadro_number / 3.7e4
        # 1/(cm^2 s) (Bq s) (barn cm^2/barn) g / (g/mol) ((1/mol) / (Bq/uCi))
        # TODO: use isotope.mass rather than isotope.isotope
        # Using isotope number rather than isotope mass to estimate the number
        # atoms introduces some error, overestimating Li activation by 0.25%
        # and underestimating Fe activation by 0.11%
        root = flux * (initialXS * 1e-24) * mass / isotope.isotope * Bq_to_uCi
        # Column M: 0.69/t1/2  [1/h] lambda of produced nuclide
        lam = LN2/ai.Thalf_hrs
        #print(ai.thermalXS, ai.resonance, env.epithermal_reduction_factor)
        #print(isotope, "D", mass, "F", ai.daughter, "G", ai.Thalf_str,
        #      "H", initialXS, "I", ai.reaction, "J", ai.fast, "K", flux,
        #      "L", root, "M", lam)

        # Column Y: activity at the end of irradiation [uCi]
        if ai.reaction == 'b':
            # Column N: 0.69/t1/2 [1/h] lambda of parent nuclide
            parent_lam = LN2 / ai.Thalf_parent
            # Column O: Activation if "b" mode production
            # 2022-05-18 PAK: addressed the following
            #    Note: problems resulting from precision limitation not addressed
            #    in "b" mode production
            # by rewriting equation to use expm1:
            #    activity = root*(1 - exp(-lam*exposure)/(1 - (lam/parent_lam))
            #                     + exp(-parent_lam*exposure)/((parent_lam/lam)-1))
            # Let x1=-lam*exposure x2=-parent_lam*exposure a=x1/x2=lam/parent_lam
            #    activity = root * (1 - exp(x1)/(1-a) + exp(x2)/(1/a - 1))
            #             = root * (1 - x2*exp(x1)/(x2-x1)) + x1*exp(x2)/(x2-x1))
            #             = root * ((x2-x1) - x2*exp(x1) + x1*exp(x2)) / (x2-x1)
            #             = root * (x2*(1-exp(x1)) + x1*(exp(x2)-1)) / (x2-x1)
            #             = root * (x1*expm1(x2) - x2*expm1(x1)) / (x2-x1)
            #             = root * (lam*expm1(x2) - parent_lam*expm1(x1)) / (parent_lam - lam)
            # Checked for each b-mode production that small halflife results are
            # unchanged to four digits and Eu[151] => Gd[152] no longer fails.
            activity = root/(parent_lam - lam) * (
                lam*expm1(-parent_lam*exposure) - parent_lam*expm1(-lam*exposure))
            #print("N", parent_lam, "O", activity)
        elif ai.reaction == '2n':
            # Column N: 0.69/t1/2 [1/h] lambda of parent nuclide
            parent_lam = LN2 / ai.Thalf_parent
            # Column P: effective cross-section 2n product and n, g burnup [barn]
            # Note: This cross-section always uses the total thermal flux
            effectiveXS = ai.thermalXS_parent + env.epithermal_reduction_factor*ai.resonance_parent
            # Column Q: 2n mode effective lambda of stable target [1/h]
            lam_2n = (flux*3600)*(initialXS*1e-24)
            # Column R: radioactive parent [1/h]
            parent_activity = (env.fluence*3600)*(effectiveXS*1e-24) + parent_lam
            # Column S: resulting product [1/h]
            product_2n = lam if ai.reaction == '2n' else 0
            # Column T: activity if 2n mode
            activity = root*lam*(parent_activity-parent_lam)*(
                (exp(-lam_2n*exposure)
                 / ((parent_activity-lam_2n)*(product_2n-lam_2n)))
                + (exp(-parent_activity*exposure)
                   / ((lam_2n-parent_activity)*(product_2n-parent_activity)))
                + (exp(-product_2n*exposure)
                   / ((lam_2n-product_2n)*(parent_activity-product_2n)))
                )
            #print("N", parent_lam, "P", effectiveXS, "Q", lam_2n,
            #      "R", parent_activity, "S", product_2n, "T", activity)
        else:
            # Provide the fix for the limitied precision (15 digits) in the
            # floating point calculation.  For neutron fluence rates above
            # 1e16 the precision in certain cells needs to be improved to
            # avoid erroneous results.  Also, burnup for single capture
            # reactions (excluding 'b') is included here.
            # See README file for details.

            # Column P: effective cross-section 2n product and n, g burnup [barn]
            # Note: This cross-section always uses the total thermal flux
            effectiveXS = ai.thermalXS_parent + env.epithermal_reduction_factor*ai.resonance_parent
            # Column U: nv1s1t [neutrons]
            U = (flux*3600)*(initialXS*1e-24)*exposure
            # Column V: nv2s2t+L2*t [neutrons]
            V = ((env.fluence*3600)*(effectiveXS*1e-24)+lam)*exposure
            # Column W: L/(L-nvs1+nvs2) [unitless]
            W = lam/(lam-(flux*3600)*(initialXS*1e-24)+(env.fluence*3600)*(effectiveXS*1e-24))
            # Column X: W*(exp(-U)-exp(V)) if U,V > 1e-10 else W*(V-U+(V-U)*(V+U)/2)
            # [PAK 2024-02-28] Rewrite the exponential difference using expm1()
            #Xp = exp(-V)*expm1(V-U) if U>V else -exp(-U)*expm1(U-V)
            X = W*exp(-V)*expm1(V-U) if U>V else -W*exp(-U)*expm1(U-V)
            # Column Y: O if "b" else T if "2n" else L*X
            activity = root*X
            #print(f"{ai.isotope}=>{ai.daughter} {U=} {V=} {W=} {Xp=} {X=} {activity=} {lam=} {flux*initialXS*3600*1e-24+env.fluence*effectiveXS*3600*1e-24}")
            #print(f"{ai.isotope}=>{ai.daughter} {U=} {V=} {W=} {X=} {activity=}")

            if activity < 0:
                msg = "activity %g less than zero for %g"%(activity, isotope)
                raise RuntimeError(msg)
            #print(ai.thermalXS_parent, ai.resonance_parent, exposure)
            #print("P", effectiveXS, "U", U, "V", V, "W", W, "X",
            #      precision_correction, "Y", activity)
            # columns: F32 H K L U V W X
            #data = env.fluence, initialXS, flux, root, U, V, W, precision_correction
            #print(" ".join("%.5e"%v for v in data))

        # TODO: chained activity (e.g., )
        result[ai] = [activity*exp(-lam*Ti) for Ti in rest_times]
        #print([(Ti, Ai) for Ti, Ai in zip(rest_times, result[ai])])

    return result

class ActivationResult:
    r"""
    *isotope* :

        Activation target for this result ({symbol}-{A})

    *abundance* : float | %

        IAEA 1987 isotopic abundance of activation target

    *symbol* :

        Element symbol for isotope

    *A* : int

        Number of protons plus neutrons in isotope

    *Z* : int

        Number of protons in isotope

    *reaction* :

        Activation type

        - "act" for (n,gamma)
        - "n,p" for (n,proton)
        - "n,a" for (n,alpha)
        - "2n" for activation of a daughter (e.g., 95Mo + n => 95Nb + n => 96Nb)
        - "n,2n" for neutron catalyzed release of a neutron
        - "b" for beta decay of a daughter (e.g., 98Mo + n => 99Mo => Tc-99)

    *comments* :

        Notes relating to simplifications or assumptions in the
        database. For most situations these do not affect your results.

    *daughter* :

        Daughter product from activation ({symbol}-{A}{isomer})

    *isomer* :

        Daughter product isotope annotation

        m, m1, m2: indicate metastable states.  Decay may be to the ground
        state or to another nuclide.

        \+: indicates radioactive daughter production already included in
        daughter listing several parent t1/2's required to acheive calculated
        daughter activity.  All activity assigned at end of irradiation.  In
        most cases the added activity to the daughter is small.

        \*: indicates radioactive daughter production NOT calculated,
        approximately secular equilibrium

        s: indicates radioactive daughter of this nuclide in secular equilibrium
        after several daughter t1/2's

        t: indicates transient equilibrium via beta decay.  Accumulation of that
        nuclide during irradiation is separately calculated.

    *Thalf_hrs* : float | hours

        Half-life of daughter in hours

    *Thalf_str* :

        Half-life of daughter as string, such as "29.0 y"

    *Thalf_parent* : float | hours

        Half-life of parent in chained "2n" or "b" reaction

    *fast* : bool

        Indicates whether the reaction is fast or thermal. If fast then the
        fluence is reduced by the specified fast/thermal ratio. When
        *fast_ratio* is zero in the environment this activation will not appear.

    *thermalXS*, *resonance*, *thermalXS_parent*, *resonance_parent* : float | barns

        Activation database values for computing the reaction cross section.

    *gT*, *percentIT* :

        Unused.

    Database notes:

    Most of the cross-section data is from IAEA 273. This is an excellent
    compilation. I highly recommend its purchase. The fast neutron cross-section
    data entered in the spreadsheet is weighted by an U-235 maxwellian
    distributed fission spectrum. Allthermal reactions producing nuclides with a
    half-life in excess of 1 second, and some less than 1 second, are included.
    Also included are radioactive daughters with halflives substantially longer
    than the parent produced by the neutron induced reaction, i.e., those
    daughters not in secular equilibrium. See the database in for more detailed
    notes relating to the database entries.

    [Note: 150 fission spectrum averaged fast neutron reactions added from
    IAEA156 on 2/1/94.  All those that are tabulated as measured have been
    entered.  Those estimated by calculation have not been entered.  In practice
    this means that the convenient, observable products are in this database.
    Fast reactions included are n,p; n,alpha; n,2n; and n,n'.

    Reaction = b  : This is the beta produced daughter of an activated parent.
    This is calculated only for the cases where the daughter is long lived
    relative to the parent. In the reverse case the daughter activity is
    reasonably self evident and the parent nuclide is tagged to indicate a
    radioactive daughter. The calculated activity of the beta produced
    daughter is through the end of irradiation. Contributions from the
    added decay of the parent after the end of irradiation are left for the
    user to determine, but are usually negligible for irradiations that are
    long relative to the parent halflife.

    During irradiation parent activity (A1) is:

        A1 = K [1 - exp(-L1*t)]

    where K is the parent saturation activity and L1 is the decay constant.

    For the beta produced daughter the activity (A2) is:

        A2 = K [1 - exp(-L1*t) * L2/(L2-L1) + exp(-L2*t) * L1/(L2-L1)]

    """
    isotope: Isotope
    abundance: float
    symbol: str
    A: int
    Z: int
    reaction: str
    comments: str
    daughter: str
    isomer: str
    Thalf_hrs: float
    Thalf_str: str
    Thalf_parent: float
    fast: bool
    thermalXS: float
    resonance: float
    thermalXS_parent: float
    resonance_parent: float
    percentIT: float|None
    gT: float|None

    def __init__(self, **kw):
        self.__dict__ = kw
    def __repr__(self):
        return f"ActivationResult({self.isotope},{self.reaction},{self.daughter})"
    def __str__(self):
        return f"{self.isotope}={self.reaction}=>{self.daughter}"

def init(table, reload=False):
    """
    Add neutron activation levels to each isotope.
    """
    # TODO: importlib.reload does not work for iso.neutron_activation attribute
    if 'neutron_activation' not in table.properties:
        table.properties.append('neutron_activation')
    elif not reload:
        return
    else:
        # Reloading activation table so clear the existing data.
        for el in table:
            for iso in el.isotopes:
                if hasattr(el[iso], 'neutron_activation'):
                    del el[iso].neutron_activation

    # We are keeping the table as a simple export of the activation data
    # from the ncnr health physics excel spreadsheet so that it is easier
    # to validate that the table contains the same data. Unfortunately some
    # of the cells involved formulas, which need to be reproduced when loading
    # in order to match full double precision values.
    activations = {}
    path = os.path.join(get_data_path('.'), 'activation.dat')
    with open(path, 'r') as fd:
        for row in fd:
            #print(row, end='')
            columns = row.split('\t')
            if columns[0].strip() in ('', 'xx'):
                continue
            columns = [c[1:-1] if c.startswith('"') else c for c in columns]
            #print columns
            for c in INT_COLUMNS:
                columns[c] = int(columns[c])
            for c in BOOL_COLUMNS:
                columns[c] = (columns[c] == 'y')
            for c in FLOAT_COLUMNS:
                s = columns[c].strip()
                columns[c] = float(s) if s else 0.
            # clean up comment column
            columns[-1] = columns[-1].replace('"', '').strip()
            kw = dict(zip(COLUMN_NAMES, columns))
            iso = (kw['Z'], kw['A'])
            act = activations.setdefault(iso, [])

            # TODO: use NuBase2020 for halflife
            # NuBase2020 uses different isomer labels.
            # Strip the (+, *, s, t) tags
            # Convert m1/m2 to m/n for
            #   In-114 In-116 Sb-124 Eu-152 Hf-179 Ir-192
            # Convert Eu-152m2 to Eu-152r
            # Convert m to n for
            #   Sc-46m Ge-73m Kr-83m Pd-107m Pd-109m Ag-110m Ta-182m Ir-194m Pb-204m
            # Convert m to p for
            #   Sb-122m Lu-177m

            # Recreate Thalf_hrs column using double precision.
            # Note: spreadsheet is not converting half-life to hours in cell AW1462 (186-W => 188-W)
            kw['Thalf_hrs'] = float(kw['_Thalf']) * UNITS_TO_HOURS[kw['_Thalf_unit']]
            #print(f"T1/2 {kw['Thalf_hrs']} +/- {kw['Thalf_hrs_unc']}")

            # Half-lives use My, Gy, Ty, Py
            value, units = float(kw['_Thalf']), kw['_Thalf_unit']
            if units == 'y':
                if value >= 1e12:
                    value, units = value/1e12, 'Ty'
                elif value >= 1e9:
                    value, units = value/1e9, 'Gy'
                elif value >= 200e3: # above 200,000 years use My
                    value, units = value/1e6, 'My'
                elif value >= 20e3: # between 20,000 and 200,000 use ky
                    value, units = value/1e3, 'ky'
            formatted = f"{value:g} {units}"
            #if formatted.replace(' ', '') != kw['Thalf_str'].replace(' ', ''):
            #    print(f"{kw['_index']}: old {kw['Thalf_str']} != {formatted} new")
            kw['Thalf_str'] = formatted

            # Override half-lives using NUBASE2020.
            # Note: to use the original values, import periodictable.activation and set
            # use_nubase2020_halflife = False before accessing any activation data.
            if use_nubase2020_halflife:
                daughter = kw['daughter']
                # NUBASE2020 lists 123Te as stable, so drop it from the reaction table
                if daughter == "Te-123":
                    continue
                Thalf_hrs = nubase2020_halflife[daughter][0] / 3600.0
                # Show change for a specific element or isotope
                #if daughter.startswith('Co'):
                #    print(f"{daughter}: {kw['Thalf_hrs']} => {Thalf_hrs}")
                kw['Thalf_hrs'] = Thalf_hrs
                kw['Thalf_str'] = nubase2020_halflife[daughter][2]

            # Recreate Thalf_parent by fetching from the new Thalf_hrs
            # e.g., =IF(OR(AR1408="2n",AR1408="b"),IF(AR1407="b",AW1406,AW1407),"")
            # This requires that the parent is directly before the 'b' or 'nb'
            # with its activation list already entered into the isotope.
            # Note: 150-Nd has 'act' followed by two consecutive 'b' entries.
            if kw['reaction'] in ('b', '2n'):
                parent = act[-2] if act[-1].reaction == 'b' else act[-1]
                kw['Thalf_parent'] = parent.Thalf_hrs
            else:
                #assert kw['Thalf_parent'] == 0
                kw['Thalf_parent'] = None

            # Strip columns whose names start with underscore
            kw = dict((k, v) for k, v in kw.items() if not k.startswith('_'))

            # Create an Activation record and add it to the isotope
            act.append(ActivationResult(**kw))

            # Check abundance values
            #if abs(iso.abundance - kw['abundance']) > 0.001*kw['abundance']:
            #    percent = 100*abs(iso.abundance - kw['abundance'])/kw['abundance']
            #    print "Abundance of", iso, "is", iso.abundance, \
            #        "but activation.dat has", kw['abundance'], "(%.1f%%)"%percent

    # Plug the activation products into the table
    for (Z, A), daughters in activations.items():
        table[Z][A].neutron_activation = tuple(daughters)

# Halflife from NUBASE2020 in seconds
# F.G. Kondev, M. Wang, W.J. Huang, S. Naimi, and G. Audi, Chin. Phys. C45, 030001 (2021)
#   isomer: (halflife, uncertainty, formatted)
use_nubase2020_halflife = True # Set this to False before accessing activation data
nubase2020_halflife = {
    "H-3": (3.8878e+08, 6.3e+05, "12.32 y"),
    "Li-8": (8.3870e-01, 3.0e-04, "838.7 ms"),
    "Be-10": (4.377e+13, 3.8e+11, "1.387 My"),
    "He-6": (8.0692e-01, 2.4e-04, "806.92 ms"),
    "B-12": (2.0200e-02, 2.0e-05, "20.20 ms"),
    "C-11": (1.2204e+03, 3.2e-01, "20.3402 m"),
    "C-14": (1.7987e+11, 9.5e+08, "5.70 ky"),
    "N-16": (7.13e+00, 2.0e-02, "7.13 s"),
    "O-15": (1.2227e+02, 4.3e-02, "122.266 s"),
    "N-17": (4.173e+00, 4.0e-03, "4.173 s"),
    "O-19": (2.6470e+01, 6.0e-03, "26.470 s"),
    "F-20": (1.1006e+01, 8.0e-03, "11.0062 s"),
    "F-18": (6.5840e+03, 4.8e-01, "109.734 m"),
    "Ne-23": (3.715e+01, 3.0e-02, "37.15 s"),
    "Na-24m+": (2.0180e-02, 1.0e-04, "20.18 ms"),
    "Na-24": (5.38416e+04, 5.4e+00, "14.9560 h"),
    "Na-22": (8.2108e+07, 1.9e+04, "2.6019 y"),
    "Na-25": (5.91e+01, 6.0e-01, "59.1 s"),
    "Mg-27": (5.661e+02, 1.6e+00, "9.435 m"),
    "Mg-28": (7.5294e+04, 3.2e+01, "20.915 h"),
    "Al-28": (1.347e+02, 3.0e-01, "2.245 m"),
    "Al-29": (3.936e+02, 3.6e+00, "6.56 m"),
    "Si-31": (9.430e+03, 1.2e+01, "157.16 m"),
    "Si-32": (4.95e+09, 2.2e+08, "157 y"),
    "P-32": (1.23284e+06, 6.0e+02, "14.269 d"),
    "P-33": (2.1902e+06, 9.5e+03, "25.35 d"),
    "S-35": (7.5488e+06, 3.5e+03, "87.37 d"),
    "P-34": (1.243e+01, 1.0e-01, "12.43 s"),
    "S-37": (3.030e+02, 1.2e+00, "5.05 m"),
    "Cl-36": (9.508e+12, 4.7e+10, "301.3 ky"),
    "Cl-38m+": (7.150e-01, 3.0e-03, "715 ms"),
    "Cl-38": (2.2338e+03, 8.4e-01, "37.230 m"),
    "Ar-37": (3.0250e+06, 1.6e+03, "35.011 d"),
    "Ar-39": (8.46e+09, 2.5e+08, "268 y"),
    "Ar-41": (6.5766e+03, 2.4e+00, "109.61 m"),
    "Ar-42": (1.038e+09, 3.5e+07, "32.9 y"),
    "K-40": (3.9383e+16, 9.5e+13, "1.248 Gy"),
    "K-42": (4.4478e+04, 2.5e+01, "12.355 h"),
    "Ca-41": (3.137e+12, 4.7e+10, "99.4 ky"),
    "K-43": (8.028e+04, 3.6e+02, "22.3 h"),
    "Ca-45": (1.40495e+07, 7.8e+03, "162.61 d"),
    "Ca-47*": (3.9191e+05, 2.6e+02, "4.536 d"),
    "Ca-49t": (5.231e+02, 3.6e-01, "8.718 m"),
    "Sc-49": (3.4308e+03, 7.8e+00, "57.18 m"),
    "Sc-46m+": (1.875e+01, 4.0e-02, "18.75 s"),
    "Sc-46": (7.2366e+06, 1.2e+03, "83.757 d"),
    "Sc-47": (2.89371e+05, 5.2e+01, "3.3492 d"),
    "Ti-45": (1.1088e+04, 3.0e+01, "184.8 m"),
    "Sc-48": (1.5721e+05, 3.2e+02, "43.67 h"),
    "Ti-51": (3.456e+02, 6.0e-01, "5.76 m"),
    "V-52": (2.246e+02, 3.0e-01, "3.743 m"),
    "Cr-51": (2.393410e+06, 9.5e+01, "27.7015 d"),
    "Cr-49": (2.5380e+03, 6.0e+00, "42.3 m"),
    "V-53": (9.26e+01, 8.4e-01, "1.543 m"),
    "Cr-55": (2.098e+02, 1.8e-01, "3.497 m"),
    "Mn-56": (9.2840e+03, 3.6e-01, "2.5789 h"),
    "Mn-54": (2.69638e+07, 2.8e+03, "312.081 d"),
    "Fe-55": (8.6977e+07, 1.3e+04, "2.7562 y"),
    "Fe-53": (5.106e+02, 1.2e+00, "8.51 m"),
    "Fe-59": (3.8448e+06, 1.0e+03, "44.500 d"),
    "Co-60m+": (6.280e+02, 3.6e-01, "10.467 m"),
    "Co-61": (5.936e+03, 1.8e+01, "1.649 h"),
    "Co-60": (1.66349e+08, 1.9e+04, "5.2714 y"),
    "Co-58": (6.1209e+06, 1.7e+03, "70.844 d"),
    "Ni-57": (1.2816e+05, 2.2e+02, "35.60 h"),
    "Ni-59": (2.56e+12, 1.6e+11, "81 ky"),
    "Co-58m+": (3.1871e+04, 8.3e+01, "8.853 h"),
    "Ni-63": (3.194e+09, 4.7e+07, "101.2 y"),
    "Ni-65": (9.0630e+03, 1.8e+00, "2.5175 h"),
    "Cu-62": (5.803e+02, 4.8e-01, "9.672 m"),
    "Cu-64": (4.57214e+04, 4.7e+00, "12.7004 h"),
    "Cu-66": (3.072e+02, 8.4e-01, "5.120 m"),
    "Cu-67": (2.2259e+05, 4.3e+02, "61.83 h"),
    "Ni-66s": (1.966e+05, 1.1e+03, "54.6 h"),
    "Zn-65": (2.10764e+07, 3.5e+03, "243.94 d"),
    "Zn-69ms": (4.9489e+04, 4.0e+01, "13.747 h"),
    "Zn-69": (3.384e+03, 5.4e+01, "56.4 m"),
    "Zn-71m": (1.4933e+04, 4.3e+01, "4.148 h"),
    "Zn-71": (1.440e+02, 3.0e+00, "2.40 m"),
    "Ga-70": (1.2684e+03, 3.0e+00, "21.14 m"),
    "Ga-72m+": (3.968e-02, 1.3e-04, "39.68 ms"),
    "Ga-72": (5.0490e+04, 3.6e+01, "14.025 h"),
    "Ga-73": (1.750e+04, 1.1e+02, "4.86 h"),
    "Ge-69": (1.4058e+05, 3.6e+02, "39.05 h"),
    "Ge-71m+": (2.041e-02, 1.8e-04, "20.41 ms"),
    "Ge-71": (9.876e+05, 2.6e+03, "11.43 d"),
    "Ge-73m": (4.99e-01, 1.1e-02, "499 ms"),
    "Ge-75m+": (4.77e+01, 5.0e-01, "47.7 s"),
    "Ge-75": (4.9668e+03, 2.4e+00, "82.78 m"),
    "Zn-71ms": (1.4933e+04, 4.3e+01, "4.148 h"),
    "Ge-77m+": (5.37e+01, 6.0e-01, "53.7 s"),
    "Ge-77": (4.0360e+04, 1.1e+01, "11.211 h"),
    "Ge-78s": (5.280e+03, 6.0e+01, "88.0 m"),
    "As-74": (1.5353e+06, 1.7e+03, "17.77 d"),
    "As-76": (9.446e+04, 3.3e+02, "1.0933 d"),
    "As-77": (1.3964e+05, 1.8e+02, "38.79 h"),
    "Se-75": (1.03490e+07, 2.6e+03, "119.78 d"),
    "Se-77m": (1.736e+01, 5.0e-02, "17.36 s"),
    "Se-79m+": (2.340e+02, 1.1e+00, "3.900 m"),
    "Se-79": (1.032e+13, 8.8e+11, "327 ky"),
    "Se-79m": (2.340e+02, 1.1e+00, "3.900 m"),
    "Se-81m*": (3.4368e+03, 1.2e+00, "57.28 m"),
    "Se-81": (1.1070e+03, 7.2e+00, "18.45 m"),
    "Se-83m+": (7.01e+01, 4.0e-01, "70.1 s"),
    "Se-83t": (1.3350e+03, 2.4e+00, "22.25 m"),
    "Br-83": (8.546e+03, 1.4e+01, "2.374 h"),
    "Br-80m*": (1.59138e+04, 2.9e+00, "4.4205 h"),
    "Br-80": (1.0608e+03, 1.2e+00, "17.68 m"),
    "Br-82m+": (3.678e+02, 3.0e+00, "6.13 m"),
    "Br-82": (1.27015e+05, 2.5e+01, "35.282 h"),
    "Kr-79m+": (5.00e+01, 3.0e+00, "50 s"),
    "Kr-79": (1.2614e+05, 3.6e+02, "35.04 h"),
    "Kr-81m+": (1.310e+01, 3.0e-02, "13.10 s"),
    "Kr-81": (7.23e+12, 3.5e+11, "229 ky"),
    "Kr-83m": (6.588e+03, 4.7e+01, "1.830 h"),
    "Kr-85m+": (1.6128e+04, 2.9e+01, "4.480 h"),
    "Kr-85": (3.3854e+08, 2.2e+05, "10.728 y"),
    "Kr-87": (4.578e+03, 3.0e+01, "76.3 m"),
    "Kr-88": (1.0170e+04, 6.8e+01, "2.825 h"),
    "Rb-84": (2.8356e+06, 6.0e+03, "32.82 d"),
    "Rb-86m+": (6.10e+01, 1.8e-01, "1.017 m"),
    "Rb-86": (1.61093e+06, 6.9e+02, "18.645 d"),
    "Rb-88": (1.0668e+03, 1.8e+00, "17.78 m"),
    "Rb-89*": (9.192e+02, 6.0e+00, "15.32 m"),
    "Sr-85m+": (4.0578e+03, 2.4e+00, "67.63 m"),
    "Sr-85": (5.60269e+06, 5.2e+02, "64.846 d"),
    "Sr-87m": (1.0098e+04, 3.2e+01, "2.805 h"),
    "Sr-89": (4.3686e+06, 2.2e+03, "50.563 d"),
    "Sr-90": (9.1231e+08, 9.5e+05, "28.91 y"),
    "Y-88": (9.2127e+06, 2.1e+03, "106.629 d"),
    "Y-89m": (1.5663e+01, 5.0e-03, "15.663 s"),
    "Y-90m+": (1.1614e+04, 4.0e+01, "3.226 h"),
    "Y-90": (2.3058e+05, 1.8e+02, "64.05 h"),
    "Y-91": (5.0553e+06, 5.2e+03, "58.51 d"),
    "Zr-89": (2.82096e+05, 8.3e+01, "78.360 h"),
    "Zr-93": (5.08e+13, 1.6e+12, "1.61 My"),
    "Sr-90s": (9.1231e+08, 9.5e+05, "28.91 y"),
    "Zr-95s": (5.53236e+06, 5.2e+02, "64.032 d"),
    "Zr-97s": (6.0296e+04, 2.9e+01, "16.749 h"),
    "Nb-92mt": (8.740e+05, 1.1e+03, "10.116 d"),
    "Nb-92": (1.095e+15, 7.6e+13, "34.7 My"),
    "Nb-94m+": (3.758e+02, 2.4e-01, "6.263 m"),
    "Nb-94": (6.44e+11, 1.3e+10, "20.4 ky"),
    "Nb-93m": (5.087e+08, 3.8e+06, "16.12 y"),
    "Mo-93m+": (2.466e+04, 2.5e+02, "6.85 h"),
    "Mo-93": (1.26e+11, 2.5e+10, "4.0 ky"),
    "Nb-95": (3.02322e+06, 5.2e+02, "34.991 d"),
    "Nb-96": (8.406e+04, 1.8e+02, "23.35 h"),
    "Mo-99st": (2.37355e+05, 1.8e+01, "65.932 h"),
    "Tc-99": (6.662e+12, 3.8e+10, "211.1 ky"),
    "Zr-95": (5.53236e+06, 5.2e+02, "64.032 d"),
    "Mo-101s": (8.766e+02, 1.8e+00, "14.61 m"),
    "Tc-99m+": (2.16238e+04, 7.2e-01, "6.0066 h"),
    "Ru-97t": (2.4512e+05, 1.2e+02, "2.8370 d"),
    "Tc-97": (1.329e+14, 5.0e+12, "4.21 My"),
    "Ru-103": (3.39077e+06, 6.9e+02, "39.245 d"),
    "Ru-105t": (1.5980e+04, 4.0e+01, "4.439 h"),
    "Rh-105": (1.27228e+05, 6.8e+01, "35.341 h"),
    "Rh-103m": (3.3668e+03, 5.4e-01, "56.114 m"),
    "Rh-104m": (2.604e+02, 1.8e+00, "4.34 m"),
    "Rh-104": (4.23e+01, 4.0e-01, "42.3 s"),
    "Pd-103": (1.4680e+06, 1.6e+03, "16.991 d"),
    "Pd-107m+": (2.13e+01, 5.0e-01, "21.3 s"),
    "Pd-107": (2.051e+14, 9.5e+12, "6.5 My"),
    "Pd-109m+": (2.822e+02, 5.4e-01, "4.703 m"),
    "Pd-109": (4.892e+04, 4.3e+02, "13.59 h"),
    "Pd-111m*": (2.0027e+04, 4.7e+01, "5.563 h"),
    "Pd-111t": (1.4136e+03, 5.4e+00, "23.56 m"),
    "Ag-111": (6.4221e+05, 8.6e+02, "7.433 d"),
    "Ag-106ms": (7.154e+05, 1.7e+03, "8.28 d"),
    "Ag-108m*": (1.385e+10, 2.8e+08, "439 y"),
    "Ag-108": (1.429e+02, 6.6e-01, "2.382 m"),
    "Ag-110m*": (2.15882e+07, 2.1e+03, "249.863 d"),
    "Ag-110": (2.46e+01, 1.1e-01, "24.56 s"),
    "Cd-107": (2.3400e+04, 7.2e+01, "6.50 h"),
    "Ag-106": (1.4376e+03, 2.4e+00, "23.96 m"),
    "Cd-109": (3.9856e+07, 4.3e+04, "461.3 d"),
    "Cd-111m": (2.9100e+03, 5.4e+00, "48.50 m"),
    "Cd-113m": (4.383e+08, 3.5e+06, "13.89 y"),
    "Cd-113": (2.537e+23, 1.6e+21, "8.04 Py"),
    "Cd-115ms": (3.850e+06, 2.1e+04, "44.56 d"),
    "Cd-115s": (1.9246e+05, 1.8e+02, "53.46 h"),
    "Cd-117ms": (1.2388e+04, 3.2e+01, "3.441 h"),
    "Cd-117s": (9.011e+03, 1.8e+01, "2.503 h"),
    "In-114m2+": (4.310e-02, 6.0e-04, "43.1 ms"),
    "In-114m1*": (4.27766e+06, 8.6e+02, "49.51 d"),
    "In-114": (7.190e+01, 1.0e-01, "71.9 s"),
    "In-116m2+": (2.18e+00, 4.0e-02, "2.18 s"),
    "In-116m1": (3.257e+03, 1.0e+01, "54.29 m"),
    "In-116": (1.410e+01, 3.0e-02, "14.10 s"),
    "In-115m": (1.6150e+04, 1.4e+01, "4.486 h"),
    "Sn-113m+": (1.284e+03, 2.4e+01, "21.4 m"),
    "Sn-113": (9.9429e+06, 3.5e+03, "115.08 d"),
    "Sn-117m": (1.2043e+06, 2.1e+03, "13.939 d"),
    "Sn-119m": (2.5324e+07, 6.0e+04, "293.1 d"),
    "Sn-121m": (1.385e+09, 1.6e+07, "43.9 y"),
    "Sn-121": (9.731e+04, 1.4e+02, "27.03 h"),
    "Sn-123m": (2.4036e+03, 6.0e-01, "40.06 m"),
    "Sn-123": (1.1163e+07, 3.5e+04, "129.2 d"),
    "Sn-125t": (8.324e+05, 1.3e+03, "9.634 d"),
    "Sn-125m": (5.86e+02, 1.5e+01, "9.77 m"),
    "Sb-125": (8.7021e+07, 3.5e+04, "2.7576 y"),
    "Sn-126": (7.26e+12, 4.4e+11, "230 ky"),
    "Sb-122m+": (2.515e+02, 1.8e-01, "4.191 m"),
    "Sb-122": (2.35336e+05, 1.7e+01, "2.7238 d"),
    "Sb-124m2*": (1.212e+03, 1.2e+01, "20.2 m"),
    "Sb-124m1+": (9.30e+01, 5.0e+00, "93 s"),
    "Sb-124": (5.2013e+06, 2.6e+03, "60.20 d"),
    "Te-121m*": (1.4230e+07, 4.3e+04, "164.7 d"),
    "Te-121": (1.6684e+06, 6.0e+03, "19.31 d"),
    "Te-123m*": (1.02989e+07, 8.6e+03, "119.2 d"),
    "Te-125m": (4.959e+06, 1.3e+04, "57.40 d"),
    "Te-127m*": (9.167e+06, 6.0e+04, "106.1 d"),
    "Te-127": (3.366e+04, 2.5e+02, "9.35 h"),
    "Te-129ms": (2.9030e+06, 8.6e+03, "33.6 d"),
    "Te-129t": (4.176e+03, 1.8e+01, "69.6 m"),
    "I-129": (5.093e+14, 3.8e+12, "16.14 My"),
    "Te-131m": (1.1693e+05, 4.0e+02, "32.48 h"),
    "Te-132": (2.768e+05, 1.1e+03, "3.204 d"),
    "Te-131t": (1.5000e+03, 6.0e+00, "25.0 m"),
    "I-131": (6.93351e+05, 5.2e+01, "8.0249 d"),
    "I-126": (1.1172e+06, 4.3e+03, "12.93 d"),
    "I-128": (1.4994e+03, 1.2e+00, "24.99 m"),
    "Xe-125m+": (5.69e+01, 9.0e-01, "56.9 s"),
    "Xe-125t": (6.073e+04, 2.9e+02, "16.87 h"),
    "I-125": (5.13147e+06, 6.9e+02, "59.392 d"),
    "Xe-127m+": (6.92e+01, 9.0e-01, "69.2 s"),
    "Xe-127": (3.13995e+06, 2.6e+02, "36.342 d"),
    "Xe-129m": (7.672e+05, 1.7e+03, "8.88 d"),
    "Xe-131m": (1.0323e+06, 1.0e+03, "11.948 d"),
    "Xe-133m+": (1.899e+05, 1.1e+03, "2.198 d"),
    "Xe-133": (4.53375e+05, 4.3e+01, "5.2474 d"),
    "Xe-135m+": (9.174e+02, 3.0e+00, "15.29 m"),
    "Xe-135": (3.2904e+04, 7.2e+01, "9.14 h"),
    "Xe-137t": (2.291e+02, 7.8e-01, "3.818 m"),
    "Cs-137s": (9.480e+08, 1.3e+06, "30.04 y"),
    "Cs-134m+": (1.04832e+04, 7.2e+00, "2.912 h"),
    "Cs-134": (6.5165e+07, 1.3e+04, "2.0650 y"),
    "Cs-135": (4.20e+13, 6.0e+12, "1.33 My"),
    "I-130": (4.4496e+04, 3.6e+01, "12.36 h"),
    "Ba-131m+": (8.556e+02, 5.4e+00, "14.26 m"),
    "Ba-131s": (9.9533e+05, 8.6e+02, "11.52 d"),
    "Ba-133m+": (1.4004e+05, 2.2e+02, "38.90 h"),
    "Ba-133": (3.32544e+08, 5.0e+04, "10.5379 y"),
    "Cs-132": (5.5987e+05, 5.2e+02, "6.480 d"),
    "Ba-135m": (1.01196e+05, 7.2e+01, "28.11 h"),
    "Ba-136m": (3.084e-01, 1.9e-03, "308.4 ms"),
    "Ba-137m": (1.5312e+02, 6.0e-02, "2.552 m"),
    "Cs-136": (1.1241e+06, 4.3e+03, "13.01 d"),
    "Cs-137": (9.480e+08, 1.3e+06, "30.04 y"),
    "Ba-139": (4.9758e+03, 5.4e+00, "82.93 m"),
    "Ba-140": (1.10189e+06, 1.8e+02, "12.7534 d"),
    "La-140": (1.45040e+05, 1.4e+01, "40.289 h"),
    "La-141*": (1.411e+04, 1.1e+02, "3.92 h"),
    "Ce-137m*": (1.238e+05, 1.1e+03, "34.4 h"),
    "Ce-137t": (3.24e+04, 1.1e+03, "9.0 h"),
    "La-137": (1.89e+12, 6.3e+11, "60 ky"),
    "Ce-139m+": (5.76e+01, 3.2e-01, "57.58 s"),
    "Ce-139": (1.18923e+07, 1.7e+03, "137.642 d"),
    "Ce-141": (2.80843e+06, 8.6e+02, "32.505 d"),
    "Ce-143t": (1.18940e+05, 2.2e+01, "33.039 h"),
    "Pr-143": (1.1724e+06, 1.7e+03, "13.57 d"),
    "Ce-144s": (2.46142e+07, 2.2e+03, "284.886 d"),
    "Pr-142m+": (8.76e+02, 3.0e+01, "14.6 m"),
    "Pr-142": (6.883e+04, 1.4e+02, "19.12 h"),
    "Nd-147t": (9.4867e+05, 8.6e+02, "10.98 d"),
    "Pm-147": (8.27864e+07, 6.3e+03, "2.6234 y"),
    "Nd-149t": (6.2208e+03, 3.6e+00, "1.728 h"),
    "Pm-149": (1.9109e+05, 1.8e+02, "53.08 h"),
    "Nd-151t": (7.464e+02, 4.2e+00, "12.44 m"),
    "Pm-151": (1.0224e+05, 1.4e+02, "28.40 h"),
    "Sm-151": (2.985e+09, 1.9e+07, "94.6 y"),
    "Sm-145t": (2.938e+07, 2.6e+05, "340 d"),
    "Pm-145": (5.59e+08, 1.3e+07, "17.7 y"),
    "Sm-153": (1.666246e+05, 8.3e+00, "46.2846 h"),
    "Sm-155": (1.3308e+03, 3.6e+00, "22.18 m"),
    "Eu-152m2+": (5.748e+03, 2.4e+01, "95.8 m"),
    "Eu-152m1": (3.35218e+04, 4.7e+00, "9.3116 h"),
    "Gd-152": (3.41e+21, 2.5e+20, "108 Ty"),
    "Eu-152t": (4.2655e+08, 1.9e+05, "13.517 y"),
    "Eu-154": (2.71137e+08, 9.5e+04, "8.592 y"),
    "Eu-155": (1.4964e+08, 2.5e+05, "4.742 y"),
    "Gd-153": (2.0788e+07, 6.0e+04, "240.6 d"),
    "Gd-159": (6.6524e+04, 1.4e+01, "18.479 h"),
    "Gd-161t": (2.188e+02, 1.8e-01, "3.646 m"),
    "Tb-161": (6.0031e+05, 4.3e+02, "6.948 d"),
    "Tb-160": (6.247e+06, 1.7e+04, "72.3 d"),
    "Dy-157t": (2.930e+04, 1.4e+02, "8.14 h"),
    "Tb-157": (2.24e+09, 2.2e+08, "71 y"),
    "Dy-159": (1.2476e+07, 1.7e+04, "144.4 d"),
    "Dy-165m": (7.54e+01, 3.6e-01, "1.257 m"),
    "Dy-165": (8.395e+03, 1.4e+01, "2.332 h"),
    "Dy-166": (2.9376e+05, 3.6e+02, "81.6 h"),
    "Ho-166m": (3.574e+10, 1.2e+08, "1.1326 ky"),
    "Ho-166": (9.6523e+04, 2.5e+01, "26.812 h"),
    "Er-163t": (4.500e+03, 2.4e+01, "75.0 m"),
    "Ho-163": (1.4422e+11, 7.9e+08, "4.570 ky"),
    "Er-165": (3.730e+04, 1.4e+02, "10.36 h"),
    "Er-167m": (2.269e+00, 6.0e-03, "2.269 s"),
    "Er-169": (8.115e+05, 1.6e+03, "9.392 d"),
    "Er-171t": (2.70576e+04, 7.2e+00, "7.516 h"),
    "Tm-171": (6.059e+07, 3.2e+05, "1.92 y"),
    "Er-172*": (1.775e+05, 1.8e+03, "49.3 h"),
    "Tm-170": (1.1111e+07, 2.6e+04, "128.6 d"),
    "Yb-169m+": (4.60e+01, 2.0e+00, "46 s"),
    "Yb-169": (2.76601e+06, 4.3e+02, "32.014 d"),
    "Yb-175m+": (6.820e-02, 3.0e-04, "68.2 ms"),
    "Yb-175": (3.61584e+05, 8.6e+01, "4.185 d"),
    "Yb-177m+": (6.41e+00, 2.0e-02, "6.41 s"),
    "Yb-177t": (6.880e+03, 1.1e+01, "1.911 h"),
    "Lu-177": (5.74068e+05, 7.8e+01, "6.6443 d"),
    "Lu-176m": (1.3190e+04, 6.8e+01, "3.664 h"),
    "Lu-176": (1.1679e+18, 5.4e+15, "37.01 Gy"),
    "Lu-177m*": (1.3859e+07, 2.6e+04, "160.4 d"),
    "Hf-175": (6.104e+06, 1.6e+04, "70.65 d"),
    "Hf-178m": (4.0e+00, 2.0e-01, "4.0 s"),
    "Hf-179m2": (2.160e+06, 1.5e+04, "25.00 d"),
    "Hf-179m1": (1.867e+01, 4.0e-02, "18.67 s"),
    "Hf-180m": (1.9908e+04, 7.2e+01, "5.53 h"),
    "Hf-181": (3.6625e+06, 5.2e+03, "42.39 d"),
    "Hf-182s": (2.809e+14, 2.8e+12, "8.90 My"),
    "Ta-182m+": (9.504e+02, 6.0e+00, "15.84 m"),
    "Ta-182": (9.914e+06, 1.0e+04, "114.74 d"),
    "Ta-183": (4.406e+05, 8.6e+03, "5.1 d"),
    "Lu-178": (1.704e+03, 1.2e+01, "28.4 m"),
    "W-181": (1.04506e+07, 1.6e+03, "120.956 d"),
    "W-185m+": (9.58e+01, 2.4e-01, "1.597 m"),
    "W-185": (6.489e+06, 2.6e+04, "75.1 d"),
    "W-187": (8.5712e+04, 9.0e+01, "23.809 h"),
    "W-188": (6.0281e+06, 4.3e+03, "69.77 d"),
    "Re-184ms": (1.460e+07, 6.9e+05, "169 d"),
    "Re-184": (3.059e+06, 6.0e+04, "35.4 d"),
    "Re-186": (3.21278e+05, 4.3e+01, "3.7185 d"),
    "Re-188m+": (1.1154e+03, 2.4e+00, "18.59 m"),
    "Re-188": (6.1218e+04, 1.1e+01, "17.005 h"),
    "Os-185": (8.0309e+06, 7.8e+03, "92.95 d"),
    "Os-189m": (2.092e+04, 3.6e+02, "5.81 h"),
    "Os-190m": (5.916e+02, 1.8e+00, "9.86 m"),
    "Os-191m+": (4.716e+04, 1.8e+02, "13.10 h"),
    "Os-191": (1.2951e+06, 1.7e+03, "14.99 d"),
    "Os-193": (1.07388e+05, 6.5e+01, "29.830 h"),
    "Os-194s": (1.893e+08, 6.3e+06, "6.0 y"),
    "Ir-192m2*": (7.61e+09, 2.8e+08, "241 y"),
    "Ir-192m1+": (8.70e+01, 3.0e+00, "1.45 m"),
    "Ir-192": (6.3780e+06, 1.2e+03, "73.820 d"),
    "Ir-193m": (9.098e+05, 3.5e+03, "10.53 d"),
    "Ir-194m": (1.477e+07, 9.5e+05, "171 d"),
    "Ir-194": (6.966e+04, 2.5e+02, "19.35 h"),
    "Pt-191": (2.445e+05, 1.7e+03, "2.83 d"),
    "Pt-193m": (3.741e+05, 2.6e+03, "4.33 d"),
    "Pt-193": (1.58e+09, 1.9e+08, "50 y"),
    "Pt-195m": (3.4646e+05, 4.3e+02, "4.010 d"),
    "Pt-197m+": (5.725e+03, 1.1e+01, "95.41 m"),
    "Pt-197": (7.16094e+04, 6.8e+00, "19.8915 h"),
    "Pt-199m+": (1.35e+01, 1.6e-01, "13.48 s"),
    "Pt-199t": (1.848e+03, 1.3e+01, "30.80 m"),
    "Au-199": (2.7121e+05, 6.0e+02, "3.139 d"),
    "Pt-200s": (4.54e+04, 1.1e+03, "12.6 h"),
    "Au-196": (5.3266e+05, 9.5e+02, "6.165 d"),
    "Au-197m": (7.73e+00, 6.0e-02, "7.73 s"),
    "Au-198": (2.32817e+05, 1.2e+01, "2.69464 d"),
    "Hg-197m*": (8.575e+04, 1.4e+02, "23.82 h"),
    "Hg-197": (2.3375e+05, 2.5e+02, "64.93 h"),
    "Hg-199m": (2.5602e+03, 5.4e+00, "42.67 m"),
    "Hg-203": (4.02710e+06, 8.6e+02, "46.610 d"),
    "Hg-205": (3.084e+02, 5.4e+00, "5.14 m"),
    "Tl-202": (1.0636e+06, 6.9e+03, "12.31 d"),
    "Tl-204": (1.1938e+08, 3.8e+05, "3.783 y"),
    "Tl-206": (2.521e+02, 6.6e-01, "4.202 m"),
    "Pb-203": (1.86926e+05, 5.4e+01, "51.924 h"),
    "Pb-204m": (4.0158e+03, 6.0e+00, "66.93 m"),
    "Pb-205": (5.36e+14, 2.8e+13, "17.0 My"),
    "Pb-209": (1.1646e+04, 1.8e+01, "3.235 h"),
    "Bi-210ms": (9.59e+13, 1.9e+12, "3.04 My"),
    "Bi-210t": (4.3304e+05, 4.3e+02, "5.012 d"),
    "Po-210": (1.195569e+07, 1.7e+02, "138.376 d"),
    "Bi-211": (1.284e+02, 1.2e+00, "2.14 m"),
}

def demo():  # pragma: nocover
    import sys
    import argparse
    import periodictable as pt

    parser = argparse.ArgumentParser(description='Process some data with mass, flux, exposure, and decay options.')
    parser.add_argument('-m', '--mass', type=float, default=1, help='Specify the mass value (default: 1)')
    parser.add_argument('-f', '--flux', type=float, default=1e8, help='Specify the flux value (default: 1e8)')
    parser.add_argument('-e', '--exposure', type=float, default=10, help='Specify the exposure value (default: 10)')
    parser.add_argument('-d', '--decay', type=float, default=5e-4, help='Specify the decay value (default: 5e-4)')
    parser.add_argument('--cd-ratio', type=float, default=70, help='Specify the Cd ratio value (default: 70)')
    parser.add_argument('--fast-ratio', type=float, default=50, help='Specify the fast ratio value (default: 50)')

    parser.add_argument('formula', nargs='?', type=str, default=None, help='Specify the formula as a positional argument')
    args = parser.parse_args()

    formula = args.formula
    if formula is None:
        # Make sure all elements compute
        #formula = "".join(str(el) for el in pt.elements)[1:]
        formula = build_formula([(1/el.mass, el) for el in pt.elements][1:])
        # Use an enormous mass to force significant activation of rare isotopes
        mass, fluence = 1e15, 1e8
    env = ActivationEnvironment(fluence=args.flux, Cd_ratio=args.cd_ratio, fast_ratio=args.fast_ratio, location="BT-2")
    sample = Sample(formula, mass=args.mass)
    abundance = IAEA1987_isotopic_abundance
    #abundance=NIST2001_isotopic_abundance,
    sample.calculate_activation(
        env, exposure=args.exposure, rest_times=(0, 1, 24, 360),
        abundance=abundance,
        )
    decay_time = sample.decay_time(args.decay)
    print(f"{args.mass} g {formula} for {args.exposure} hours at {args.flux} n/cm^2/s")
    print(f"Time to decay to {args.decay} uCi is {decay_time} hours.")
    sample.calculate_activation(
        env, exposure=args.exposure, rest_times=(0, 1, 24, 360, decay_time),
        abundance=abundance,
        )
    sample.show_table(cutoff=0.0)

    ## Print a table of flux vs. activity so we can debug the
    ## precision_correction value in the activity() function.
    ## Note that you also need to uncomment the print statement
    ## at the end of activity() that shows the column values.
    #import numpy as np
    #sample = Sample('Co', mass=10)
    #for fluence in np.logspace(3, 20, 20-3+1):
    #    env = ActivationEnvironment(fluence=fluence)
    #    sample.calculate_activation(
    #        env, exposure=exposure, rest_times=[0],
    #        abundance=IAEA1987_isotopic_abundance)

if __name__ == "__main__":
    demo()  # pragma: nocover
