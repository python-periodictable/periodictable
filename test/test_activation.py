import numpy as np

import periodictable as pt
from periodictable.activation import Sample, ActivationEnvironment
from periodictable.activation import IAEA1987_isotopic_abundance, table_abundance

def test():
    # This is not a very complete test of the activation calculator.
    # Mostly just a smoke test to see that things run and produce the
    # same answers as before.  The target values have been checked
    # against the NCNR internal activation calculator spreadsheet.  The
    # values herein differ slightly from the spreadsheet since we are using
    # different tables for materials and different precision on our
    # constants. There has also been some small corrections to the
    # formulas over time.

    def _get_Au_activity(fluence=1e5):
        sample = Sample('Au', mass=1)
        env = ActivationEnvironment(fluence=fluence)
        sample.calculate_activation(env, rest_times=[0])
        for product, activity in sample.activity.items():
            if str(product.daughter) ==  'Au-198':
                return activity[0]
        else:
            raise RuntimeError("missing activity from Au-198")

    # Activity scales linearly with fluence and mass
    act1 = _get_Au_activity(fluence=1e5)
    act2 = _get_Au_activity(fluence=1e8)
    #print(f"Au: {act1} x 1000 = {act2}")
    assert (act2 - 1000*act1) < 1e-8

    # Smoke test: does every element run to completion?
    formula = "".join(str(el) for el in pt.elements)
    # Use an enormous mass to force significant activation of rare isotopes
    mass, fluence, exposure = 1e15, 1e8, 10
    env = ActivationEnvironment(fluence=fluence, Cd_ratio=70, fast_ratio=50, location="BT-2")
    sample = Sample(formula, mass=mass)
    sample.calculate_activation(
        env, exposure=exposure, rest_times=(0, 1, 24),
        abundance=IAEA1987_isotopic_abundance,
        #abundance=table_abundance,
        )
    #sample.show_table(cutoff=0)

    sample = Sample('Co', mass=10)
    env = ActivationEnvironment(fluence=1e8)
    sample.calculate_activation(
        env, exposure=1, rest_times=[0, 1, 24],
        abundance=IAEA1987_isotopic_abundance,
        #abundance=table_abundance,
        )
    #sample.show_table(cutoff=0)
    # ACT_2N_X.xls for 10g Co at Io=1e8 for t = [0, 1, 24]
    # There are duplicate entries for Co-61 because there are different
    # activation paths for 59Co + n -> 61Co.
    # Results are limited to 3 digits because formulas use ln(2) = 0.693. There
    # is also precision loss on the calculation of differences in exponentials,
    # with exp(a) - exp(b) converted to exp(b)(exp(a-b) - 1) = exp(b)expm1(a-b)
    # in activation.py.
    # results = {
    #     "C0-60m+":  [5.08800989419543E+03,	9.69933141298983E+01,	2.69898447909379E-38],
    #     "Co-61":    [7.30937687202485E-09,	4.80260282859366E-09,	3.06331725449002E-13],
    #     "Co-60":    [1.55042700053951E-01,	1.55040373560730E-01,	1.54986873850802E-01],
    #     "Co-61":    [1.36469792582999E-09,	8.96670432174800E-10,	5.71936948464344E-14],
    # }
    # Results from R2.1.0-pre, with halflife updated to NUBASE2020
    results = [
        ("Co-60m+", [5.08870610591228E+03,	9.57100873523175E+01,	1.95440969864725E-38]),
        ("Co-61", [7.29613649387052E-09,	4.79211894136011E-09,	3.03056972241711E-13]),
        ("Co-60", [1.54971558843742E-01,	1.54969234199432E-01,	1.54915777003567E-01]),
        # Stable two-capture result, checked against a 100-digit calculation.
        ("Co-61", [1.3655818015418685E-09, 8.969172140136206E-10, 5.672167543347562E-14]),
    ]
    check_table(sample, results)



    # 129-I has a long half-life (16 My) so any combination of exposure
    # fluence and mass that causes significant activation will yield a product
    # that is radioactive for a very long time.
    sample = Sample('Te', mass=1e13)
    env = ActivationEnvironment(fluence=1e8)
    sample.calculate_activation(env, rest_times=[1, 10, 100])
    #sample.show_table(cutoff=0)
    target = 1e-5
    t_decay = sample.decay_time(target)
    #print(f"{t_decay=}")
    sample.calculate_activation(env, rest_times=[t_decay])
    total = sum(v[-1] for v in sample.activity.values())
    assert abs(total - target) < 1e-10, f"total activity {total} != {target}"

    # Al and Si daughters have short half-lives
    sample = Sample('AlSi', mass=1e3)
    env = ActivationEnvironment(fluence=1e8)
    sample.calculate_activation(env, rest_times=[100, 200])
    #sample.show_table(cutoff=0)
    target = 1e-5
    t_decay = sample.decay_time(target)
    #print(f"{t_decay=}")
    sample.calculate_activation(env, rest_times=[t_decay])
    total = sum(v[-1] for v in sample.activity.values())
    assert abs(total - target) < 1e-10, f"total activity {total} != {target}"

    #print()
    if 0: # TODO: doesn't work for high fluence
        # Test high fluence
        sample = Sample('Al2O3', mass=1e3)
        sample = Sample('Al', mass=1e3)
        flow, scale = 1e16, 1e2
        rest_times = [0, 10]
        env = ActivationEnvironment(fluence=flow*scale)
        sample.calculate_activation(env, rest_times=rest_times)
        sample.show_table(cutoff=0)
        ahigh = list(sample.activity.values())
        env = ActivationEnvironment(fluence=flow)
        sample.calculate_activation(env, rest_times=rest_times)
        sample.show_table(cutoff=0)
        alow = np.asarray(list(sample.activity.values()))
        for alow_k, ahigh_k in zip(alow, ahigh):
            #print(f"{alow_k=} {ahigh_k=}")
            assert np.allclose(alow_k*scale, ahigh_k), f"scaling fails at {alow_k=} {ahigh_k=}"

def test_bateman():
    sample = Sample('Nd', mass=10)
    env = ActivationEnvironment(fluence=1e8, Cd_ratio=70, fast_ratio=50)
    sample.calculate_activation(
        env, exposure=1, rest_times=[0, 1, 24],
        abundance=IAEA1987_isotopic_abundance,
        #abundance=table_abundance,
        )
    # Spreadsheet results
    #    ("Nd-147t", [7.1846042E-02	7.1658714E-02	6.7482432E-02],
    #    ("Pm-147",  [1.0837323E-06	1.0836997E-06	1.0829483E-06].
    #    ("Nd-149t", [5.4279655E+00	3.6934626E+00	5.2693381E-04],
    #    ("Pm-149",  [3.7518439E-02	3.7031973E-02	2.7429311E-02],
    #    ("Nd-151t", [8.2839488E+00	2.5906397E-01	6.3431983E-36],
    #    ("Pm-151",  [1.4885091E-01	1.4526270E-01	8.2873173E-02],
    #    ("Sm-151",  [5.4151602E-06	5.4151554E-06	5.4150460E-06],
    # Results match v1.3.5 (initial conversion)
    #    ("Nd-147t", [7.18460819622254E-02,	7.16587537625008E-02,	6.74824694565637E-02]),
    #    ("Pm-147",  [1.08373233644659E-06,	1.08369965659330E-06,	1.08294829189398E-06]),
    #    ("Nd-149t", [5.42796550374046E+00,	3.69346258031675E+00,	5.26933810918582E-04]),
    #    ("Pm-149",  [3.75184393106657E-02,	3.70319731847190E-02,	2.74293110532566E-02]),
    #    ("Nd-151t", [8.28394877194627E+00,	2.59063974896533E-01,	6.34319831491801E-36]),
    #    ("Pm-151",  [1.48850911418719E-01,	1.45262696324095E-01,	8.28731733878711E-02]),
    #    ("Sm-151",  [5.41516018833055E-06,	5.41515542843156E-06,	5.41504595190962E-06]),
    # Results from v2.1.0 (switch to NuBase 2020 for half-life)
    #    ("Nd-147t", [7.23754310236518E-02,	7.21853087887262E-02,	6.79477209334824E-02]),
    #    ("Pm-147",  [1.09122954601348E-06,	1.09119665490308E-06,	1.09044043292477E-06]), # decreasing
    #    ("Nd-149t", [5.61218053312601E+00,	3.75772037482054E+00,	3.69964003980992E-04]),
    #    ("Pm-149",  [3.89111773898684E-02,	3.84063628404178E-02,	2.84423696213705E-02]), # decreasing
    #    ("Nd-151t", [8.24828184506010E+00,	2.91373687385306E-01,	1.17610906034480E-34]),
    #    ("Pm-151",  [1.47015198605385E-01,	1.43470491863699E-01,	8.18409551527891E-02]), # decreasing
    #    ("Sm-151",  [5.08521250925607E-06,	5.08520825824191E-06,	5.08511048589707E-06]), # decreasing
    # Results from PR #130 (switch to 3-stage Bateman equation for Nd[150] decay)
    results = [
        ("Nd-147t", [7.23754310236518E-02,	7.21853087887262E-02,	6.79477209334824E-02]),
        ("Pm-147",  [1.09122954601348E-06,	3.26982190261643E-06,	5.18101337117162E-05]), # increasing
        ("Nd-149t", [5.61218053312601E+00,	3.75772037482054E+00,	3.69964003980992E-04]),
        ("Pm-149",  [3.89111773898684E-02,	9.83584332412775E-02,	1.66469845209970E-01]), # increasing
        ("Nd-151t", [8.24828184506010E+00,	2.91373687385306E-01,	1.17610906034480E-34]),
        ("Pm-151",  [1.47015198605385E-01,	2.00524291103734E-01,	1.15608922581670E-01]), # up then down
        ("Sm-151",  [4.97572645610762E-08,	2.06630403513325E-07,	3.18791122026443E-06]), # lower initial, increasing
    ]
    check_table(sample, results) # use results=[] to print table


def check_table(sample, results):
    #print(list(sample.activity.keys()))
    #print(" ".join(k for k in dir(list(sample.activity.keys())[0]) if k[0] != '_'))
    #print(list(sample.activity.keys())[0].__dict__)
    for k, (product, activity) in enumerate(sample.activity.items()):
        #print(product)
        #print(dir(product))
        if results: # Check
            assert product.daughter == results[k][0], f"Expected {product.daughter} but got {results[k][0]}"
            # Test that results haven't changed since last update
            assert np.allclose(results[k][1], activity, atol=0, rtol=1e-12)
        else: # Show results table
            activity_str = ",\t".join(f"{Ia:.14E}" for Ia in activity)
            print(f'        ("{product.daughter}", [{activity_str}]),')

test()
test_bateman()
