"""TEXONO finite-t certificate and LEP 207 GeV angular shape with free luminosity."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
GF = 1.1663787e-5
ME = .00051099895  # GeV
G = .63
TMAX = .010
TABLE = np.array([
    (1.440,.196,1.339), (1.426,.163,1.517), (1.889,.177,1.745),
    (2.156,.198,2.240), (3.215,.233,3.194), (4.434,.357,4.380),
    (6.393,.463,5.729), (6.951,.481,7.972), (11.221,.615,11.019),
    (15.933,.739,16.053), (25.676,.923,25.254), (42.075,1.188,42.456),
    (77.611,1.569,79.639), (173.825,3.002,178.042), (573.637,6.024,576.688),
])
EDGES = np.array([-.9,-.72,-.54,-.36,-.18,0,.09,.18,.27,.36,.45,.54,.63,.72,.81,.9])
ENERGY = 207.
E2 = 4 * math.pi / 128
S2 = .2312
Z2 = E2 / (S2 * (1-S2))
ZTERM = 1 / complex(1-(91.1876/ENERGY)**2, 91.1876*2.4952/ENERGY**2)
HL, HR = -.5+S2, S2
MEASURED, ERROR, BHWIDE = TABLE.T


def reactor():
    upper = 1.08 + 2 * math.hypot(.21, .16)
    transfer, retain = .8, .90
    a, b = .76, .25
    norm = 4 * math.sqrt(2) * GF
    target = math.sqrt(upper / transfer)
    delta_limit = max(b*(1+target), a*(target-1)) / retain
    min_mass = math.sqrt(G**2 / (norm*delta_limit) - 2*ME*TMAX)
    fractions = {}
    for m in (0., .001, 1., 32., 64.):
        q2m2 = m*m + 2*ME*TMAX
        d = retain * G**2 / (norm*q2m2)
        assert d > b
        lower = min((d/b-1)**2, (1+d/a)**2)
        fractions[str(m)] = dict(delta_retained=d, min_bin_ratio=lower,
                                 min_fitted_ratio=transfer*lower,
                                 maximum_response_weight=upper/lower,
                                 coupling_ceiling=math.sqrt(delta_limit*norm*q2m2))
    return dict(limit_two_widths=upper, delta_max=delta_limit,
                mass_floor_GeV=min_mass, cases=fractions)


def hull(a, b):
    return (np.minimum(a,b), np.maximum(a,b))


def scaled(pair, multiplier):
    return hull(pair[0]*multiplier, pair[1]*multiplier)


def add(a, b):
    return (a[0]+b[0], a[1]+b[1])


def square(pair):
    lo, hi = pair
    return (np.maximum(np.maximum(lo,-hi),0.)**2, np.maximum(lo*lo,hi*hi))


def angular(mlo, mhi, g, n, width_max, midpoint=False):
    """Positive interval integrals in angle and mass, stripping one common constant."""
    mesh = EDGES[:-1,None] + np.diff(EDGES)[:,None]*np.arange(n+1)/n
    lo, hi = mesh[:,:-1], mesh[:,1:]
    if midpoint:
        lo = hi = (lo+hi)/2
    # All propagators are dimensionless (multiplied by s).
    t0, t1 = (1-hi)/2, (1-lo)/2
    gamma_t = hull(-E2/t0,-E2/t1)
    z_t = hull(-1/(t0+(91.1876/ENERGY)**2),-1/(t1+(91.1876/ENERGY)**2))
    ml, mh = (mlo/ENERGY)**2, (mhi/ENERGY)**2
    family_t = hull(-1/(t0+ml),-1/(t1+mh))
    omega = width_max*mh/(1-mh)
    family_s = (1/((1-ml)*(1+omega**2)),1/(1-mh))
    # |Im s/(s-M²+iMΓ)|≤omega/(1-mh).
    im_budget = g*g*omega/(2*(1-mh))
    result = np.zeros((2,15,n))
    for i, j in [(HL,HL),(HR,HR),(HL,HR)]:
        product = i*j*Z2
        strength = g*g/2 if i == j else -g*g/2
        s_real = add((E2+product*ZTERM.real,)*2, scaled(family_s,strength))
        s_im2 = square((product*ZTERM.imag-im_budget, product*ZTERM.imag+im_budget))
        t_real = add(add(gamma_t,scaled(z_t,product)),scaled(family_t,strength))
        if i == j:
            total = add(square(add(s_real,t_real)),s_im2)
            result += np.array([total[0]*(1+lo)**2,total[1]*(1+hi)**2])
        else:
            sch = add(square(s_real),s_im2)
            tch = square(t_real)
            result += np.array([2*sch[0]*(1-hi)**2+8*tch[0],
                                2*sch[1]*(1-lo)**2+8*tch[1]])
    return result.mean(axis=2)


def norm_interval(ratio_low, ratio_high):
    lower = np.max((MEASURED-2*ERROR)/(BHWIDE*ratio_high))
    upper = np.min((MEASURED+2*ERROR)/(BHWIDE*ratio_low))
    return float(lower), float(upper)


def required_shape(lower, upper):
    return max(0.,(lower-upper)/(lower+upper))


def collider():
    n = 256
    base = angular(0, 0, 0, n, 0)
    bound = 1.
    for i in range(640):
        lo, hi = i/10, (i+1)/10
        low, high = angular(lo,hi,G,n,.1)
        intervals = norm_interval(low/base[1]-1e-9,high/base[0]+1e-9)
        bound = min(bound,required_shape(*intervals))
    # Full interval integration rather than point quadrature for the positive witness.
    low, high = angular(25,25,G,n,.015)
    witness_low = np.max((MEASURED-2*ERROR)/(BHWIDE*(low/base[1]-1e-9)*1.05))
    witness_high = np.min((MEASURED+2*ERROR)/(BHWIDE*(high/base[0]+1e-9)*.95))
    # If intervals prove too wide, refine the fixed-mass integration.
    if witness_low > witness_high:
        fine_base = angular(0,0,0,3072,0)
        low, high = angular(25,25,G,3072,.015)
        witness_low = np.max((MEASURED-2*ERROR)/(BHWIDE*(low/fine_base[1]-1e-9)*1.05))
        witness_high = np.min((MEASURED+2*ERROR)/(BHWIDE*(high/fine_base[0]+1e-9)*.95))
    assert witness_low < witness_high
    examples = {}
    ref = angular(0,0,0,1536,0,True)[0]
    for m in (0.,25.,64.):
        ratio = angular(m,m,G,1536,0,True)[0]/ref
        l,u = norm_interval(ratio,ratio)
        examples[str(m)] = dict(common_norm=(l+u)/2, needed_shape=required_shape(l,u),
                                ratio_backward=ratio[0], ratio_forward=ratio[-1])
    return dict(angular_bins=len(TABLE), angle_subdivisions=n, mass_subdivisions=640,
                minimum_required_shape_certified=bound,
                global_witness_25_GeV_common_norm=[witness_low,witness_high], examples=examples)


def stable(obj):
    if isinstance(obj, dict):
        return {k:stable(v) for k,v in obj.items()}
    if isinstance(obj, list):
        return [stable(v) for v in obj]
    return float(f"{obj:.9g}")


if __name__ == "__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check", action="store_true")
    args = cli.parse_args()
    out = stable(dict(TEXONO=reactor(), LEP207=collider()))
    assert out["LEP207"]["minimum_required_shape_certified"] > .03
    path = HERE / "independent_results.json"
    if args.check:
        assert json.loads(path.read_text(encoding="utf-8-sig")) == out
    else:
        with path.open("x", encoding="utf-8") as f:
            json.dump(out, f, indent=2, ensure_ascii=False)
            f.write("\n")
    print("TEXONO and LEP angular-shape independent checks passed.")
