"""Exotic chiral character, oriented pair production and finite ATLAS count."""
import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
SW2 = 0.23153  # same illustrative tree scheme as round 1775
MZ, GZ = 91.1876, 2.4952  # nominal GeV parameters, fixed-width propagator
QF, TF = -7, -0.5  # down component of the k=2 exotic lepton doublet
BEAMS = {"up": (2/3, 0.5), "down": (-1/3, -0.5)}
TEMPLATE = {500: 2300, 1100: 49, 2000: 0.31}  # ATLAS 2303.13613, Table 2


def exotic_character():
    z = cmath.exp(2j*math.pi/11)
    weights = [(6, 13), (3, -16), (3, -10), (2, -39), (1, 36), (1, 42)]
    char = sum(dim*z**power for dim, power in weights)
    assert abs(char.imag) > 1 and sum(dim for dim, _ in weights) == 16
    return {"real": char.real, "imaginary": char.imag, "mirror_imaginary": -char.imag}


def amplitudes(s, charge, t3):
    weak = np.array([t3, 0.]) - charge*SW2
    chi = s/(s-MZ**2+1j*MZ*GZ)/(SW2*(1-SW2))
    v0 = charge*QF + weak*chi*(-QF*SW2)
    shift = weak*chi*TF/2
    return v0, shift


def collider_lower(mass, charge, t3):
    delta = MZ*math.hypot(MZ, GZ)/(4*mass**2-MZ**2)
    t = np.array([t3, 0.])
    g = t-charge*SW2
    vinf = abs(QF*(charge-t)/(1-SW2))
    zpart = abs(g*QF/(1-SW2))
    delta_axial = (1+delta)*abs(g*TF)/(2*SW2*(1-SW2))
    numerator = sum(np.maximum(0, vinf - delta*zpart - 2*delta_axial)**2)
    denominator = sum((vinf+delta*zpart)**2)
    return numerator/denominator


def scheme_envelope(mass, charge, t3):
    lo, hi = .22, .24
    delta = 92*math.hypot(92,3)/(4*mass**2-92**2)
    t = np.array([t3,0.])
    gmax = np.maximum(abs(t-charge*lo), abs(t-charge*hi))
    zbound = delta*gmax*abs(QF)/(1-hi)
    vinf_lower = abs(QF*(charge-t))/(1-lo)
    vinf_upper = abs(QF*(charge-t))/(1-hi)
    twice_axial = (1+delta)*gmax*abs(TF)/(lo*(1-hi))
    return sum(np.maximum(0,vinf_lower-zbound-twice_axial)**2)/sum((vinf_upper+zbound)**2)


def angular(s, m, charge, t3, z, orientation):
    beta = math.sqrt(max(0, 1-4*m*m/s))
    v0, d = amplitudes(s, charge, t3)
    v, a = v0+d, orientation*d
    kv = 2-beta**2 + beta**2*z**2
    ka = beta**2*(1+z**2)
    odd = 4*beta*z*np.real(v*a.conjugate())
    return sum(kv*abs(v)**2+ka*abs(a)**2) + odd[0]-odd[1], sum(kv*abs(v0)**2)


def afb(s, m, orientation):
    # Three-point Gaussian quadrature on each half integrates the tree polynomial exactly.
    nodes, weights = np.polynomial.legendre.leggauss(3)
    forward = sum(w*angular(s,m,-1,-0.5,(1+z)/2,orientation)[0] for z,w in zip(nodes,weights))
    backward = sum(w*angular(s,m,-1,-0.5,(z-1)/2,orientation)[0] for z,w in zip(nodes,weights))
    return (forward-backward)/(forward+backward)


def calculate():
    rows = []
    for m, n in TEMPLATE.items():
        bounds = {label:collider_lower(m,*beam) for label,beam in BEAMS.items()}
        rho = min(*bounds.values(), 1.)
        wide = min(1.,*(scheme_envelope(m,*beam) for beam in BEAMS.values()))
        assert 0 < wide <= rho
        sampled = 1.
        for b in BEAMS.values():
            for x in [1,1.02,1.4,2,10,100]:
                for angle in np.linspace(-1,1,17):
                    for chirality in [-1,1]:
                        num, den = angular((2*m*x)**2,m,*b,angle,chirality)
                        assert num/den + 1e-12 >= rho
                        sampled = min(sampled, num/den)
        rows.append(dict(mass_GeV=m, unpolarized_DY_guarantees=bounds, overall_guarantee=rho,
                         scan_min_including_photons=sampled, reference_yield=n,
                         guaranteed_tree_yield=n*rho, sufficient_retention=-math.log(.05)/(n*rho),
                         scheme_envelope=wide, envelope_yield=n*wide,
                         envelope_sufficient_retention=-math.log(.05)/(n*wide)))
    example = afb(1200**2,500,1)
    assert abs(example + afb(1200**2,500,-1)) < 1e-12
    assert rows[0]["guaranteed_tree_yield"]>1000 and rows[1]["guaranteed_tree_yield"]>30
    return {"character_at_11th_root":exotic_character(),"Qf":QF,
            "e_minus_beam_1200_GeV_m500_AFB":{"family":example,"mirror":-example},
            "observed_absQ_above_2":0,"poisson_upper_95_count":-math.log(.05),"points":rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    result = json.dumps(calculate(), indent=2, allow_nan=False)+"\n"
    path = HERE/"independent_results.json"
    if args.check:
        assert json.loads(result) == json.loads(path.read_bytes())
    else:
        path.write_text(result,encoding="utf-8")
    print("1779 exotic Weyl character and finite MCP sensitivity",hashlib.sha256(result.encode()).hexdigest()[:12])
