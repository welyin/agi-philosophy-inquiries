"""Single-event heat and local clock/propagation readouts on a finite Dirac cycle."""

import argparse
import hashlib
import json
import math
import runpy
from fractions import Fraction
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
PARAM = json.loads((HERE / "inputs.json").read_text(encoding="utf-8"))
N = PARAM["cycle_sites"]
ELL = float(Fraction(PARAM["spacing"]))
EPS = float(Fraction(PARAM["clock_gain"]))
NU = float(Fraction(PARAM["bond_gain"]))
T = float(Fraction(PARAM["time"]))
ORDER = PARAM["taylor_order"]
MASSES = [float(Fraction(s)) for s in PARAM["clock_masses"]]
MOM = 2*np.pi/N
WAVE = np.exp(1j * MOM * np.arange(N))
POSITION = (1+np.sin(MOM*np.arange(N)))/2
Z = np.array([1, -1])[:, None, None]
SPIN_Z = np.array([1, -1])[None, :, None]
SQRT2 = np.sqrt(2)


def generator(psi, mass, d, clock_gain=EPS, bond_gain=NU):
    """Coordinates are relation qubit, Dirac spin and location. One tagged bond."""
    kinetic = -0.5j/ELL * (np.roll(psi, -1, axis=2)-np.roll(psi, 1, axis=2))
    out = -psi[::-1] + d*Z*psi + mass*SPIN_Z*psi + kinetic[:, ::-1]
    local_mass = clock_gain * mass * (Z*SPIN_Z*psi)[:, :, 0]
    scale = -0.5j * bond_gain/ELL
    out[:, :, 0] += local_mass + scale * (Z*psi[:, ::-1])[:, :, 1]
    out[:, :, 1] -= scale * (Z*psi[:, ::-1])[:, :, 0]
    return out


def series(mass, d, initial, clock_gain, bond_gain):
    terms = [initial.copy()]
    for n in range(1, ORDER + 1):
        terms.append(-1j/n * generator(terms[-1], mass, d, clock_gain, bond_gain))
    return np.array(terms)


def evaluate(terms, t):
    out = np.zeros_like(terms[0])
    for term in reversed(terms):
        out = term + t*out
    return out


def prepare(clock_gain=EPS, bond_gain=NU):
    # G starts in |+X>, clock spin in |+beta>, massless spin in |+alpha_x>.
    clock = np.zeros((2, 2, N), complex)
    clock[:, 0] = 1/np.sqrt(2*N)
    position = np.ones((2, 2, 1))*(1+WAVE)/(2*np.sqrt(2*N))
    return {d: [series(mass, d, clock, clock_gain, bond_gain) for mass in MASSES] +
               [series(0, d, position, clock_gain, bond_gain)] for d in (-1, 0, 1)}


def event(branch, t):
    low, high, packet = [evaluate(b, t) for b in branch]
    for psi in (low, high, packet):
        np.testing.assert_allclose(np.vdot(psi, psi), 1, atol=3e-14, rtol=0)
    clock = (1-np.vdot(low, high).imag)/2
    space = np.sum(np.abs(packet)**2 * POSITION)
    return np.array([clock, space])


def aggregate(spec, weights, t):
    return sum(w*event(spec[d], t) for d, w in weights.items() if w)


def weights(r, p, mix):
    if r:
        return {-1: (1-p)*(1-mix), 0:p, 1:(1-p)*mix}
    return {-1:p*mix, 0:1-p, 1:p*(1-mix)}


def compact(o):
    if isinstance(o, (list, np.ndarray)):
        return [compact(v) for v in o]
    if isinstance(o, dict):
        return {str(k):compact(v) for k,v in o.items()}
    if isinstance(o, (float,np.floating)):
        return float(f"{o:.13g}")
    return o


def run():
    old = runpy.run_path(str(HERE.parent/'1818'/'independent_check.py'))
    src_t = HERE.parents[1]/"archive_1234_1280"/"1278"/"check.py"
    previous = runpy.run_path(str(src_t))
    assert callable(previous["series"]) and previous["T"] == Fraction(1,8)
    old_source = old['run']()
    assert old_source['conditional_maps_checked_on_four_matrix_units'] == 24
    spec = prepare()
    branch = {d:event(s,T) for d,s in spec.items()}
    offset = branch[1]-branch[-1]
    p = old_source["initial_thermal_population"]
    witness = {}
    for r in range(2):
        reference = aggregate(spec, weights(r,p,0), T)
        rows = []
        for x in (0,.5,1):
            state, ledger = old["reset"](math.asin(math.sqrt(x)), r, p)
            # These are physical B,A occupation probabilities in the previous reset.
            ba = old["reduced"](state,[1,2], (2,2,2))
            measured = {-1:ba[1,1].real,0:(ba[0,0]+ba[3,3]).real,1:ba[2,2].real}
            np.testing.assert_allclose(list(measured.values()),list(weights(r,p,x).values()),atol=1e-14)
            readouts = aggregate(spec, measured, T)
            change = readouts-reference
            heat = ledger['heat_into_bath_gap_units']
            np.testing.assert_allclose(change,heat*offset,atol=7e-16,rtol=0)
            rows.append(dict(mix=x, bath_heat=heat, clock_and_position=readouts,
                             clock_and_position_change=change))
        witness[r]=rows
    # Equal-clock comparison: bracket one physical time, interpolate baseline by bisection.
    # In the tight bracket use a first-derivative lower bound from the H-norm curvature.
    baseline=weights(1,p,0)
    target=np.array(witness[1][-1]["clock_and_position"])
    lo,hi=T-.00024,T+.00024
    left,right=aggregate(spec,baseline,lo),aggregate(spec,baseline,hi)
    assert left[0] < target[0] < right[0]
    # Each full Hamiltonian norm <= Hcap. For any probability |p''|<=4Hcap².
    hcap=2+max(MASSES)+1/ELL+EPS*max(MASSES)+NU/(2*ELL)
    local_base = [(d,w) for d,w in baseline.items() if w]
    local_derivative = 0.
    for d,w in local_base:
        low,high = [evaluate(b,T) for b in spec[d][:2]]
        low_d,high_d = [-1j*generator(b,m,d) for b,m in zip((low,high),MASSES)]
        local_derivative -= .5*w*(np.vdot(low_d,high)+np.vdot(low,high_d)).imag
    monotonic_lower=local_derivative - 4*hcap**2*(hi-T)
    assert monotonic_lower>0.04
    for _ in range(44):
        mid=(lo+hi)/2
        if aggregate(spec,baseline,mid)[0] < target[0]:
            lo=mid
        else:
            hi=mid
    tmatch=(lo+hi)/2
    matched=aggregate(spec,baseline,tmatch)
    free=event(prepare(0,0)[0],T)
    free_expected=[.5+.5*math.sin((MASSES[1]-MASSES[0])*T),
                   .5+.25*math.sin(math.sin(MOM)*T/ELL)]
    np.testing.assert_allclose(free,free_expected,atol=4e-15,rtol=0)
    tail=(hcap*.126)**(ORDER+1)/math.factorial(ORDER+1)
    assert tail<1e-29
    return compact(dict(conditional_readouts_by_D=branch,
                        slopes_per_unit_bath_heat=offset,
                        position_to_clock_delta_ratio=offset[1]/offset[0],
                        conditional_histories=witness,
                        equal_clock=dict(baseline_time=tmatch,target_time=T,
                            central_derivative=local_derivative,minimum_derivative_in_bracket=monotonic_lower,
                            delta_position=target[1]-matched[1],delta_clock=target[0]-matched[0]),
                        neutral_free_reference=free,
                        hamiltonian_norm_bound=hcap, analytic_state_tail_at_time_0126=tail))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    generated=run()
    output=HERE/'results.json'
    if args.check:
        expected=json.loads(output.read_text(encoding='utf-8'))
        def verify(a,b):
            if isinstance(a,dict):
                assert a.keys()==b.keys()
                for k in a: verify(a[k],b[k])
            elif isinstance(a,list):
                assert len(a)==len(b)
                for u,v in zip(a,b): verify(u,v)
            else:
                np.testing.assert_allclose(a,b,rtol=2e-8,atol=8e-15)
        verify(generated,expected)
    else:
        output.write_text(json.dumps(generated,indent=2)+"\n",encoding='utf-8')
    print('1820',hashlib.sha256(output.read_bytes()).hexdigest()[:12])
