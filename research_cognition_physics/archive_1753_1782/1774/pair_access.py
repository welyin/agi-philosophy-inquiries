"""Accessible two-line minor and a two-copy, forward-delay flavour-rank witness."""
import argparse
import itertools
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def accessible_minor(v, pair):
    a = v[:, pair]
    return sum(abs(np.linalg.det(a[[i,j]]))**2
               for i,j in itertools.combinations(range(len(v)),2))


def high_pair(weights, omega, time):
    return sum(weights[i]*weights[j]*2*math.sin((omega[i]-omega[j])*time/2)**2
               for i,j in itertools.combinations(range(len(weights)),2))


def antisymmetric_flag(a, b):
    return np.linalg.norm((np.kron(a,b)-np.kron(b,a))/2)**2


def mixed_flag(state, phases):
    other=phases[:,None]*state*phases.conj()[None,:]
    return (np.trace(state).real**2-np.trace(state@other).real)/2


def exposure(success, miss=.05):
    if not 0 < success < 1:
        raise ValueError("A nontrivial single-trial probability is required.")
    return math.ceil(math.log(miss)/math.log1p(-success))


def run():
    fit = json.loads((HERE.parent/"1766/independent_results.json").read_bytes())["refits"]["c1h"]
    v3 = runpy.run_path(str(HERE.parent/"1773/independent_check.py"))["unitary_from_tree"](fit)
    eta = .001  # Candidate mixing probability, not measured.
    v = np.eye(4,dtype=complex)
    v[:3,:3] = v3
    r = np.eye(4)
    r[2:,2:]=[[math.sqrt(1-eta), math.sqrt(eta)],[-math.sqrt(eta), math.sqrt(1-eta)]]
    v = v@r
    visible=v[:3]
    blind = accessible_minor(visible,[2,3])
    assert blind < 1.e-28 and np.min(abs(visible[:,2:])) > 0
    for pair in itertools.combinations(range(4), 2):
        prob = accessible_minor(visible,list(pair))
        gram = visible[:,pair].conj().T@visible[:,pair]
        np.testing.assert_allclose(prob,np.linalg.det(gram).real,atol=1.e-14)
        np.testing.assert_allclose(prob,1-np.sum(abs(v[3,list(pair)])**2),atol=1.e-14)
    obs={}
    for i,key in enumerate(("u","c","t")):
        modes=visible[i,2:]
        probs=abs(modes)**2
        peak=high_pair(probs,[0.,1.],math.pi)
        np.testing.assert_allclose(peak,2*eta*(1-eta)*abs(v3[i,2])**4,rtol=1.e-12)
        np.testing.assert_allclose(peak,antisymmetric_flag(modes,modes*[1,-1]),rtol=1.e-12)
        obs[key]=dict(weak_source_b_probability=abs(v3[i,2])**2,
                      peak_flag_per_prepared_pair=peak,
                      ideal_trials_for_95pct_at_least_one=exposure(peak))
    # Arbitrary many high modes; projected flag is supported entirely above target.
    high=np.array([.31,.12j,-.05,.16j])
    omega=np.array([1.4,2.1,5.3,6.2])
    w=abs(high)**2
    t=1.73
    observed=antisymmetric_flag(high,high*np.exp(-1j*omega*t))
    np.testing.assert_allclose(observed,high_pair(w,omega,t),rtol=1.e-13)
    # Independent midpoint numerical integration of a random forward wait.
    duration=2.9
    average=sum(w[j]*w[k]*(1-np.sinc((omega[j]-omega[k])*duration/math.pi))
                for j,k in itertools.combinations(range(4),2))
    n=5000
    mid=(np.arange(n)+.5)*duration/n
    overlaps=np.exp(-1j*np.outer(mid,omega))@w
    checked=np.mean((w.sum()**2-abs(overlaps)**2)/2)
    assert abs(average-checked)<1.e-9
    mixed=np.diag(w)
    flag_mixed=(np.trace(mixed)**2-np.trace(mixed@mixed))/2
    assert flag_mixed>0 and high_pair(w,omega,0)==0
    # Tight delayed-minus-undelayed bound with a single true high flavour.
    total,leak=.8,.012
    pure=np.sqrt([total-leak,leak])
    pure_state=np.outer(pure,pure)
    phases=np.array([1.,-1.])
    leakage_cap=2*leak*(total-leak)
    np.testing.assert_allclose(mixed_flag(pure_state,phases)-mixed_flag(pure_state,np.ones(2)),leakage_cap)
    diag_state=np.diag(pure**2)
    assert mixed_flag(diag_state,phases)==mixed_flag(diag_state,np.ones(2))
    # One true high flavour with leakage: determinant in each output pair containing
    # a low flavour is bounded by the raw probability of at least one such entry.
    rng=np.random.default_rng(1774)
    for _ in range(96):
        a,b=(rng.normal(size=(2,3))+1j*rng.normal(size=(2,3)))
        a*=.87/np.linalg.norm(a)
        b*=.77/np.linalg.norm(b)
        la,lb=np.linalg.norm(a[1:])**2,np.linalg.norm(b[1:])**2
        ceiling=la*np.vdot(b,b).real+lb*np.vdot(a,a).real-la*lb
        assert antisymmetric_flag(a,b) <= ceiling+1.e-13
        mixture=.65*np.outer(a,a.conj())+.35*np.outer(b,b.conj())
        second=.15*np.outer(a,a.conj())+.85*np.outer(b,b.conj())
        h_a,h_b=mixture[0,0].real,second[0,0].real
        total_a,total_b=np.trace(mixture).real,np.trace(second).real
        l_a,l_b=total_a-h_a,total_b-h_b
        raw=(total_a*total_b-np.trace(mixture@second).real)/2
        lower_overlap=max(0.,math.sqrt(h_a*h_b)-math.sqrt(l_a*l_b))**2
        assert raw<=(total_a*total_b-lower_overlap)/2+1.e-13
        diag=mixture.diagonal().real
        high=np.argmax(diag)
        low=diag.sum()-diag[high]
        phase=np.exp(1j*rng.uniform(-math.pi,math.pi,3))
        diff=mixed_flag(mixture,phase)-mixed_flag(mixture,np.ones(3))
        if low<diag[high]:
            assert -1.e-13<=diff<=2*low*diag[high]+1.e-13
    return dict(hypothetical_eta=eta, single_pass_high_pair_all_visible_row_minor_zero=bool(blind<1.e-28),
                common_raw_state="two successful weak flavour preparations; high-band losses included",
                delayed_sources=obs, many_high_modes=dict(raw_pure_delayed=observed,
                    random_wait_mean=average, mixed_incoherent_no_delay=flag_mixed),
                one_high_mode_example=dict(total_pass_weight=total,low_mode_leakage=leak,
                    saturable_false_modulation=leakage_cap))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    output=json.loads(json.dumps(run()),parse_float=lambda s:float(f"{float(s):.11g}"))
    path=HERE/"pair_access_results.json"
    if args.check:
        assert output==json.loads(path.read_bytes())
    else:
        with path.open("x",encoding="utf8") as f:
            json.dump(output,f,indent=2,allow_nan=False)
            f.write("\n")
    print("1774: restricted two-column minors and finite-delay high-band witness verified")
