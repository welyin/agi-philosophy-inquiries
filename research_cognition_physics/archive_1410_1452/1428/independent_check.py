"""Two eigenphases: bounds on amplitude moduli, reciprocal profiles, and open counterexample."""
import argparse
import json
from math import pi, sin, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE / name).read_text(encoding="utf-8-sig"))


def pairs(unitary):
    return abs(unitary[..., 1, 0])**2, abs(unitary[..., 0, 1])**2


def pmns(x, y, z, delta):
    a, b, c = np.sqrt([x, y, z])
    A, B, C = np.sqrt(1-np.array([x,y,z]))
    d = np.exp(1j*delta)
    return np.array([[A*B,a*B,b/d],[-a*C-A*b*c*d,A*C-a*b*c*d,c*B],
                     [a*c-A*b*C*d,-A*c-a*b*C*d,C*B]])


def path(vectors, m2, density, distances, data, energies, reverse=False):
    base = (vectors * m2) @ vectors.conj().T
    h_vac = base / energies[:,None,None]
    hc = data["hbarc_GeV_cm"]
    pot = (2e18*sqrt(2)*data["GF"]*data["NA"]*hc**3*data["electron_fraction"])
    result = np.broadcast_to(np.eye(3), h_vac.shape).copy()
    seq = list(zip(density, distances))
    for rho, length in (reversed(seq) if reverse else seq):
        h = h_vac.copy()
        h[:,0,0] += rho*pot
        vals, basis = np.linalg.eigh(h)
        phases = np.exp(-1j*length*(1e-13/(2*hc))*vals)
        step = (basis*phases[:,None,:]) @ basis.conj().transpose(0,2,1)
        result = step @ result
    return result


def profile_test(s, window=None):
    data = read("../1423/inputs.json")
    pos = np.array(data["position_km"])
    distances = np.diff(np.r_[pos[0],(pos[1:]+pos[:-1])/2,pos[-1]])
    ev = np.linspace(*(window or data["energy_GeV"]),1001)
    y = (1-sqrt(1-data["theta13_amplitude"]))/2
    central_phase, min_phase, max_phase = data["alpha_rad"]
    rows = {}
    for name, density in data["density_g_cm3"].items():
        for ordering,(x,z,solar,atm) in data["mass_and_angles"].items():
            m2 = np.array([0,solar,atm])
            def signals(es, angle):
                vectors = pmns(x,y,z,angle)
                f = path(vectors,m2,density,distances,data,es)
                r = path(vectors,m2,density,distances,data,es,True)
                p, q = pairs(f)
                v, w = pairs(r)
                return (p-q+v-w)/2, 1-abs(f[:,0,0])**2, np.stack([p,q,v,w],axis=-1)
            diff, loss, _ = signals(ev, pi/2)
            derivative = (1e-13/(2*data["hbarc_GeV_cm"]))*sum(distances)*np.ptp(m2)/(2*ev[0]**2)
            halfstep = (ev[1]-ev[0])/2
            assert np.all(diff>0)
            d_lo = (np.min(diff)-4*derivative*halfstep)*min(sin(min_phase),sin(max_phase))
            loss_hi = np.max(loss)+2*derivative*halfstep
            phases = (central_phase,-central_phase,pi-central_phase,pi+central_phase)
            for angle in phases:
                sampled = ev[::100]
                dd, ll, pair = signals(sampled, angle)
                assert np.allclose(dd, np.sin(angle)*diff[::100], atol=1e-13)
                assert np.allclose(ll, loss[::100], atol=1e-13)
                assert np.all(pair <= ll[:,None]+1e-14)
            small_rate = 4*s*sqrt(loss_hi)
            assert d_lo > small_rate
            # If every competing channel lies within epsilon of the projection, then
            # its absolute average difference >= d_lo-2 eps and all four <= loss_hi+eps.
            lo, hi = 0., d_lo/2
            for _ in range(70):
                mid = (lo+hi)/2
                if d_lo-2*mid > 4*s*sqrt(loss_hi+mid):
                    lo=mid
                else:
                    hi=mid
            rows[name+"_"+ordering] = {
                "absolute_D_mean_lower":float(d_lo),"all_four_probabilities_upper":float(loss_hi),
                "two_eigenphase_budget_using_rates":float(small_rate),
                "sufficient_per_probability_error_strictly_below":float(lo),
                "central_four_probabilities_alpha": signals(ev[500:501],central_phase)[2][0].tolist()}
    return rows


def checks():
    z = read("../1422/inputs.json")
    u = z["Z_reanalysis"][0]-z["matching_widths"]*z["Z_reanalysis"][1]
    minimum=sqrt(u-2)
    s=(1-minimum)/(1+minimum)
    # Rank-one projector onto (|e>+i |e_perp>)/sqrt(2) attains the signed bound.
    e = np.array([1.,0.])
    mu = np.array([s,sqrt(1-s*s)])
    v = np.array([1.,1j])/sqrt(2)
    for phase in (.2,pi/2,2.3):
        mat = np.eye(2)+np.expm1(1j*phase)*np.outer(v,v.conj())
        a,b = mu@mat@e,e@mat@mu
        assert np.isclose(abs(abs(a)**2-abs(b)**2),2*s*sqrt(1-s*s)*abs(sin(phase)))
        half = np.exp(-.5j*phase)
        assert abs(half*a+(half*b).conj()-2*s*np.cos(phase/2)) < 1e-14
    had = np.array([[1,1],[1,-1]])/sqrt(2)
    second = had.T@np.diag([1,-1])
    f,r=had@second,second@had
    pf,qf=pairs(f)
    pr,qr=pairs(r)
    assert np.isclose(abs(pf-qr),1) and np.isclose(pf-qf+pr-qr,0)
    # Ordered basis e0, e1, mu0, mu1. Fix mu0 and transfer probability e0 -> mu1.
    v = np.array([1.,0.,0.,1.])/sqrt(2)
    proj = np.outer(v,v)
    phase=2*np.arcsin(sqrt(.024))
    mat = np.eye(4)+np.expm1(1j*phase)*proj
    inclusive = [float(sum(abs(mat[2:,0])**2)),float(sum(abs(mat[:2,2])**2))]
    rank_one = [float(abs(mat[2,0])**2),float(abs(mat[0,2])**2)]
    assert np.allclose(inclusive,[.024,0]) and np.allclose(rank_one,[0,0])
    return {"port_overlap_bound":s,"same_path_D_absolute_budget":2*s*sqrt(1-s*s),
            "single_crossed_difference_two_state_example":float(abs(pf-qr)),
            "inclusive_bath_probabilities":inclusive,"reciprocal_rank_one_probabilities":rank_one,
            "trace_distance_of_unbiased_mixture_after_bath":(inclusive[0]-inclusive[1])/2,
            "finite_window":profile_test(s)}


def compare(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(compare(a[k],b[k]) for k in a)
    return np.allclose(a,b,rtol=1e-9,atol=1e-12)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    modes=p.add_mutually_exclusive_group()
    modes.add_argument("--check",action="store_true")
    modes.add_argument("--write",action="store_true")
    args=p.parse_args()
    result=checks()
    if args.write:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as h:
            json.dump(result,h,indent=2)
            h.write("\n")
    else:
        assert compare(result,read("independent_results.json"))
        author = read("results.json")["windows"]
        for lo,hi in read("inputs.json")["finite_windows_GeV"]:
            rows = profile_test(result["port_overlap_bound"], [lo,hi])
            for label, row in rows.items():
                other = author[f"{lo}-{hi}_{label}"]
                assert row["absolute_D_mean_lower"] > other["absolute_four_rate_lower"]
                assert row["sufficient_per_probability_error_strictly_below"] > other["rate_dependent_per_probability_error_budget"]
        print("1428 independent bounded four-channel response and two-phase examples verified")
