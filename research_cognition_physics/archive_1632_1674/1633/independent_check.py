"""Independent pole-readout, representation and mass-derivative tests for round 1633."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# Tables 2.13, 3.6, 5.10 of the final LEP/SLD report, hep-ex/0509008.
POLE = {"R_b": (.21629, .00066), "R_e": (20.804, .050),
        "A_e": (.1516, .0021), "A_FB_b": (.0992, .0016), "A_b_direct": (.923, .020)}


def pole_readouts():
    out = {}
    def width_factor(rb, re, ae, ab):
        return 3*rb*re*(1-ab)/(1-ae)
    for widths in (0, 1, 2):
        interval = {k: (m-widths*d, m+widths*d) for k, (m,d) in POLE.items()}
        fb, sl, predictions = [], [], []
        for rb,re,ae,afb,ab in itertools.product(*interval.values()):
            fb.append(width_factor(rb,re,ae,4*afb/(3*ae)))
            sl.append(width_factor(rb,re,ae,ab))
            predictions.append(.75*ae*(1-1.10*(1-ae)/(3*rb*re)))
        out[str(widths)] = {"required_K_using_forward_backward": [min(fb),max(fb)],
                            "required_K_using_polarised_direct": [min(sl),max(sl)],
                            "AFB_lower_for_K_at_most_1p10": min(predictions)}
    # Table 7.9, lepton-universal extraction on the published sign branch.
    def ladder(grb,grl,gll):
        return (grb-grl/3)/(grl-gll)
    lower = min(ladder(rb,re,le) for rb,re,le in itertools.product(
        [.0962-.0126,.0962+.0126], [.23186-.00046,.23186+.00046],
        [-.26939-.00044,-.26939+.00044]))
    out["right_current_difference"] = {
        "positive_central":ladder(.0962,.23186,-.26939),
        "positive_lower_two_marginal_widths":lower,
        "sign_flipped_central":ladder(-.0962,.23186,-.26939)}
    assert out["2"]["required_K_using_forward_backward"][0]>1.1545
    assert out["2"]["AFB_lower_for_K_at_most_1p10"]>.0992+2*.0016
    assert out["2"]["required_K_using_polarised_direct"][0]<1.04
    return out


def representation_check():
    # Q=-1/3 states from arbitrary SU(2) irreps, including triplet/quartet.
    # Y is chosen separately in each representation to give this same charge.
    multiplets = [(0,0),(.5,.5),(.5,-.5),(1,1),(1,-1),(1.5,.5)]
    rng = np.random.default_rng(1633)
    for _ in range(64):
        v=rng.normal(size=len(multiplets))+1j*rng.normal(size=len(multiplets))
        weights=abs(v)**2/np.vdot(v,v).real
        down=sum(w*(j+m)*(j-m+1) for w,(j,m) in zip(weights,multiplets))
        up=sum(w*(j-m)*(j+m+1) for w,(j,m) in zip(weights,multiplets))
        avg=sum(w*m for w,(_,m) in zip(weights,multiplets))
        variance=sum(w*m*m for w,(_,m) in zip(weights,multiplets))-avg**2
        assert abs((down-up)-2*avg)<3e-15
        assert variance>=abs(avg)*(.5-abs(avg))-1e-14
    return 64


def scalar_witness():
    # One ordinary down-type doublet component, two extra singlets, one exotic
    # doublet component. Bare masses and Yukawa insertion reproduced separately.
    rng=np.random.default_rng(510)
    def random(shape):
        return rng.uniform(-1,1,size=shape)+1j*rng.uniform(-1,1,size=shape)
    competed=[]
    for _ in range(32):
        bare=np.zeros((4,4),dtype=complex)
        bare[1:3,:3]=.12*random((2,3))
        bare[1,1]+=1.5
        bare[2,2]+=2.
        bare[3,3]=2.5+.11j
        yuk=np.zeros((4,4),dtype=complex)
        yuk[0,:3]=.06*random(3)
        yuk[3,:3]=.30*random(3)
        for open_edge in (False,True):
            if open_edge:
                yuk[1:3,3]=.4*random(2)
            u,s,vh=np.linalg.svd(bare+yuk)
            l,r=u[:,-1],vh.conj().T[:,-1]
            p_s=float(np.vdot(l[1:3],l[1:3]).real)
            p_el, p_er=abs(l[3])**2,abs(r[3])**2
            tl=-.5*abs(l[0])**2+.5*p_el
            kappa=np.vdot(l,yuk@r)/s[-1]
            # Feynman-Hellmann from physical singular values.
            step=1e-5
            plus,minus=[np.linalg.svd(bare+(1+h)*yuk,compute_uv=False)[-1] for h in (step,-step)]
            assert abs((plus-minus)/(2*step*s[-1])-kappa.real)<1e-8
            predicted=1-2*(tl+.5)-p_er+2*(s[-1]/abs(bare[3,3]))**2*p_er
            if not open_edge:
                assert abs(kappa-(1-p_s-p_er))<2e-12
                assert abs(kappa-predicted)<2e-12
            else:
                competed.append(abs(kappa-predicted))
    assert max(competed)>.01
    return {"restricted_complex_mass_witnesses":32,
            "ordinary_allowed_extra_Yukawa_yields_counterexample":True}


def result():
    return {"pole_readouts":pole_readouts(), "ladder_witnesses":representation_check(),
            "Higgs_mass_derivative":scalar_witness()}


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    out=result()
    saved=HERE/"independent_results.json"
    if args.check:
        if out!=json.loads(saved.read_text(encoding="utf-8")):
            raise SystemExit("1633 independent results changed")
        print("1633 independent PASS: pole readouts, weak ladders and complex mass derivatives")
    else:
        with saved.open("x",encoding="utf-8") as f:
            f.write(json.dumps(out,indent=2)+"\n")
