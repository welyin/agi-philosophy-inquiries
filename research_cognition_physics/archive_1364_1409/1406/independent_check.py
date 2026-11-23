"""Chirality comparison using a second, vector-meson final state."""
import argparse
from math import sqrt
from pathlib import Path
import json
import numpy as np

HERE=Path(__file__).resolve().parent


def norm_certificate(k,star,k_sm,star_sm):
    u=sqrt(k/k_sm)
    cap=sqrt(star/star_sm)
    return dict(required_R_K=u*u,vector_rate_cap=cap*cap,
                pure_left_vector_prediction_1e6=u*u*star_sm,
                right_over_left_norm_lower=max(0,(u-cap)/(u+cap)))


def calculate():
    # Branching fractions scaled by 10^6. Belle limit assumes its nominal spectrum.
    one=norm_certificate(23-7,18,4.29+.13+.2,8.24-.99)
    two=norm_certificate(23-14,18,4.29+2*(.13+.2),8.24-2*.99)
    assert one["pure_left_vector_prediction_1e6"]>18>two["pure_left_vector_prediction_1e6"]
    # Matrix identity for flavor-inclusive vector/axial components.
    rng=np.random.default_rng(1406)
    for _ in range(50):
        L,R=rng.normal(size=(2,3,3))+1j*rng.normal(size=(2,3,3))
        assert np.isclose(np.linalg.norm(L+R)**2-np.linalg.norm(L-R)**2,4*np.vdot(L,R).real)
        u,v=np.linalg.norm(L+R),np.linalg.norm(L-R)
        assert abs(u-v)/(u+v)<=np.linalg.norm(R)/np.linalg.norm(L)+1e-14
    return {"central":norm_certificate(23,18,4.29,8.24),"one_width_favourable":one,"two_width_favourable":two}


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    result=calculate()
    if p.parse_args().check:
        old=json.loads((HERE/"independent_results.json").read_text(encoding="utf8"))
        assert old.keys()==result.keys()
        for k,section in result.items():
            assert old[k].keys()==section.keys() and all(np.isclose(a,old[k][b]) for b,a in section.items())
        print("1406 second invisible mode, chirality norm and precision window reproduced")
    else:
        print(json.dumps(result,indent=2))
