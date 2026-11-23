"""Symmetric-square filter, five fixed singlet tests, and independently tagged weak currents."""
import argparse
import json
from math import asin, degrees, log, pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
Q = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])
embed = np.array([[1,0,0],[0,1/sqrt(2),0],[0,1/sqrt(2),0],[0,0,1]])


def read(path):
    return json.loads((HERE / path).read_text(encoding="utf-8-sig"))


def symmetric_square(k):
    return embed.T @ np.kron(k,k) @ embed


def probe(z):
    return np.array([1,sqrt(2)*z,z*z])/(1+abs(z)**2)


def probability(kraus, state):
    outputs = np.array([k @ state for k in kraus])
    return float(np.sum(abs(outputs @ Q @ outputs.T)**2)/3)


def ckm(s12,s13,s23,phase):
    c12,c13,c23=np.sqrt(1-np.array([s12,s13,s23])**2)
    e=np.exp(1j*phase)
    return np.array([[c12*c13,s12*c13,s13/e],
                     [-s12*c23-c12*s23*s13*e,c12*c23-s12*s23*s13*e,s23*c13],
                     [s12*s23-c12*c23*s13*e,-c12*s23-s12*c23*s13*e,c23*c13]])


def run():
    data=read("inputs.json")
    rng=np.random.default_rng(1441)
    n=data["probe_count"]
    roots=np.exp(2j*pi*np.arange(n)/n)
    for _ in range(32):
        k=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
        k/=np.linalg.norm(k,ord=2)
        t=symmetric_square(k)
        assert np.allclose(t.T@Q@t,np.linalg.det(k)**2*Q,atol=1e-13)
        expected=np.linalg.svd(k,compute_uv=False)
        assert np.allclose(np.linalg.svd(t,compute_uv=False),
                           [expected[0]**2,np.prod(expected),expected[1]**2])
        assert max(probability([t],probe(z)) for z in roots)<1e-26
        mats=rng.normal(size=(3,3,3))+1j*rng.normal(size=(3,3,3))
        mats/=sqrt(np.linalg.norm(sum(m.conj().T@m for m in mats),ord=2))
        cap=16*max(probability(mats,probe(z)) for z in roots)
        for z in [0,*(rng.normal(size=15)+1j*rng.normal(size=15))]:
            assert probability(mats,probe(z))<=cap+1e-14
        u=np.linalg.qr(k)[0]
        v=abs(symmetric_square(u))
        assert np.isclose(v[0,1],v[1,2]) and np.isclose(v[0,0],v[2,2])

    s=data["counterfeit_strength"]
    false_four=np.array([[1,0,s],[0,sqrt(1-s*s),0],[-s,0,1]])/sqrt(1+s*s)
    assert max(probability([false_four],probe(z)) for z in (1,1j,-1,-1j))<1e-26
    four=probability([false_four],probe(0))
    K=np.diag(data["isospectral_single_filter_singular_values"])
    good=symmetric_square(K)
    batched=np.diag([1,1j,1])@good
    assert np.allclose(np.linalg.svd(batched,compute_uv=False),np.linalg.svd(good,compute_uv=False))
    equal_spectrum=probability([batched],probe(1))
    # A unipotent map perturbs only the z^4 term. It attains factor 16 at the south pole.
    sharp=np.eye(3)
    sharp[0,2]=s
    sharp/=np.linalg.norm(sharp,ord=2)
    ratio=probability([sharp],np.array([0,0,1]))/max(probability([sharp],probe(z)) for z in roots)
    assert np.isclose(ratio,16)
    replacement=[np.outer(probe(.25j),r) for r in np.eye(3)]
    assert max(probability(replacement,probe(z)) for z in roots)<1e-26

    shared=read(data["product_at_MZ"])
    up,down=[np.array(shared["yukawa"][f])[:,0] for f in ("up","down")]
    V=ckm(*np.array(shared["quark_sines_and_phase"])[:,0])
    singular=np.linalg.svd(up[:,None]*V*down,compute_uv=False)
    geometric_defect=float(singular[0]*singular[2]/singular[1]**2)
    measurements=read(data["separate_semileptonic_inputs"])
    kb=measurements["K_product"]
    kf=measurements["K_form"]
    cb=measurements["F_times_cb"]
    two=[]
    for widths in (1,2):
        us=(kb[0]-widths*kb[1])/(kf[0]+widths*kf[1])
        cbmax=(cb[0]+widths*sum(cb[1:]))/measurements["F_lower"]
        two.append({"widths":widths,"Vus_min":us,"Vcb_max":cbmax,
                    "element_gap_min":us-cbmax,
                    "sum_of_sector_ray_errors_min_deg":degrees(asin(us)-asin(cbmax))/2,
                    "F1_max_to_restore_exact_equality":(cb[0]+widths*sum(cb[1:]))/us})
    return {"four_dark_pole_bright":four,"same_spectra_distinct_joint_bright":equal_spectrum,
            "sharp_probability_amplification":round(ratio,10),
            "weak_product_central_singular_values":singular.tolist(),
            "endpoint_product_over_middle_squared":geometric_defect,
            "least_spectral_log_change_central":abs(log(geometric_defect))/4,
            "common_left_tests":two}


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return np.isclose(a,b,rtol=3e-10,atol=2e-14)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    opt=parser.parse_args()
    result=run()
    if opt.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    elif opt.check:
        assert equal(result,read("results.json"))
        print("1441 conic filter, optimal five tests, geometric spectrum and weak overlap verified")
    else:
        print(json.dumps(result,indent=2))
