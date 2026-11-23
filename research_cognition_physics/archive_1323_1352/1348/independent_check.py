"""Independent finite SU(2) saturation family for the historical bottom-current bound."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np


def saturate(split, raise_strength, cb2, ub2):
    # One physical doublet plus arbitrary singlets, in the mass-labelled basis.
    visible = float(cb2+ub2)
    strength = float(raise_strength)
    hidden = strength-visible
    u = np.concatenate(([math.sqrt(float(ub2)),1j*math.sqrt(float(cb2))],
                        np.full(split, math.sqrt(hidden/split),dtype=complex)))
    u /= math.sqrt(strength)
    d = np.array([math.sqrt(strength),1j*math.sqrt(1-strength)])
    n = split+2
    tp = np.zeros((n+2,n+2),complex)
    tp[:n,n:] = np.outer(u,d.conj())
    tm = tp.conj().T
    t3 = np.zeros_like(tp)
    t3[:n,:n] = np.outer(u,u.conj())/2
    t3[n:,n:] = -np.outer(d,d.conj())/2
    charge=np.diag([2/3]*n+[-1/3]*2)
    err = max(np.linalg.norm(tp@tm-tm@tp-2*t3),
              np.linalg.norm(t3@tp-tp@t3-tp),
              np.linalg.norm(charge@tp-tp@charge-tp))
    raised = tp[:,n]
    hidden_components = np.abs(raised[2:n])**2
    assert err<5e-14 and abs(t3[n,n]+strength/2)<3e-15
    assert np.linalg.norm(tm[:,n])<3e-15
    assert abs(hidden_components.sum()-hidden)<4e-15
    assert np.allclose(hidden_components,hidden/split)
    return dict(partners=split,each_squared=hidden/split,total_squared=hidden)


def block_witness():
    # One doublet, one triplet, one singlet. All three tested states have
    # electric charge -1/3. Declare only the triplet's upper state calibrated.
    blocks = (2,3,1)
    size=sum(blocks)
    tp=np.zeros((size,size))
    t3=np.zeros_like(tp)
    offset=0
    for dim in blocks:
        j=(dim-1)/2
        for i in range(dim):
            m=-j+i
            t3[offset+i,offset+i]=m
            if i+1<dim:
                tp[offset+i+1,offset+i]=math.sqrt((j-m)*(j+m+1))
        offset+=dim
    z=np.eye(size)[:,[0,3,5]]
    angle=.43
    rotation=np.array([[math.cos(angle),math.sin(angle),0],[-math.sin(angle),math.cos(angle),0],[0,0,1]])
    z=z@rotation
    full=tp@z
    observed=np.diag([0,0,0,0,1,0])@full
    h=full-observed
    d=-2*z.T@t3@z-observed.T@observed
    lower=tp.T@z
    assert np.linalg.norm(h.T@h-d-lower.T@lower)<2e-14
    eigen=np.linalg.eigvalsh(d)
    assert np.allclose(eigen,[-2,0,1]) and np.linalg.matrix_rank(h)==1
    assert abs(np.trace(h.T@h)-max(eigen))<1e-14
    return {"observed_energy_rank":1,"unseen_subspace_rank":1,
            "certificate_eigenvalues":[-2,0,1],"certificate_trace":-1,
            "positive_trace_bound":1,"saturation":True}


def main():
    # Endpoints of the stated working rectangle, no unitarity normalization.
    strength=F(95,100)
    p=F(396,10000)
    f=F(4,5)
    ratio2=F(1,20)
    cb2=(p/f)**2
    ub2=cb2*ratio2
    h2=strength-cb2-ub2
    critical_f2=p*p*(1+ratio2)/strength
    trials=[saturate(k,strength,cb2,ub2) for k in (1,2,7,19)]
    # An operator identity error, a neutral-current expectation error, and a
    # joint vector mismatch for the two calibrated charged amplitudes.
    dz,dv,da=.02,.01,.005
    robust=float(strength)-2*dz-da-(math.sqrt(float(cb2+ub2))+dv)**2
    assert robust>.89
    return dict(status="PASS",round=1348,historical_upper_visible_squared=str(cb2+ub2),
                lower_unseen_squared=str(h2),lower_unseen_amplitude=math.sqrt(float(h2)),
                minimal_form_factor_for_positive_lower=math.sqrt(float(critical_f2)),
                explicit_su2_examples=trials,matrix_inertia_example=block_witness(),
                example_mismatch_bound=dict(neutral_expectation=dz,visible_vector=dv,
                                            algebra_expectation=da,lower_unseen_squared=robust))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=main()
    path=Path(__file__).resolve().with_name("independent_results.json")
    if args.write_new:
        with path.open("x",encoding="utf-8",newline="\n") as f:
            f.write(json.dumps(out,indent=2)+"\n")
    else:
        assert json.loads(path.read_text(encoding="utf-8"))==out
    print(json.dumps(out))

