"""Independent one-direction bound on an arbitrary complex linear mass matrix."""
import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
A=np.diag([1.,1.,0.])
B=np.diag([2.,1.,1.])


def witness(limit):
    mass=1.
    z=1-limit/2
    left=np.array([sqrt(z/2),0,sqrt(1-z/2)])
    right=np.array([sqrt(z/2),sqrt(1-z/2),0])
    left_basis=np.column_stack([left,[-left[2],0,left[0]],[0,1,0]])
    right_basis=np.column_stack([right,[-right[1],right[0],0],[0,0,1]])
    target=left_basis @ np.diag([mass,.07,.003]) @ right_basis.T
    gradient=np.outer(A@left,right)+np.outer(left,B@right)
    u,sv,vh=np.linalg.svd(gradient)
    floor=sqrt(5-2*limit)/2
    direct=floor*(u@vh)
    high=target-direct
    powers=1+2*(np.diag(A)[:,None]+np.diag(B)[None,:])
    def combined(x):
        return high*x**powers+direct*x
    slope=left@(combined(1+1e-5)-combined(1-1e-5))@right/(2e-5)
    assert np.isclose(slope,limit,rtol=1e-8)
    assert np.isclose(np.linalg.norm(direct,ord=2),floor)
    assert np.isclose(sv.sum(),sqrt(1+4*z))
    assert np.allclose(np.linalg.svd(combined(1),compute_uv=False),[1,.07,.003])
    # Uncorrected factorizable source alone: nine diagonal matrix derivatives, real.
    rng=np.random.default_rng(1435)
    for _ in range(20):
        coeff=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
        u0,s0,v0=np.linalg.svd(coeff)
        diag=np.diag(u0.conj().T@(coeff*powers)@v0.conj().T)/s0
        assert max(abs(diag.imag))<1e-10 and min(diag.real)>3-1e-10
        assert np.isclose(diag.sum(),15)
    return {"max_abs_kappa":limit,"least_linear_mass_operator_norm_over_measured_mass":floor,
            "saturating_projector_weight_each_side":z/2}


def run():
    widths=[(20.5,8.4,22.,.51),(2.04,.62,2.369,.017)]
    return {label:witness(sqrt((v+2*e)/(.9*(sm-2*se))))
            for label,(v,e,sm,se) in zip(["bottom","tau"],widths)}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as output:
            json.dump(result,output,indent=2)
            output.write("\n")
    else:
        saved=json.loads(path.read_text(encoding="utf-8"))
        assert saved.keys()==result.keys()
        for group, row in result.items():
            assert row.keys()==saved[group].keys()
            assert all(np.isclose(v,saved[group][k],rtol=1e-10) for k,v in row.items())
        print("1435 published ratio endpoints and sharp complex linear mass repair checked")
