"""Independent mass-eigenstate CP test and entry-cofactor test, round 1638."""
import argparse
import cmath
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent

def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))

def unpack(x):
    y=np.asarray(x)
    return y[...,0]+1j*y[...,1]

def real_down_basis(w):
    """In the down mass eigenbasis make all four distinguished gauge directions real."""
    l,s,r,t=[unpack(w[key]) for key in ("left_D","left_S","right_D","right_S")]
    rephase=np.ones(len(r),dtype=complex)
    rephase[1]=r[1].conjugate()/abs(r[1])
    real_vectors=[]
    for v in (s,l,t,r):
        phased=v*rephase
        assert max(abs(phased.imag))<1e-13
        real_vectors.append(phased.real)
    bases=[]
    for side in (real_vectors[:2],real_vectors[2:]):
        distinguished=np.column_stack(side)
        qr,_=np.linalg.qr(distinguished,mode="complete")
        basis=np.column_stack((qr[:,2:],distinguished))
        np.testing.assert_allclose(basis.T@basis,np.eye(5),atol=2e-11)
        bases.append(basis)
    diagonal=np.diag(w["physical_masses_GeV"])
    mass=bases[0].T@diagonal@bases[1]
    assert max(abs(mass[:4,4]))<1e-7 and abs(mass[4,3])<1e-7
    return max(abs(mass[:4,4]))

def test_charges():
    count=0
    # In Z_n every surviving path q-A-d, q-B-S, D-C-d forces D=S-h.
    for n in (2,3,4,5):
        for h,q,d,s,D in itertools.product(range(n),repeat=5):
            if (q-d-h)%n or (q-s-h)%n or (D-d+h)%n:
                continue
            assert (s-D-h)%n==0
            count+=1
    return count

def signed_window(c,z):
    t=math.tan(c)
    return [-t/(z.imag+t*z.real),t/(z.imag-t*z.real)]

def point(w,vev,angle,neutron):
    real_down_basis(w)
    M=unpack(w["weak_basis_mass_GeV"]); S,D=3,4
    m=M[:3,:3]
    # Admissible vertex obtained with the up-sector top doublet direction.
    tagged=(M[D,:3]@M[2,:3].conjugate())*M[2,S]/vev**3
    down=sum((M[D,:3]@M[i,:3].conjugate())*M[i,S] for i in range(3))/vev**3
    phase=M[S,S]*M[D,D]/abs(M[S,S]*M[D,D])
    minor=np.linalg.det(np.delete(np.delete(M,S,axis=0),D,axis=1))
    cofactor=(-1)**(S+D)*minor/np.linalg.det(M)
    slope=vev*phase*tagged.conjugate()*cofactor
    down_slope=vev*phase*down.conjugate()*cofactor
    assert abs(down_slope.imag)<2e-13 and abs(slope.imag)>1e-7
    ratio=lambda b:1+b*slope
    # Direct determinant of the filled five-state matrix at finite amplitudes.
    for b in (-.7,.2,.8):
        modified=M.copy()
        modified[S,D]=modified[S,D]+b*vev*phase*tagged.conjugate()
        np.testing.assert_allclose(np.linalg.det(modified)/np.linalg.det(M),ratio(b),atol=1e-13,rtol=1e-13)
    bounds={}
    for name,cap in (("working",angle),("neutron",neutron)):
        lo,hi=signed_window(cap,slope)
        assert lo<0<hi and np.isclose(cmath.phase(ratio(lo)),-cap,rtol=3e-9,atol=0)
        assert np.isclose(cmath.phase(ratio(hi)),cap,rtol=3e-9,atol=0)
        bounds[name]=[lo,hi]
    theta=abs(cmath.phase(ratio(.001)))
    return {"ordinary_mass_min_singular":np.linalg.svd(m,compute_uv=False)[-1],
            "down_only_CP_basis":True,"top_slope":[slope.real,slope.imag],
            "signed_response_bounds":bounds,"theta_at_0p001":theta,
            "passes_theta_benchmark":theta<=angle,"passes_conditional_neutron":theta<=neutron}

def run():
    inputs=read("inputs.json"); states=read(inputs["matrices"])
    neutron=read("../1579/independent_results.json")["neutron_theta_upper_in_declared_matching"]
    return {"allowed_ab_charge_paths_scanned":test_charges(),"neutron_theta_ceiling":neutron,
            "fixed_points":{row:point(states[row],inputs["v_GeV"]/math.sqrt(2),inputs["working_angle"],neutron)
                            for row in inputs["keys"]}}

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    res=json.loads(json.dumps(run()),parse_float=lambda f:float(f"{float(f):.10g}"))
    if args.check:
        assert res==read("independent_results.json")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf8") as out:
            json.dump(res,out,indent=2)
            out.write("\n")
    print("1638 independent PASS: real down axes, cofactor, signed response bounds")
