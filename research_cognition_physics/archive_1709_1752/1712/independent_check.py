"""Extract the physical family identification required by a reciprocal down matrix."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def svd(a):
    l,s,rh=np.linalg.svd(a)
    return l[:,::-1],s[::-1],rh.conj().T[:,::-1]


def tagged_matrices(flow, rep, step):
    y=np.zeros((4,3,4),complex)
    y[0,:,:3]=np.diag(flow["central"]("up"))
    y[1,:,:3]=flow["low_ckm"](rep["gamma_MZ_deg"],flow["angles"]("direct"))@np.diag(flow["central"]("down"))
    le=rep["charged_left_MZ"]
    y[2,:,:3]=(np.array(le["real"])+1j*np.array(le["imag"]))@np.diag(flow["central"]("charged_leptons"))
    poles=flow["fields"]["neutral_poles_GeV"]
    for a,susy,dt,g in flow["segments"](step):
        if a==poles[2]:
            at_cluster=y[:2,:,:3].copy()
        if a==flow["MS"]:
            tan=flow["TAN"]
            y/=np.array([tan,1.,1.,tan])[:,None,None]/math.hypot(tan,1)
        for j,m in enumerate(poles):
            if m==a:
                y[3,:,j]=flow["physics"]["Y0"][:,j]
        for i in range(0,len(g)-1,2):
            beta=flow["beta"]
            k1=beta(y,g[i],susy)
            k2=beta(y+dt*k1/2,g[i+1],susy)
            k3=beta(y+dt*k2/2,g[i+1],susy)
            k4=beta(y+dt*k3,g[i+2],susy)
            y+=dt*(k1+2*k2+2*k3+k4)/6
    Lu,_,_=svd(y[0,:,:3])
    Ld,_,Rd=svd(y[1,:,:3])
    physical=rep["down_texture_10PeV"]
    target=np.array(physical["real"])+1j*np.array(physical["imag"])
    T,_,R=svd(target)
    # Fix row and column phases by the nonzero CKM first row and first column.
    phases=(Lu.conj().T@Ld)/T
    np.testing.assert_allclose(abs(phases),1.,atol=3e-8)
    phases/=abs(phases)
    left=Lu*(phases[:,0]/phases[0,0])
    right=(Rd*phases[0,:].conj())@R.conj().T
    np.testing.assert_allclose(left.conj().T@y[1,:,:3]@right,target,atol=3e-10)
    pu=left.conj().T@svd(at_cluster[0])[0][:,0]
    down=left.conj().T@at_cluster[1]@right
    return down,pu,target


def inspect(down,pu):
    # Flow uses the holomorphic Weyl texture. Taking its complex conjugate as the
    # four-component down Yukawa conjugates J but leaves all moduli and |quartets| invariant.
    L,s,R=svd(down)
    qbar=np.array([2,0,1])[:,None]
    qten=np.array([0,0,1])[None,:]
    summed=qbar+qten
    out=[]
    for strange,bottom in itertools.product((1,-1),repeat=2):
        signs=np.array([-1,strange,bottom])
        j=(R*signs)@L.conj().T
        t=down@j
        np.testing.assert_allclose(j.conj().T@j,np.eye(3),atol=2e-14)
        np.testing.assert_allclose(t,t.conj().T,atol=2e-16)
        w=abs(j)**2
        perms=itertools.permutations(range(3))
        best=max(sum(abs(j[r,col]) for r,col in enumerate(p)) for p in perms)
        quartet=(j[0,0]*j[1,1]*j[0,1].conj()*j[1,0].conj()).imag
        out.append(dict(signs=signs.tolist(),return_over_yd=np.vdot(pu,t@pu).real/s[0],
                        abs_J=abs(j).tolist(),physical_cross_charge_fraction=w[(summed%12)!=0].sum()/3,
                        min_frobenius_distance_to_rephased_permutation=math.sqrt(max(0.,6-2*best)),
                        min_single_insertion_Wilson_max=float(np.max(abs(j)/.5**summed)),
                        abs_quartet=abs(quartet)))
    return out


def examine(case):
    generator=runpy.run_path(str(HERE/"check.py"))
    fl=generator["FL"]
    solver=fl["solve"]
    g={}
    def capture(*args,**kwargs):
        rep=solver(*args,**kwargs)
        g["rep"]=rep
        return rep
    fl["solve"]=capture
    test=generator["run_point"](**case)
    down,pu,target=tagged_matrices(fl,g["rep"],.20)
    rows=inspect(down,pu)
    physical=min(rows,key=lambda row:abs(row["return_over_yd"]))
    assert physical["signs"]==[-1,1,case["sign"]] and abs(physical["return_over_yd"])<1e-7
    assert all(row["physical_cross_charge_fraction"]>=2/3-1e-12 for row in rows)
    # A left-only Gram cannot see a relative right-source sign. With both
    # chiral ports and three resolved singular values the signed polar factor is explicit.
    L,s,R=svd(down)
    signs=np.array(physical["signs"])
    j=(R*signs)@L.conj().T
    alt=(L*(s*[-1,1,1]))@R.conj().T
    alt_j=(R*(signs*[-1,1,1]))@L.conj().T
    np.testing.assert_allclose(alt@alt.conj().T,down@down.conj().T,atol=1e-17)
    np.testing.assert_allclose(np.linalg.norm(j-alt_j,2),2.,atol=1e-13)
    np.testing.assert_allclose(np.linalg.norm(down-alt,2),2*s[0],rtol=5e-12)
    gram=down@down.conj().T
    poly=np.zeros((3,3),complex)
    for i in range(3):
        pr=np.eye(3,dtype=complex)
        for k in range(3):
            if i != k:
                pr=pr@(gram-s[k]**2*np.eye(3))/(s[i]**2-s[k]**2)
        poly+=signs[i]*pr/s[i]
    residual=np.linalg.norm(down.conj().T@poly-j,2)
    assert residual<2e-8
    uv=inspect(target,np.array([1.,0.,0.]))
    return dict(input=case,predicted_yd=test["y_d_MZ"],
                required_identification=physical,all_roots=rows,
                polar_polynomial_error=residual,
                signed_polar_pulse_lower_for_error_1e_minus4=(1-1e-4)/s[0],
                least_common_spectral_resolution_factor=1/s[0],
                same_signs_refitted_10PeV_J=next(row for row in uv if row["signs"]==physical["signs"]),
                original_texture_10PeV_abs_23_and_32=[abs(target[1,2]),abs(target[2,1])])


def compute():
    return [examine(row) for row in CASES]


def stable(x):
    if isinstance(x,list):
        return [stable(v) for v in x]
    if isinstance(x,dict):
        return {k:stable(v) for k,v in x.items()}
    if isinstance(x,(float,np.floating)):
        if abs(x)<1e-12:
            return 0.
        return float(f"{x:.9g}")
    return x


CASES=[dict(sign=-1),dict(sign=1),
       dict(sign=-1,us=-.42,strange=-.93,bottom=.93,k=1,ub=-.93,cb=.93),
       dict(sign=-1,us=.75,strange=-1.8,bottom=.4,k=1,ub=-.95,cb=.8)]
if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    value=stable(compute())
    target=HERE/"independent_results.json"
    if cli.parse_args().check:
        assert json.loads(target.read_bytes())==value
        print("Verified reciprocal unitary ports, transported texture, flavour charges and quartets.")
    else:
        with target.open("x",encoding="utf8") as fh:
            fh.write(json.dumps(value,indent=2)+"\n")
