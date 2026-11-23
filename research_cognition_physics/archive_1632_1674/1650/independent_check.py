"""Independent exact pair cancellation and least-squares lower bound for all flavor pairs."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
I,J=np.nonzero(~np.eye(3,dtype=bool))


def read(path):
    return json.loads((HERE/path).read_bytes())


def mass():
    src=read("../1489/results.json")["reconstruction"]
    a,b,c,e,f=src["normalized_real_edges_A_B_C_E_F"]
    real,amp=src["two_real_coefficients"]
    ph=amp*np.exp(1j*math.radians(src["eta_example_deg"]))
    return np.array([[0,a,0],[b,c,real+ph],[0,e,f]]),real,ph


def embed(y):
    L,m,Rh=np.linalg.svd(y)
    idx=np.argsort(m)
    return L[:,idx],m[idx],Rh.conj().T[:,idx]


def flatten(z):
    v=z[I,J]
    return np.r_[v.real,v.imag]


def recover_real_endpoints(U,target,i,j):
    # Fix the irrelevant constant port by setting the third coefficient to zero.
    row=U[:2,i].conjugate()*U[:2,j]
    mat=np.vstack([row.real,row.imag])
    v=np.linalg.solve(mat,[target.real,target.imag])
    np.testing.assert_allclose(row@v,target,rtol=3e-12,atol=2e-13)
    return np.r_[v,0.],np.linalg.cond(mat),np.linalg.det(mat)


def widths(mat,m):
    return {name:(abs(mat[i,j])**2+abs(mat[j,i])**2)/m[2]**2
            for (name,i,j) in [("ds",0,1),("db",0,2),("sb",1,2)]}


def texture_test(y,real,part,qr,qd):
    L,m,R=embed(y)
    dl=np.diag(m)
    x=np.zeros_like(y)
    x[1,2]=qr*real+qd*part
    g=L.conj().T@x@R
    p0=g.diagonal().imag/m
    i,j=0,1
    mi,mj=m[i],m[j]
    den=mj**2-mi**2
    tl=(-mj*g[i,j]+mi*g[j,i].conjugate())/den
    tr=(mi*g[i,j]-mj*g[j,i].conjugate())/den
    l,cl,jl=recover_real_endpoints(L,tl,i,j)
    r,cr,jr=recover_real_endpoints(R,tr,i,j)
    def shift(l,r):
        return L.conj().T@(np.diag(l)@y+y@np.diag(r))@R
    pair=g+shift(l,r)
    assert max(abs(pair[i,j]),abs(pair[j,i]))<2e-11
    np.testing.assert_allclose(pair.diagonal().imag/m,p0,atol=2e-10)
    # QR projection on real and imaginary pieces, with a real port convention.
    # Only four coefficients affect the off-diagonal observations; constant L/R adds a real mass.
    basis=np.eye(3)
    zeros=np.zeros(3)
    directions=[shift(a,zeros) for a in basis[:2]]+[shift(zeros,a) for a in basis[:2]]
    design=np.column_stack([flatten(a) for a in directions])
    assert np.linalg.matrix_rank(design)==4
    qmat,_=np.linalg.qr(design)
    yvec=flatten(g)
    opt=-np.linalg.lstsq(design,yvec,rcond=None)[0]
    optimal=g+shift(np.r_[opt[:2],0],np.r_[opt[2:],0])
    resid=flatten(optimal)
    np.testing.assert_allclose(resid,yvec-qmat@(qmat.T@yvec),atol=2e-13)
    np.testing.assert_allclose(design.T@resid,np.zeros(4),atol=2e-13)
    np.testing.assert_allclose(optimal.diagonal().imag/m,p0,atol=2e-10)
    # Direct discriminant of the right Gram matrix provides an independent orientation check.
    left,right=y@y.conj().T,y.conj().T@y
    oriented=lambda h:(h[0,1]*h[1,2]*h[2,0]).imag
    a,b,c,e,f=y[0,1].real,y[1,0].real,y[1,1].real,y[2,1].real,y[2,2].real
    ring=c*e*f*y[1,2].imag
    np.testing.assert_allclose([oriented(left),oriented(right)],[a*a*ring,-b*b*ring],rtol=1e-11,atol=1e-15)
    assert abs(jl)>1e-12 and abs(jr)>1e-12
    sources=[]
    for coef in (-real,real+part):
        ins=np.zeros_like(y);ins[1,2]=coef
        sources.append(flatten(L.conj().T@ins@R))
    T=np.column_stack(sources)
    projected=T-qmat@(qmat.T@T)
    gram=projected.T@projected/m[2]**2
    cvec=np.array([qd-qr,qd])
    minimum=np.dot(resid,resid)/m[2]**2
    np.testing.assert_allclose(cvec@gram@cvec,minimum,atol=1e-13)
    eig=np.linalg.eigvalsh(gram)
    assert eig[0]>0
    pair_width=widths(pair,m)
    yb0,dyb=read("../1386/inputs.json")["yukawa"]["down"][2]
    hbs_lower=(yb0-2*dyb)*math.sqrt(pair_width["sb"])/2
    assert hbs_lower>math.sqrt(4*math.pi)
    return {"qR_qD":[qr,qd],"p":p0.tolist(),"ordered_left_right_quartet_minors":[jl,jr],
            "condition_numbers_for_ds_cancellation":[cl,cr],
            "pair_cancel_left_right_ports":[l.tolist(),r.tolist()],
            "cancelled_pair_ij_ji_abs":[abs(pair[i,j]),abs(pair[j,i])],
            "widths_base":widths(g,m),"widths_after_ds_cancel":pair_width,
            "conditional_max_of_physical_hbs_hsb_abs_lower":hbs_lower,
            "sum_after_ds_cancel_over_SM_bb":sum(pair_width.values()),
            "optimal_ports_l1_l2_r1_r2":opt.tolist(),"widths_at_joint_minimum":widths(optimal,m),
            "minimum_all_FC_over_SM_bb":minimum,"all_contrasts_minimum_Gram":gram.tolist(),
            "all_contrasts_minimum_eigenvalues":eig.tolist()}


def stable(x):
    if isinstance(x,dict):return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):return [stable(v) for v in x]
    if isinstance(x,float):return 0. if abs(x)<1e-16 else float(f"{x:.12g}")
    return x


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    y,r,d=mass()
    out=stable({"phase_source":texture_test(y,r,d,0,2),"real_source":texture_test(y,r,d,2,0)})
    if args.check:assert out==read("independent_results.json")
    else:(HERE/"independent_results.json").write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    for key,v in out.items():print(key,v["widths_after_ds_cancel"],v["minimum_all_FC_over_SM_bb"])

if __name__=="__main__":
    main()
