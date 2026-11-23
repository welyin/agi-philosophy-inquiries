"""Cofactor certificate against any exact I+II tensor cancellation."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE/"inputs.json").read_bytes())


def triangle(lower,upper):
    return max(0.,max(lower+upper)-sum(upper)),sum(upper)


def pencil_enclosure(data, width):
    u,d,e,v,_=data
    mag=abs(v)
    low=1-width
    high=1+width
    ml,mh=mag*low,np.minimum(1.,mag*high)
    measures={}
    for i in range(3):
        j,k=[x for x in range(3) if i!=x]
        lp=ml[:,j]*ml[:,k]
        hp=mh[:,j]*mh[:,k]
        # A=V^dagger diag(u) V^*, adj(A) has the reciprocal phases.
        amin,amax=triangle(u*low*lp,u*high*hp)
        cof=np.array([np.prod(np.delete(u,l)) for l in range(3)])
        cmin,cmax=triangle(cof*low**2*lp,cof*high**2*hp)
        measures[i]=(amin,amax,cmin,cmax)
    A=measures[1]
    B=measures[2]
    s_upper=A[3]/(d[1]*low*A[0]) if A[0]>0 else None
    b_lower=B[2]/(d[2]*high*B[1])
    # |minor| <= ||M||2 sigma2(M), with j,k,i distinct.
    cut=CFG["split_modulus_z"]
    low_z=max(0.,B[2]-cut*d[2]*high*B[1])/(u[2]*high+cut*d[2]*high)
    high_z=max(0.,cut*d[1]*low*A[0]-A[3])/(u[2]*high+cut*d[2]*high)
    ratio_min=(B[2]*A[0])/(A[3]*B[1])
    # Exceptional E=aD+bU. In U's diagonal basis bound D's light compression
    # and heavy-light column separately, excluding all light/small-|a| roots.
    ta=min(1.,math.hypot(mh[0,2],mh[1,2]))
    tb=min(1.,math.hypot(mh[2,0],mh[2,1]))
    comp=d[1]*high + d[2]*high*ta**2
    mix=d[2]*high*ta*mh[2,2]+d[1]*high*tb
    def mu_cap(a):
        beta=(e[2]*high+a*d[2]*high)/(u[2]*low)
        aa=a*comp+beta*u[1]*high
        col=a*mix
        c=e[2]*low-aa-col
        return aa+col**2/c if c>0 else math.inf
    left,right=0.,8.
    for _ in range(64):
        mid=(left+right)/2
        if mu_cap(mid)<e[1]*low:
            left=mid
        else:
            right=mid
    dct_min,_=triangle(d*low*ml[1,:]*ml[2,:],d*high*mh[1,:]*mh[2,:])
    abound=(e[2]*high)/dct_min if dct_min>0 else None
    bbound=(e[2]*high+abound*d[2]*high)/(u[2]*low) if abound is not None else None
    dsfrac=d[1]*high/(d[2]*low)
    bfloor=max(0.,e[1]*low-e[2]*high*dsfrac)/(u[2]*high*(1+dsfrac))
    # Homogenize the rank-one-pencil gap. For r >= far, Weyl suffices.
    # For cut<=r<=far, (h*r-k)/((T+B*r)*(1+r)) has no positive internal minimum.
    far=2*u[2]*high/(d[1]*low)
    def mid_gap(r):
        return max(0.,r*d[1]*low*A[0]-A[3])/((u[2]*high+r*d[2]*high)*(1+r))
    proj_gap=min(low_z/(1+cut),mid_gap(cut),mid_gap(far),
                 d[1]*low*.5*far/(1+far))
    wedge=0. if abound is None or left<=1 else min((left-1)/(left+2),bfloor/(abound+3))/4
    F_frob=high*(np.linalg.norm(d)+np.linalg.norm(e))/4
    F_op=high*(d[2]+e[2])/4
    explicit=proj_gap**2*wedge**2/(2*F_frob*F_op)
    return dict(relative_width=width,minor_input=[list(measures[i]) for i in range(3)],
                modulus_z_required_by_s_at_most=s_upper,modulus_z_required_by_b_at_least=b_lower,
                d_b_over_s_required_at_least=ratio_min,d_b_over_s_input_at_most=d[2]*high/(d[1]*low),
                lower_sigma2_U_minus_zD_over_all_complex_z=min(low_z,high_z),
                second_lepton_if_abs_a_1_at_most=mu_cap(1),muon_at_least=e[1]*low,
                every_fit_abs_a_greater_than=left,compact_a_upper=abound,b_interval=[bfloor,bbound],
                homogeneous_second_singular_floor=proj_gap,coefficient_wedge_floor=wedge,
                normalized_cancellation_floor=explicit)


def shape_at_point(old,data,x):
    U,D,E,F,_,_,_,_=old["matrices"](x,data)
    _,_,_,_,Y,I=old["symmetric_sources"](x,data)
    J=Y@np.linalg.solve(F,Y)*np.linalg.det(F)
    np.testing.assert_allclose(J,I*np.linalg.det(F),rtol=1e-10,atol=1e-20)
    f=F/np.linalg.norm(F)
    i=I/np.linalg.norm(I)
    rho=np.vdot(f,i)
    best=-rho.conjugate()/abs(rho)
    # Use direct residual for small sine values rather than loss of precision in 1-|rho|.
    floor=np.linalg.norm(f+best*i)/2
    np.testing.assert_allclose(floor**2,(1-abs(rho))/2,atol=3e-16)
    z=complex(x[7],x[8])
    actual=np.linalg.norm(f+z*i)/(1+abs(z))
    assert floor<=actual+1e-13
    return dict(smallest_attainable_relative_residual=float(floor),
                physical_relative_residual=float(actual),phase_overlap=float(abs(rho)),
                optimal_normalized_I_coefficient=[float(best.real),float(best.imag)])


def algebra_tests(data):
    rng=np.random.default_rng(23812)
    u,d,_,v,_=data
    min_sigma=math.inf
    for _ in range(32):
        tagged=np.exp(1j*rng.uniform(-math.pi,math.pi,3))[:,None]*v*np.exp(1j*rng.uniform(-math.pi,math.pi,3))[None,:]
        a=tagged.conj().T@np.diag(u)@tagged.conj()
        z=rng.uniform(-25,25)+1j*rng.uniform(-25,25)
        t=a-z*np.diag(d)
        spec=np.linalg.svd(t,compute_uv=False)
        min_sigma=min(min_sigma,spec[1])
        co=np.linalg.det(a)*np.linalg.inv(a)
        for i in range(3):
            j,k=[q for q in range(3) if q!=i]
            const=a[i,i]*a[j,k]-a[i,j]*a[i,k]
            np.testing.assert_allclose(const,-co[j,k],rtol=1e-10,atol=1e-17)
            minor=t[i,i]*t[j,k]-t[i,j]*t[i,k]
            np.testing.assert_allclose(minor,const-z*d[i]*a[j,k],rtol=1e-12,atol=1e-17)
            assert abs(minor)<=spec[0]*spec[1]+1e-14
    return min_sigma


def calculate():
    old=runpy.run_path(str(HERE/CFG["charged_parent"]))
    data=old["backend"]()
    orig=json.loads((HERE/CFG["frozen_coordinates"]).read_bytes())
    boxes=[pencil_enclosure(data,w) for w in CFG["relative_widths"]]
    for row in boxes[:4]:
        assert row["modulus_z_required_by_b_at_least"]>row["modulus_z_required_by_s_at_most"]
        assert row["every_fit_abs_a_greater_than"]>1
    assert boxes[0]["lower_sigma2_U_minus_zD_over_all_complex_z"]>0
    return dict(phase_independent_boxes=boxes,
                crosscheck_random_example_minimum_second_singular=algebra_tests(data),
                exact_charged_and_mixed_witnesses={k:shape_at_point(old,data,np.array(x)) for family in
                  ["mixed_type_I_II","fixed_delta_minus70_mass_readout"] for k,x in orig[family].items()})


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--check",action="store_true");args=ap.parse_args()
    out=json.dumps(calculate(),indent=2,allow_nan=False)+"\n"
    dest=HERE/"results.json"
    if args.check:
        assert json.loads(dest.read_bytes())==json.loads(out)
    else:
        dest.write_text(out,encoding="utf-8")
    print("1801 OK",hashlib.sha256(out.encode()).hexdigest()[:12])
