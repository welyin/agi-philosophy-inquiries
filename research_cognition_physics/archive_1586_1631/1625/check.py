"""Test an ordered diagonal plus a single flavoured column in the oscillation tensor."""
import argparse
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def matrix(x,y,z,delta):
    s,c=math.sqrt(x),math.sqrt(1-x)
    a,b=math.sqrt(z),math.sqrt(1-z)
    p=math.sqrt(y)*np.exp(1j*delta)
    v=math.sqrt(1-y)
    return np.array([[c*v,s*v,p.conjugate()],[-s*b-c*a*p,c*b-s*a*p,a*v],
                     [s*a-c*b*p,-c*a-s*b*p,b*v]])


def certificate(box):
    x0,x1=box["sin2_theta12"]
    y0,y1=box["sin2_theta13"]
    z0,z1=box["sin2_theta23"]
    r0,r1=box["solar_over_abs31"]
    tmin=math.sqrt(min(z0,1-z1))
    tmax=math.sqrt(max(z1,1-z0))
    s_min,s_max=math.sqrt(x0),math.sqrt(x1)
    c_min,c_max=math.sqrt(1-x1),math.sqrt(1-x0)
    pmax=math.sqrt(y1)
    vmin=math.sqrt(1-y1)
    # Bounds for |U_{mu/tau,1}|, |U_{mu/tau,2}|, |U_{mu/tau,3}|.
    u1=s_min*tmin-c_max*tmax*pmax
    u2=c_min*tmin-s_max*tmax*pmax
    u3=tmin*vmin
    er_low=vmin*(math.sqrt(y0)*tmin*(1-r1*x1)-r1*tmax/2)
    er_high=pmax*tmax*(1+r1*x1)+r1*tmax/2
    mut_low=(1-y1)*math.sqrt(z0*(1-z1))-r1/2
    mut_high=(1+r1)/2
    assert min(u1,u2,u3,er_low,mut_low)>0
    gap=r0*(1-r1)*u1*u2*u3/(er_high*mut_high)
    return {"lower_each_abs_row_entry_mu_tau":[u1,u2,u3],
            "abs_H_e_mu_tau_over_A":[er_low,er_high],"abs_H_mutau_over_A":[mut_low,mut_high],
            "forced_diagonal_gap_over_A_at_least":gap,
            "forced_diagonal_gap_eV2_at_least":gap*box["abs31_eV2"][0]}


def inferred_diagonal(H):
    vals=[]
    for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
        vals.append(H[i,i]-H[i,j]*H[k,i]/H[k,j])
    return np.array(vals)


def contains(box,value):
    return box[0]<=value<=box[1]


def run():
    cfg=read(HERE/"inputs.json")
    box=cfg["box"]
    reactor=read(HERE/cfg["reactor_inputs"])
    old=read(HERE/cfg["weak_inputs"])
    amp,er,_=reactor["DayaBay_sin2_2theta13"]
    solar,lo,hi=reactor["JUNO_sin2_theta12_NO"]
    for k in (-1,1):
        assert contains(box["sin2_theta13"],(1-math.sqrt(1-amp-k*2*er))/2)
    assert contains(box["sin2_theta12"],solar-2*lo) and contains(box["sin2_theta12"],solar+2*hi)
    for order in ("NO","IO"):
        val=old[order]
        gap21,eps21=val["solar_gap"]
        gapA,epsA=val["atmospheric_gap"]
        # NO uses |dm31|=dm32+dm21; IO uses the quoted |dm31| directly.
        for sign1 in [-1,1]:
            for sign2 in [-1,1]:
                dm21=gap21+sign1*2*eps21
                dm31=gapA+sign2*2*epsA+(dm21 if order=="NO" else 0.)
                assert contains(box["abs31_eV2"],dm31)
                assert contains(box["solar_over_abs31"],dm21/dm31)
        row=val["atmospheric_sine_squared"]
        assert contains(box["sin2_theta23"],row[0]-2*row[1])
        assert contains(box["sin2_theta23"],row[0]+2*row[-1])
    joint=cfg["T2K_NOvA_marginal"]
    assert contains(box["sin2_theta23"],joint["sin2_theta23"]-2*joint["lower_error"])
    assert contains(box["sin2_theta23"],joint["sin2_theta23"]+2*joint["upper_error"])
    y=(1-math.sqrt(1-amp))/2
    measure=cfg["centers"]
    solar_gap=measure["solar_gap_eV2"]
    atm=measure["abs31_eV2"]
    out={"analytic_box":certificate(box),"central_examples":{}}
    family=read(HERE/cfg["family_inputs"])["pole_MeV"]
    tau=read(HERE/cfg["family_output"])["roots_MeV"][1]
    D=np.diag(np.array([family["electron"][0],family["muon"][0],tau])/tau)
    inv_sqrt=np.diag(1/np.sqrt(np.diag(D)))
    for order,sgn in (("NO",1),("IO",-1)):
        ev=np.array([0.,solar_gap,sgn*atm])
        record={}
        for flag,phase in [("+",0.),("-",math.pi),("complex",1.29)]:
            U=matrix(solar,y,joint["sin2_theta23"],phase)
            assert np.allclose(U@U.conj().T,np.eye(3),atol=1e-13)
            h=(U*ev)@U.conj().T
            h=(h+h.conj().T)/2
            cp=(h[0,1]*h[1,2]*h[2,0]).imag
            diag=inferred_diagonal(h)
            if flag!="complex":
                assert np.max(np.abs(diag.imag))<1e-15
                d=diag.real
                for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:
                    assert abs(d[i]+ev[1]*ev[2]*U[j,0]*U[k,0].conjugate()/h[j,k])<1e-15
                assert (d[1]-d[0])*(d[2]-d[0])<0
                assert min(abs(d[1]-d[0]),abs(d[2]-d[0]))>=out["analytic_box"]["forced_diagonal_gap_eV2_at_least"]
                assert np.linalg.svd(h-np.diag(d),compute_uv=False)[1]<1e-16
                record[flag]={"unique_diagonal_eV2_relative_to_m1_sq":d.tolist()}
            else:
                assert abs(cp)>1e-12*atm**3
            H=h+(measure["lightest_eV"]**2-min(ev))*np.eye(3)
            b=float(np.linalg.eigvalsh(inv_sqrt@H@inv_sqrt)[0])
            vals,vec=np.linalg.eigh(H-b*D)
            assert b>0 and abs(vals[0])<2e-14 and vals[1]>0
            columns=vec[:,1:]*np.sqrt(vals[1:])
            error=float(np.max(np.abs(H-b*D-columns@columns.conj().T)))
            assert error<1e-13
            record.setdefault(flag,{}).update({"positive_common_scale_eV2":b,
                         "two_source_norms_eV":np.sqrt(vals[1:]).tolist(),
                         "matrix_reconstruction_error_eV2":error})
        out["central_examples"][order]=record
    return out


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-15)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    result=run()
    path=HERE/"results.json"
    if cli.parse_args().check:
        assert same(result,read(path))
        print("1625 PASS: continuous sign and gap bounds, both orders, one-column obstruction and two-column completion")
    else:
        with path.open("x",encoding="utf-8") as stream:
            stream.write(json.dumps(result,indent=2)+"\n")
