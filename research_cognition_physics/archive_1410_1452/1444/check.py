"""Two noncommuting operational records and a coherent implementation of the HS texture."""
import argparse
import json
from math import acos, asin, atan2, cos, degrees, pi, radians, sin, sqrt, tan
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
I=np.eye(2)
X=np.array([[0,1],[1,0]])
Y=np.array([[0,-1j],[1j,0]])
Z=np.diag([1,-1])
PAULI=(X,Y,Z)


def component_product(amplitude,word):
    # Forward and reverse matrix elements of the controlled two-state preparation.
    value=1+0j
    for use in word:
        value*=amplitude if use=="f" else amplitude.conjugate()
    return value


def texture(amplitude,real_c,pars):
    b,aa=pars["b"],pars["A0"]
    paths=[(0,0,real_c,"frfr"),(0,1,b,"fff"),(1,0,b,"rrr"),
           (1,1,b,"fr"),(1,2,aa,"ff"),(2,1,aa,"rr"),(2,2,1,"")]
    target=np.zeros((3,3),dtype=complex)
    for i,j,coefficient,word in paths:
        target[i,j]+=coefficient*component_product(amplitude,word)
    assert np.allclose(target,target.conj().T)
    return target


def run():
    cfg=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    a,b=Z,-Y
    assert np.allclose(a@b+b@a,0) and np.allclose(-1j*a@b,X)
    swap=(a+b)/sqrt(2)
    assert np.allclose(swap@a@swap,b) and np.allclose(swap@X@swap,-X)
    proj=(I+swap)/2
    ref=proj[:,0]/np.linalg.norm(proj[:,0])
    prep=np.column_stack([ref,[-ref[1].conjugate(),ref[0].conjugate()]])
    assert np.allclose(prep.conj().T@prep,I)
    r=ref[1]/ref[0]
    assert np.isclose(r,-1j*tan(pi/8))
    quarter=np.diag([np.exp(-1j*pi/4),np.exp(1j*pi/4)])
    qx=float(((quarter@ref).conj()@X@(quarter@ref)).real)
    # Negative result of the bisector measurement can be corrected by X.
    assert np.allclose(X@(I-proj)@X,proj)
    # Equal-rate amplitude damping toward the two positive endpoints.
    adir=np.array([0,0,1])
    bdir=np.array([0,-1,0])
    R=np.linalg.solve(np.eye(3)+(np.outer(adir,adir)+np.outer(bdir,bdir))/2,adir+bdir)
    mixed=(I+sum(w*p for w,p in zip(R,PAULI)))/2
    assert np.allclose(swap@mixed@swap,mixed)
    pur=float(np.trace(mixed@mixed).real)

    fit=cfg["HS_rounded_fit_table1"]
    amplitudes=fit["lambda0"]*prep[:,0]
    quarks=[texture(amplitudes[i],fit[key],fit) for i,key in [(1,"c_up"),(0,"c_down")]]
    for i in range(2):
        assert np.allclose(prep.conj().T[0,i],prep[i,0].conjugate())
    hu,hd=[m@m.conj().T for m in quarks]
    p1,p3=np.diag([1,0,0]),np.diag([0,0,1])
    loop=np.trace(p1@hu@p3@hd)
    assert abs(loop.real)<1e-19 and loop.imag<0
    assert np.isclose(abs(hu[0,2]/hd[0,2]),abs(r)**5,rtol=1e-13)
    left=[]
    for m in quarks:
        vals,vec=np.linalg.eigh(m)
        left.append(vec[:,np.argsort(abs(vals))])
    v=left[0].conj().T@left[1]
    J=float((v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)
    def angle(z):
        return float(abs(np.angle(z,deg=True)))
    angles=[angle(-v[2,0]*v[2,2].conjugate()/(v[0,0]*v[0,2].conjugate())),
            angle(-v[1,0]*v[1,2].conjugate()/(v[2,0]*v[2,2].conjugate())),
            angle(-v[0,0]*v[0,2].conjugate()/(v[1,0]*v[1,2].conjugate()))]
    assert np.isclose(sum(angles),180)
    raw=cfg["tolerance_example"]
    err=radians(raw["separation_error_degrees"])/2
    phase=radians(raw["azimuth_error_degrees"])
    cap=2*asin(sqrt(raw["minus_bisector_probability_upper"]))
    low,high=pi/4-err,pi/4+err
    mod_range=[tan((low-cap)/2),tan((high+cap)/2)]
    azimuth=phase+asin(sin(cap)/sin(low))
    assert low>cap and high+cap<pi/2
    # Directly sample cap boundaries at both axis extrema, including nonzero azimuth.
    for polar in (low,high):
        for phi in (-pi/2-phase,-pi/2+phase):
            radial=np.array([sin(polar)*cos(phi),sin(polar)*sin(phi),cos(polar)])
            polar_vec=np.array([cos(polar)*cos(phi),cos(polar)*sin(phi),-sin(polar)])
            phase_vec=np.array([-sin(phi),cos(phi),0])
            for off in np.linspace(0,2*pi,121):
                n=radial*cos(cap)+(polar_vec*cos(off)+phase_vec*sin(off))*sin(cap)
                ratio=tan(acos(n[2])/2)
                assert mod_range[0]-1e-13<=ratio<=mod_range[1]+1e-13
                assert abs(atan2(n[1],n[0])+pi/2)<=azimuth+1e-13
    return {"worst_endpoint_pass_probability":float((ref.conj()@((I+a)/2)@ref).real),
            "signed_quarter_X_expectation":qx, "modulus":abs(r),
            "equal_damping_purity":pur,"equal_damping_apparent_modulus":float(sqrt(mixed[1,1].real/mixed[0,0].real)),
            "tolerance_modulus":mod_range,"tolerance_phase_error_degrees":degrees(azimuth),
            "closed_loop_phase_error_degrees":5*degrees(azimuth),
            "anchored_H13_normalized_modulus_ratio":float(abs(hu[0,2]/hd[0,2])),
            "anchored_loop_imag_normalized":float(loop.imag),
            "HS_rounded_angles_alpha_beta_gamma_degrees":angles,"HS_rounded_J":J}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--write-new",action="store_true")
    mode.add_argument("--check",action="store_true")
    args=parser.parse_args()
    out=run()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf-8") as stream:
            json.dump(out,stream,indent=2)
            stream.write("\n")
    elif args.check:
        saved=json.loads((HERE/"results.json").read_text(encoding="utf-8"))
        assert saved.keys()==out.keys()
        for k in out:
            assert np.allclose(out[k],saved[k],rtol=3e-10,atol=1e-15)
        print("1444 signed comparison, mixed competitor, finite tolerance and physical loop verified")
    else:
        print(json.dumps(out,indent=2))
