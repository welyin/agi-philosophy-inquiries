"""Inverse reactor-frequency test over both source branches, at the two frozen charged bases."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE=Path(__file__).resolve().parent
ARCH=HERE.parent
CKM=runpy.run_path(str(ARCH/"1678/independent_check.py"))
IV=runpy.run_path(str(ARCH/"1675/phase_box_check.py"))
add,mul,lo,hi=[IV[k] for k in ("add","prod","down","up")]


def read(name):
    return json.loads((HERE/name).read_bytes())


def sub(a,b):
    return add(a,(-b[1],-b[0]))


def div(a,b):
    assert np.all(b[0]>0)
    return mul(a,(lo(1/b[1]),hi(1/b[0])))


def unpack(p):
    return round(float(p),7)


def inputs(w):
    cfg=read("../1397/inputs.json")["high_scale"]
    c,d=np.array(cfg["down"]+cfg["angles_theta12_theta13_theta23_delta"][:3]).T
    short,long,_,_=CKM["get_matrices"](c+d*np.array(w["quark_marginal_pulls"]),w["gamma_matching"])
    phase=-math.pi/6
    target=(short+w["mag_t"]*np.exp(1j*phase)*long).T
    basis,singular,_=np.linalg.svd(target)
    e=basis[:,-1]
    residual=np.linalg.norm(target@target.conj().T@e-singular[-1]**2*e)
    assert residual/(singular[-2]**2-singular[-1]**2)<1e-10
    h=w["neutral_modulus"]
    unit=np.exp(1j*np.r_[0.,2*phase,3*phase+np.angle(3+2*np.exp(4j*phase))])
    a=unit*np.array([0,1,h])
    b=unit*np.array([1,3,h])
    null=unit.conj()*np.array([2,-1,1/h])
    null/=np.linalg.norm(null)
    pspace=1-abs(np.vdot(null,e))**2
    A,B,C=1+h*h,10+h*h,3+h*h
    D=A*B-C*C
    eta=(-1+2j*np.sqrt(3))/np.sqrt(13)
    va,vb=a@e,b@e
    baseline=np.outer(a,a)+w["source_r"]*np.exp(2j*np.pi/3)*np.outer(b,b)
    m2,modes=np.linalg.eigh(baseline.conj().T@baseline)
    pe=abs(np.vdot(e,modes[:,-1]))**2
    np.testing.assert_allclose([math.sqrt(m2[-1]/m2[-2]),4*pe*(1-pe)],
                               [w["R"],w["three_angles"][0]],atol=2e-10)
    quad=[A*A,B*B,2*C*C*eta.real,D,A*abs(va)**2,B*abs(vb)**2,
          2*C*(eta*va.conj()*vb).real,pspace]
    # Compare the scalar projector formula to direct diagonalizations at very different strengths.
    for y in [.01,.07,.15,.24,.84]:
        m=(1-y)*np.outer(a,a)+y*eta*np.outer(b,b)
        vals,vecs=np.linalg.eigh(m.conj().T@m)
        T=quad[0]*(1-y)**2+quad[1]*y*y+quad[2]*y*(1-y)
        Cmass=D*y*(1-y)
        S=quad[4]*(1-y)**2+quad[5]*y*y+quad[6]*y*(1-y)
        diff=math.sqrt(T*T-4*Cmass*Cmass)
        p=pspace/2+(2*S-pspace*T)/(2*diff)
        np.testing.assert_allclose([p,(T+diff)/(2*Cmass)],
                                   [abs(np.vdot(e,vecs[:,-1]))**2,math.sqrt(vals[-1]/vals[1])],rtol=1e-10)
    return quad


def envelope(w):
    constants=inputs(w)
    centers=np.array([unpack(v) for v in constants])
    pad=1e-6
    # The certificate extends to the entire independent coefficient box of this radius.
    c=[(lo(q-pad),hi(q+pad)) for q in centers]
    count=2**17
    grid=np.arange(count+1)/count
    y=(grid[:-1],grid[1:])
    one=(1.,1.)
    u=sub(one,y)
    xy=mul(y,u)
    y2,u2=mul(y,y),mul(u,u)
    T=add(add(mul(c[0],u2),mul(c[1],y2)),mul(c[2],xy))
    S=add(add(mul(c[4],u2),mul(c[5],y2)),mul(c[6],xy))
    massprod=mul(c[3],xy)
    discr=sub(mul(T,T),mul((4.,4.),mul(massprod,massprod)))
    assert np.all(discr[0]>0)
    diff=(lo(np.sqrt(discr[0])),hi(np.sqrt(discr[1])))
    p=mul((.5,.5),add(c[7],div(sub(mul((2.,2.),S),mul(c[7],T)),diff)))
    p=(np.maximum(0,p[0]),np.minimum(1,p[1]))
    signal=mul((4.,4.),mul(p,sub(one,p)))
    data=read("../1659/results.json")["joint_windows"][0]
    allowed=(signal[1]>=lo(data[0]))&(signal[0]<=hi(data[1]))
    ix=np.where(allowed)[0]
    assert len(ix)>0 and ix[0]>0 and ix[-1]<count-1
    assert np.all(np.diff(ix)==1)
    pair=add(T,diff)
    def rr(mask):
        numerator=tuple(z[mask] for z in pair)
        denominator=tuple(z[mask] for z in mul((2.,2.),massprod))
        return div(numerator,denominator)
    R=rr(allowed)
    permitted_R=[float(min(R[0])),float(max(R[1]))]
    h=w["neutral_modulus"]
    r_crit=(1+h*h)/(10+h*h)
    assert y[1][ix[-1]]/(1-y[1][ix[-1]])<r_crit
    gap=CKM["NEUTRINO"]["gap_box"](2)
    lower=permitted_R[0]
    assert lower>gap[1]
    # Conversely cover any frequencies in the double marginal data box, including endpoints.
    # Avoid division by massprod in endpoint cells by cross-multiplying the inequalities.
    two=mul((2.,2.),massprod)
    below=mul((lo(gap[0]),lo(gap[0])),two)
    above=mul((hi(gap[1]),hi(gap[1])),two)
    juno=(pair[1]>=below[0])&(pair[0]<=above[1])
    sep=r_crit/(1+r_crit)
    near=juno&(y[0]<sep)
    far=juno&(y[1]>sep)
    branch_rates=[[float(signal[0][mask].min()),float(signal[1][mask].max())] for mask in (near,far)]
    assert min(v[0] for v in branch_rates)>data[1]
    return {"coeff_centers_t1_t2_t3_D_s1_s2_s3_subspace":centers.tolist(),"coeff_halfwidth":pad,
            "reactor_only_y_outer":[float(y[0][ix[0]]),float(y[1][ix[-1]])],
            "reactor_only_R_outer":permitted_R,"reactor_only_large_source_branch_empty":True,
            "at_least_fractional_increase_of_R2_above_JUNO_two_width_upper":(lower/gap[1])**2-1,
            "reactor_signal_outer_given_JUNO_two_widths_small_then_large":branch_rates}


def certificate():
    return {"source_phase_deg":math.degrees(np.angle(-1+2j*np.sqrt(3))),
            "JUNO_two_width_R":CKM["NEUTRINO"]["gap_box"](2).tolist(),
            "source_cells":2**17,
            "ordered_1682_witnesses":[envelope(w) for w in read("../1682/independent_results.json")["witnesses"]]}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    actual=certificate()
    if args.check:
        assert read("independent_results.json")==actual
        print("PASS: inverse reactor and JUNO enclosures over all nonnegative source strengths")
    else:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as stream:
            stream.write(json.dumps(actual,indent=2)+"\n")
