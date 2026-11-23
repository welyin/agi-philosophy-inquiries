"""Tree unitarity triangle, independent time-ordered holdouts and filter diagnostics."""
import argparse
import cmath
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def beta_sine(radius,gamma):
    angle=math.atan2(radius*math.sin(gamma), 1-radius*math.cos(gamma))
    return angle,math.sin(2*angle)


def tree_triangle(p):
    prod,form=p["K_product"],p["K_form"]
    k_lo=(prod[0]-prod[1])/(form[0]+form[1])
    k_hi=(prod[0]+prod[1])/(form[0]-form[1])
    rho=p["Vub_over_Vcb"]
    rlo,rhi=rho[0]-sum(rho[1:]),rho[0]+sum(rho[1:])
    w=p["Vcb_max"]
    c13_lo=math.sqrt(1-(rhi*w)**2)
    s12_hi=k_hi/c13_lo
    c12_lo,c12_hi=math.sqrt(1-s12_hi**2),math.sqrt(1-k_lo**2)
    s23_hi=w/c13_lo
    c23_lo=math.sqrt(1-s23_hi**2)
    extra=c12_hi*s23_hi*rhi*w
    vcd_lo,vcd_hi=k_lo*c23_lo-extra,s12_hi+extra
    rb=[rlo*c12_lo*c13_lo/vcd_hi,rhi*c12_hi/vcd_lo]
    arcs={}
    for offset in [0,180]:
        lo,hi=[math.radians(g+offset) for g in p["gamma_working_arc_deg"]]
        candidates=[]
        for radius in rb:
            g_values=[lo,hi]
            for g in [math.acos(radius), 2*math.pi-math.acos(radius)]:
                if lo<=g<=hi:
                    g_values.append(g)
            candidates.extend(beta_sine(radius,g)[0] for g in g_values)
        bounds=[min(candidates),max(candidates)]
        assert max(map(abs,bounds))<math.pi/4
        ss=[math.sin(2*x) for x in bounds]
        # Construct an exact PDG matrix point for the center using the tree quantities.
        gamma=math.radians(p["gamma_central_deg"]+offset)
        k,vc=prod[0]/form[0],p["central_Vcb"]
        s13=rho[0]*vc
        c13=math.sqrt(1-s13*s13)
        s12,s23=k/c13,vc/c13
        c12,c23=math.sqrt(1-s12*s12),math.sqrt(1-s23*s23)
        delta=gamma+math.asin(c12*s23*s13/(s12*c23)*math.sin(gamma))
        cd=-s12*c23-c12*s23*s13*cmath.exp(1j*delta)
        radius=rho[0]*c12*c13/abs(cd)
        beta, exact=beta_sine(radius,gamma)
        td=s12*s23-c12*c23*s13*cmath.exp(1j*delta)
        tb=c23*c13
        assert abs(cmath.phase(-cd*vc/(td*tb))-beta)<1e-12
        arcs[str(offset)]=dict(beta_deg=[math.degrees(v) for v in bounds],
            ordinary_sine=ss,DSplus=[-2*ss[1],-2*ss[0]],DSminus=[2*v for v in ss],
            center_DSplus=-2*exact, center_DSminus=2*exact)
    return rb,arcs


def ray_density(v):
    v=np.asarray(v,dtype=complex)
    v=v/np.linalg.norm(v)
    return np.outer(v,v.conj())


def filter_tests():
    rho=ray_density([1,.3+.4j])
    mixed=.67*ray_density([1,2j])+.33*ray_density([2,-1j])
    sigma_x=np.array([[0.,1.],[1.,0.]])
    def transfer(h):
        return math.cos(h/2)*np.eye(2)-1j*math.sin(h/2)*sigma_x
    def prob(a,b,u):
        return np.trace(b @ u @ a @ u.conj().T).real
    identity=np.eye(2)
    w=1-np.trace(rho@mixed).real
    trace_dist=.5*np.sum(np.abs(np.linalg.eigvalsh(rho-(identity-mixed))))
    assert trace_dist**2<=1-w+1e-13
    for phase in np.linspace(0,math.pi/3,101):
        u=transfer(phase)
        ordered=(prob(identity-rho,mixed,u)+prob(identity-mixed,rho,u))/2
        assert abs(ordered-w)<=math.sin(phase/2)**2+1e-13
    u=transfer(math.pi/2)
    z=ray_density([1,0])
    complex_ray=ray_density([1,1j])
    f,r=prob(z,complex_ray,u),prob(complex_ray,z,u)
    assert np.isclose(f,0) and np.isclose(r,1)
    return dict(complex_filter_reciprocal_probabilities=[round(f,12),round(r,12)],
                mixed_filter_equal_time=w,mixed_filter_trace_distance=trace_dist)


def run(d):
    rb,arcs=tree_triangle(d["tree"])
    boxes={k:[v[0]-2*sum(v[1:]),v[0]+2*sum(v[1:])] for k,v in d["babar"].items()}
    worst_freq=sum(d["split_ps_inverse"])
    time_leak=math.sin(worst_freq*d["window_ps"]/2)**2
    margins=[abs(boxes["DS_plus"][1]),boxes["DS_minus"][0]]
    assert arcs["0"]["DSplus"][0]<=boxes["DS_plus"][1]
    assert arcs["180"]["DSplus"][0]>0
    return dict(rb_outer=rb,branches=arcs,experimental_rectangles=boxes,
        minimum_mirror_coefficient_shift=arcs["180"]["DSplus"][0]-boxes["DS_plus"][1],
        finite_clock_bias=time_leak, required_Wlower_if_other_errors_zero=[
            1-(m/2)**2+time_leak for m in margins],
        example_filter_distance_ceiling=math.sqrt(1-d["example_lower_normalized_double_tag"]+time_leak),
        filter_tests=filter_tests())


def close(a,b):
    if isinstance(a,dict): return a.keys()==b.keys() and all(close(a[k],b[k]) for k in a)
    if isinstance(a,list): return len(a)==len(b) and all(close(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-13)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check",action="store_true")
    mode.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    answer=run(json.loads((HERE/"inputs.json").read_text(encoding="utf-8-sig")))
    target=HERE/"results.json"
    if args.write_new:
        with target.open("x",encoding="utf8") as f:
            json.dump(answer,f,indent=2)
            f.write("\n")
    else:
        assert close(answer,json.loads(target.read_text(encoding="utf8")))
        print("1427 tree and two-sided filtering verified")
