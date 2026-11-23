"""Lower holomorphic contacts, the neutron EDM and the first determinant vertex."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np
from numpy.polynomial import Polynomial as Poly


HERE=Path(__file__).resolve().parent
OMEGA=complex(-.5,-math.sqrt(3)/2)


def read(name):
    return json.loads((HERE/name).read_bytes())


def compensated(a):
    # Minimize |L| at fixed a and an exactly positive determinant.
    row=np.array([-a,1.])
    gram=np.array([[1.,-.5],[-.5,1.]])
    inv_row=np.linalg.solve(gram,row)
    m,n=inv_row*(a**3/(a*a-a+1)/(row@inv_row))
    b=(m+n)/(1+m)-a
    c=a+b-a*b
    return b,c


def determinant(a,b,c):
    t=Poly([0.,1.])
    p=(t+a)*(t*t+b*t+c)
    v=p(OMEGA)
    vertex=OMEGA*p.deriv()(OMEGA)/v-3
    formula=1+a*c-.5*(a+b+c+a*b)+1j*math.sqrt(3)/2*(a+b-c-a*b)
    assert abs(v-formula)<5e-15
    return np.angle(v),vertex


def width(a,b,near_limit):
    real=(1-b)*(a*a-a+1)
    tilting=a-.5
    tangent=math.tan(near_limit)
    return [-real*tangent/(math.sqrt(3)/2+tilting*tangent),
             real*tangent/(math.sqrt(3)/2-tilting*tangent)]


def mixing(background,ratio,a,charm):
    n,Q,W,A,arm=background
    phases=np.array([(OMEGA+a)/abs(OMEGA+a),charm/abs(charm),1.])
    C=W.T@np.diag(1/(A*phases))@W
    M=Q.T@(C-np.outer(C@n,C@n)/(n@C@n))@Q
    v=Q.T@arm
    T=np.outer(v,v)
    xnorm=np.linalg.norm(M,2)/np.linalg.norm(T,2)
    T=T*xnorm
    # The determinant of M+xT is affine since the new arm has rank one.
    det=Poly([np.linalg.det(M), M[0,0]*T[1,1]+M[1,1]*T[0,0]-2*M[0,1]*T[0,1]])
    trace=Poly([np.vdot(M,M).real,2*np.vdot(M,T).real,np.vdot(T,T).real])
    eq=trace*trace-(ratio+1/ratio)**2*det*Poly(det.coef.conj())
    roots=Poly(eq.coef.real/max(abs(eq.coef))).roots()
    rows=[]
    for root in roots:
        if abs(root.imag)>1e-7:
            continue
        x=root.real
        u,s,_=np.linalg.svd(M+x*T)
        assert abs(s[0]/s[1]-ratio)<5e-10
        mix=np.column_stack([n,Q@u[:,::-1].conj()])
        wt=abs(mix)**2
        z=wt[0,2]
        sol,atm=wt[0,1]/(1-z),wt[1,2]/(1-z)
        j=np.imag(mix[0,0]*mix[1,1]*mix[0,1].conj()*mix[1,0].conj())
        cos=(wt[1,0]-sol*(1-atm)-(1-sol)*atm*z)/(2*math.sqrt(sol*(1-sol)*atm*(1-atm)*z))
        sin=j/((1-z)*math.sqrt(sol*(1-sol)*atm*(1-atm)*z))
        delta=np.arctan2(sin,cos)
        win=np.array(read("../1659/results.json")["joint_windows"])
        angles=np.array([4*z*(1-z),sol,atm])
        contours=read("../1705/independent_results.json")["official_coordinates_rad"]
        inside=runpy.run_path(str(HERE/"../1705/independent_check.py"))["interior"]
        fc=any(inside(contours[key],atm,delta) for key in ("conf90_0","conf90_1"))
        rows.append({"raw_real_coefficient":x*xnorm,"angles":angles.tolist(),
                     "delta_deg":np.degrees(delta),"within_three_windows":bool(np.all((angles>=win[:,0])&(angles<=win[:,1]))),
                     "inside_FC90":fc})
    # Compare the observed-frequency branch with the lowest reactor double angle.
    assert len(rows)==2
    return min(rows,key=lambda r:r["angles"][0])


def compute():
    ed=read("../1389/results.json")["conditional_edm"]
    eps=ed["theta_limit_rad"]
    recorded=read("../1718/results.json")
    nu=read("../1668/inputs.json")
    ratio=math.sqrt(nu["juno_NO_dm31"][0]/nu["juno_NO_dm21"][0])
    sim=runpy.run_path(str(HERE/"../1716/check.py"))
    n,Q,_,_,_,(W,A,arms,_),_=sim["background"](.1)
    axes=(n,Q,W,A,arms[:,1])
    out={"single_source_nEDM_angle_budget":eps,"delta_c_interval_at_origin":width(0,0,eps)}
    dimless=max(abs(np.array(width(0,0,eps))))
    out["one_at_a_time_lower_Wilson_leading_bounds"]={"up_identity":A[0]*dimless,
        "charm_I":A[1]*32*dimless,"charm_identity":A[1]*dimless}
    for tag,(a,b,c) in {"pure":(0.,0.,0.),"isolated_up_direct":(.01,0.,0.),
                         "EDM_balanced_best_vertex":(.05,*compensated(.05))}.items():
        th,vertex=determinant(a,b,c)
        obs=mixing(axes,ratio,a,OMEGA**2+b*OMEGA+c)
        out[tag]={"a_b_c":[a,b,c],"determinant_phase_rad":th,
                  "excess_complex_dlogI_vertex":[vertex.real,vertex.imag],
                  "excess_vertex_modulus":abs(vertex),"phase_accepted_small_branch":abs(th)<=eps,
                  "light":obs,"allowed_transverse_c_error_if_exactly_balanced":width(a,b,eps)}
        assert obs["within_three_windows"] and obs["inside_FC90"]
    pure=min(recorded["original_basis_verification"]["direct_forward_both_roots"],key=lambda r:r["angles"][0])
    assert abs(out["pure"]["light"]["delta_deg"]-pure["delta_deg"])<2e-7
    for a in [-.18,-.11,.007,.05,.1,.25]:
        b,c=compensated(a)
        th,vertex=determinant(a,b,c)
        sharp=1.5*abs(a)**3/(a*a-a+1)**1.5
        assert abs(vertex)-sharp<1e-12 and abs(vertex)>sharp-1e-12 and abs(th)<1e-12
    a=.05
    t=Poly([0,1])
    E=t*t+t+1
    gain=a**3/(a*a-a+1)**2
    H=t**3+gain*E**2
    charm,rem=divmod(H,t+a)
    assert max(abs(rem.coef))<1e-14 and abs(H(OMEGA)-1)<1e-14
    assert abs(OMEGA*H.deriv()(OMEGA)/H(OMEGA)-3)<3e-14
    q=mixing(axes,ratio,a,charm(OMEGA))
    np.testing.assert_allclose(q["angles"],out["EDM_balanced_best_vertex"]["light"]["angles"],atol=1e-12)
    out["fourth_degree_same_vertex_competitor"]={"a":a,"gain":gain,"charm_coefficients_ascending":charm.coef.tolist(),
        "second_derivative_excess_modulus":abs(H.deriv(2)(OMEGA)-6*OMEGA),"light":q}
    return out


def tidy(x):
    if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
    if isinstance(x,list): return [tidy(v) for v in x]
    if isinstance(x,(bool,np.bool_)): return bool(x)
    return float(f"{x:.11g}")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    table=tidy(compute())
    if args.check:
        assert table==read("independent_results.json")
        print("1722 independent EDM surface, determinant derivative and PMNS witnesses verified")
    else:
        (HERE/"independent_results.json").write_text(json.dumps(table,indent=2)+"\n",encoding="utf-8")

