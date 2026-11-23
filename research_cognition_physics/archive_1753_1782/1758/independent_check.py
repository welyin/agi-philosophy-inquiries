"""Can the same unbroken-phase pseudo-Dirac pair splittings calibrate the Majorana arm?"""
import argparse
import json
import math
import runpy
from pathlib import Path
from types import FunctionType, SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
VEV = 246.22/math.sqrt(2)
SOLAR = 7.54e-5


def read(path):
    return json.loads(path.read_bytes())


def poles(D, mu):
    z = np.zeros((2,2),complex)
    M = np.block([[z,D],[D.T,mu]])
    return np.linalg.svd(M,compute_uv=False)[::-1]


def old_pairs():
    package = runpy.run_path(str(ARCH/"1698"/"check.py"))
    conf = read(ARCH/"1698"/"inputs.json")
    _,c,h = package["neutral_graph"](conf)
    point = read(ARCH/"1698"/"results.json")["shared_witness"]
    k = conf["neutral_shape_k"]
    b0 = point["neutral_b"]*k*k
    c,h = b0*k*c,b0*k**4*h
    # First-order pair splitting in the Dirac singular basis is |(V.T h V)_ii|.
    _,M,vh = np.linalg.svd(c)
    V = vh.conj().T[:,::-1]
    diagonal = abs(np.diag(V.T@h@V))
    p = poles(c,h)
    exact = np.diff(p)[::2]
    np.testing.assert_allclose(exact,diagonal,atol=8e-12,rtol=0)
    return M[::-1],diagonal,p


def physical(s,u):
    joint = runpy.run_path(str(ARCH/"1755"/"joint.py"))
    captured={}
    def loader(path):
        pkg = runpy.run_path(str(path))
        if Path(path).parent.name=="1748":
            initial=pkg["daughter"]
            def daughter(bg,F):
                a,r=initial(bg,F)
                captured.update(bg=bg,filter=F,charged=a)
                return a,r
            pkg["daughter"]=daughter
        return pkg
    p=joint["point"]
    point=FunctionType(p.__code__,dict(p.__globals__,runpy=SimpleNamespace(run_path=loader)),
                       p.__name__,p.__defaults__)
    row=point(ARCH/"1755",s,u,.04216,.125)
    bg,F,a=captured["bg"],captured["filter"],captured["charged"]
    flow,pt=bg["pt"]["flow"],bg["pt"]
    y=bg["yW"].copy()
    y[2]=a*F@y[2]
    y[3]=F@y[3]
    T=np.eye(3,dtype=complex)
    threshold=flow["fields"]["neutral_poles_GeV"]
    external={}
    for start,susy,dt,g in reversed(bg["segments"]):
        y,T=bg["model"]["advance"](y,T,-dt,g[::-1],flow)
        if start==bg["Mcolor"]:
            U,m,Vh=np.linalg.svd(y[1,:,:3])
            c=abs(np.vdot(np.linalg.svd(y[0,:,:3])[0][:,-1],U[:,-1]))
            y[1,:,:3]=U@np.diag(m*[math.sqrt((1+c)/2),1,c])@Vh
        for j,scale in enumerate(threshold):
            if start==scale:
                if j in (0,2): 
                    external[j]=y[3,:,j].copy()
                y[3,:,j]=0
    L=flow["extraction"](y)[3]
    phases=np.diag(L.conj().T@pt["low_E"])
    U=np.linalg.solve(T,L*(phases/abs(phases)))
    carrier=np.linalg.inv(bg["S"]@pt["low_E"])
    W=carrier@F@U
    C=W.T@pt["C"]@W
    b=W.T@pt["b"]
    vec=np.outer(b,b)
    M=np.array(threshold)[[0,2]]
    theta=(F@bg["N"]).conj()*VEV/M
    low=U.T@theta
    inverse=np.linalg.pinv(low)
    kernels=[inverse@v@inverse.T for v in (C,vec)]
    for v,m in zip((C,vec),kernels):
        np.testing.assert_allclose(low@m@low.T,v,atol=3e-12*np.linalg.norm(v))
    rho=bg["nu"]["FROZEN"][0]
    obs=observables(low,(kernels[0]+rho*kernels[1]))
    assert abs(obs["reactor"]-row["angle_readings"][0])<1e-7
    y2=np.array([np.vdot(external[j],external[j]).real for j in (0,2)])
    return row, M,low,kernels,y2


def observables(low,mu):
    mass=low@mu@low.T
    left,m,_=np.linalg.svd(mass)
    p=abs(left[0,0])**2
    norm=math.sqrt(SOLAR)/m[1]
    return dict(reactor=4*p*(1-p),ratio=m[0]/m[1],solar_gain=norm,
                mee_meV=1000*abs(mass[0,0])*norm,solar_angle=abs(left[0,1])**2/(1-p))


def pair_quadratic(a,b):
    return np.array([[abs(a)**2,(a.conjugate()*b).real],[(a.conjugate()*b).real,abs(b)**2]])


def joint_range(A,B,low,q,M):
    from numpy.polynomial import Polynomial as Poly
    lowA=low@A@low.T
    lowB=low@B@low.T
    i=np.argmax(abs(np.diag(lowB)))
    arm=lowB[:,i]/np.sqrt(lowB[i,i])
    np.testing.assert_allclose(np.outer(arm,arm),lowB,atol=np.linalg.norm(lowB)*5e-12)
    previous=runpy.run_path(str(ARCH/"1757"/"independent_check.py"))
    gain, forward,clock,reactor=previous["model"](lowA,arm)
    kam,sk=previous["KAM"],previous["SK"]
    lower=math.sqrt(1+(sk[0]-sk[1])/(kam[0]+kam[2]))
    upper=math.sqrt(1+(sk[0]+sk[2])/(kam[0]-kam[1]))
    endpoints=sorted(r["rho"] for rr in (clock(lower),clock(upper),reactor(.0803),reactor(.0899)) for r in rr)
    regions=[]
    for left,right in zip(endpoints,endpoints[1:]):
        obs=forward((left+right)/(2*gain))
        if lower<obs["R"]<upper and .0803<obs["reactor"]<.0899:
            regions.append([left,right])
    def quad(x):
        return Poly([x[0,0],2*x[0,1],x[1,1]])
    num,den=map(quad,q)
    extrema=(num.deriv()*den-num*den.deriv()).trim(1e-14*max(num.coef)*max(den.coef)).roots()
    values=[]
    for left,right in regions:
        spots=[left,right]+[t.real for t in extrema if abs(t.imag)<1e-10 and left<t.real<right]
        values.extend(math.sqrt(num(x)/den(x)) for x in spots)
    # Positive Gram of the Frobenius norm bounds the |22| diagonal from below,
    # uniformly over the real projective line, independently of common matching.
    frob=sum(pair_quadratic(a,b) for a,b in zip(A.flat,B.flat))
    inv=np.linalg.inv(np.linalg.cholesky(frob))
    mass_floor=math.sqrt(np.linalg.eigvalsh(inv@q[1]@inv.T)[0])
    spec_gap=min(np.diff(M)[0],2*M[0])
    cap=.0001  # GeV; explicit sterile LNV norm budget.
    rel=2*cap/(spec_gap-2*cap)/mass_floor
    assert rel<1e-4
    # Bound on all four bare singular values, via Hermitian dilation.
    r_lo,r_hi=min(values),max(values)
    return dict(relative_rho_regions=regions,leading_ratio=[r_lo,r_hi],
                lower_bound_abs_mu22_over_Frobenius=mass_floor,
                finite_100_keV_max_abs_error_per_split_over_abs_mu22=rel,
                exact_pair_ratio_outer_at_100_keV=[(r_lo-rel)/(1+rel),(r_hi+rel)/(1-rel)],
                max_mu_norm_keV_from_KamLAND_SK_shape_units=1e-3*math.sqrt(kam[0]+kam[2])*upper
                    /np.linalg.svd(low,compute_uv=False)[1]**2)


def run_case(s,u,rho0,old_split):
    row,M,low,bb,y2=physical(s,u)
    A,B=bb[0],rho0*bb[1]
    q=[pair_quadratic(A[j,j],B[j,j]) for j in (0,1)]
    inv=np.linalg.inv(np.linalg.cholesky(q[1]))
    ratio_range=np.sqrt(np.linalg.eigvalsh(inv@q[0]@inv.T))
    target=old_split[0]/old_split[1]
    delta=q[0]-target**2*q[1]
    rr=np.roots([delta[1,1],2*delta[1,0],delta[0,0]])
    assert np.min(abs(rr.imag))>1e-9 or min(abs(rr[0]-rr[1]),1)>1e-8
    G=M*y2/(16*math.pi)
    vector_phase=lambda mass:(1-(mass/M)**2)**2*(1+2*(mass/M)**2)
    physical_width=G*(2*vector_phase(80.379)+vector_phase(91.1876)+(1-(125.25/M)**2)**2)/4
    def entry(x):
        mu=A+x*B
        orig=observables(low,mu)
        K=orig["solar_gain"]*1e-9*mu
        masses=poles(np.diag(M),K)
        split=np.diff(masses)[::2]
        diag=abs(np.diag(K))
        S=np.diag([1.,-1.])
        rival=S@K@S
        pe=poles(np.diag(M),rival)
        np.testing.assert_allclose(pe,masses,atol=1e-11,rtol=0)
        reflected=observables(low,S@mu@S)
        # Same fixed physical absolute source normalization in the rival.
        reflected["mee_same_source_scale_meV"]=reflected["mee_meV"]*orig["solar_gain"]/reflected["solar_gain"]
        reflected["dm21_same_source_scale"]=SOLAR*(orig["solar_gain"]/reflected["solar_gain"])**2
        Rss=split**2/(2*physical_width**2+split**2)
        return dict(rho=rho0*x,light=orig,matrix_norm_keV=np.linalg.norm(K,2)*1e6,
                    diagonal_split_keV=diag*1e6,exact_bare_split_keV=split*1e6,
                    correction_keV=(split-diag)*1e6,
                    relative_to_old_two_splittings=diag/old_split,
                    rate_ratio_SS_OS_isolated=Rss,opposite_interference=reflected,
                    unchanged_four_bare_poles_max_GeV=max(abs(pe-masses)))
    roots=[entry(x.real) for x in rr if abs(x.imag)<1e-9]
    for p in roots:
        assert abs(p["diagonal_split_keV"][0]/p["diagonal_split_keV"][1]/target-1)<3e-12
    joint=joint_range(A,B,low,q,M)
    upper=joint["exact_pair_ratio_outer_at_100_keV"][1]
    joint["necessary_relative_tolerance_each_old_split"]=(target-upper)/(target+upper)
    assert joint["necessary_relative_tolerance_each_old_split"]>.07
    return dict(Vus=row["trial_us"],Vub=u,Vcb=row["low_moduli"]["Vcb"],M_GeV=M,weak_leg_y2=y2,
                dirac_widths_massless_GeV=G,dirac_widths_finite_WZh_GeV=physical_width,first_order_ratio_full_real=ratio_range,
                joint_free_clock_and_reactor=joint,
                old_splitting_ratio_real_roots=roots,working_gain=entry(1.),
                singular_values_low_port=np.linalg.svd(low,compute_uv=False))


def clean(x):
    if isinstance(x,dict): return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)): return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)): return float(f"{x:.10g}")
    return x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    _,d,p=old_pairs()
    old=read(ARCH/"1757"/"independent_results.json")
    point=old["combined_central_solution"]
    r=old["rows"][2]
    cases=[run_case(point["Vus"],point["Vub"],point["rho"],d),
           run_case(r["Vus"],r["Vub"],r["reactor_central_real_roots"][1]["rho"],d)]
    out=clean(dict(original_canonical_dirac_splittings_keV=d*1e6,
                   old_ratio=d[0]/d[1],original_bare_svd_poles_GeV=p,cases=cases))
    dest=HERE/"independent_results.json"
    if args.check:
        assert json.loads(dest.read_bytes())==out
        print("1758 original sterile splittings, inverse mass kernel, coherence and reflected phase verified")
    else:
        dest.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n",encoding="utf8")
