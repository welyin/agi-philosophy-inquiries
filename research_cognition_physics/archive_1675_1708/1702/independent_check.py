"""Charged flavor transport in the quark mass frame, with phase-free semileptonic calibration."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE = Path(__file__).resolve().parent


def read(file):
    return json.loads((HERE/file).read_bytes())


def old(file):
    return runpy.run_path(str(HERE/file))


END = old("../1699/check.py")
CHARGED = old("../1678/independent_check.py")
PREV = CHARGED["CKM"]
DATA = read("../1386/inputs.json")
SPECTRUM = END["spec"]
TAN = read("../1699/results.json")["ordinary_finite_matching"]["minimum_root_high"]
LOOP = 16*math.pi**2


def ckm(sines, gamma):
    a,b,c = np.arcsin(sines)
    g = math.radians(gamma)
    delta = g+math.asin(math.cos(a)*math.tan(c)*math.sin(b)/math.sin(a)*math.sin(g))
    return PREV["texture"]([.001,.02,1,a,b,c],delta)[1]


def generator(values, partner):
    diff = values[None,:]-values[:,None]
    divided = np.divide(values[None,:]+values[:,None],diff,
                        out=np.zeros_like(diff),where=~np.eye(3,dtype=bool))
    return divided*partner


def run(sines,gamma,frame,tan,step=.15):
    """Integrate all six quark singular values and the unitary CKM as separate coordinates."""
    ews = DATA["yukawa"]
    u,d,e = (np.array(ews[name])[:,0] for name in ("up","down","charged_leptons"))
    y = np.r_[u,d,ckm(sines,gamma).ravel(),(frame*e).ravel(),np.zeros(12),0].astype(complex)
    g2 = np.array([END["weak"][k][0]**2 for k in ("g1_GUT","g2","g3")])
    soft = END["soft"]
    masses = SPECTRUM["neutral_poles_GeV"]
    weights = [m for m,_ in END["thresholds"]]+masses
    edges = sorted(set([END["weak"]["mu_ref_GeV"],soft,END["cutoff"]]+weights))
    for start,end in zip(edges,edges[1:]):
        if start == soft:
            c = 1/math.hypot(1,tan)
            s = c*tan
            y[0:3] /= s
            y[3:6] /= c
            y[15:24] /= c
            y[24:36] /= s
            g2 = 1/(1/g2 - np.array([0,2,3])/(48*math.pi**2))
        for j,m in enumerate(masses):
            if start == m:
                y[24:36].reshape(3,4)[:,j] = END["Y0"][:,j]
        susy = start >= soft
        slope = np.array([33/5,1,-3.] if susy else [4.1,-19/6,-7.])
        for mass,coeff in END["thresholds"]:
            if start >= mass:
                slope+=coeff
        def rhs(time,z):
            g1,w,strong = 1/(1/g2-2*slope*time)
            u,d = z[:3].real,z[3:6].real
            q,r = u*u,d*d
            v = z[6:15].reshape(3,3)
            le = z[15:24].reshape(3,3)
            n = z[24:36].reshape(3,4)
            p,nu = le@le.conj().T,n@n.conj().T
            tle,tn = np.trace(p).real,np.trace(nu).real
            tr = 3*(sum(q)+sum(r))+tle+tn
            if susy:
                a,b = 3.,1.
                gu = 3*sum(q)+tn-13*g1/15-3*w-16*strong/3
                gd = 3*sum(r)+tle-7*g1/15-3*w-16*strong/3
                ge,gn = 3*sum(r)+tle-9*g1/5-3*w,3*sum(q)+tn-3*g1/5-3*w
            else:
                a,b = 1.5,-1.5
                gu,gd = tr-17*g1/20-9*w/4-8*strong,tr-g1/4-9*w/4-8*strong
                ge,gn = tr-9*g1/4-9*w/4,tr-9*g1/20-9*w/4
            rd = (v*r)@v.conj().T
            qu = (v.conj().T*q)@v
            dv = b*(-generator(q,rd)@v+v@generator(r,qu))
            return np.r_[u*(gu+a*q+b*np.diag(rd).real),d*(gd+a*r+b*np.diag(qu).real),
                         dv.ravel(),((a*p+b*nu)@le+ge*le).ravel(),
                         ((a*nu+b*p)@n+gn*n).ravel(),2*a*(sum(q)+sum(r))]
        width = math.log(end/start)
        steps = max(1,math.ceil(width/step))
        dt = width/(LOOP*steps)
        for i in range(steps):
            t = i*dt
            z1 = rhs(t,y)
            z2 = rhs(t+dt/2,y+dt*z1/2)
            z3 = rhs(t+dt/2,y+dt*z2/2)
            z4 = rhs(t+dt,y+dt*z3)
            y += dt*(z1+2*z2+2*z3+z4)/6
        g2 = 1/(1/g2-2*slope*width/LOOP)
    v = y[6:15].reshape(3,3)
    np.testing.assert_allclose(v.conj().T@v,np.eye(3),atol=5e-12)
    left,sv,_ = np.linalg.svd(y[15:24].reshape(3,3))
    v0 = ckm(sines,gamma)
    def cp(su,sd,v):
        j = np.imag(v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate())
        def relative(s):
            p = s**2
            return np.prod([p[k]-p[i] for i,k in [(0,1),(0,2),(1,2)]])/np.prod(p)
        return j*relative(su)*relative(sd)
    cp_ratio=cp(y[:3].real,y[3:6].real,v)/cp(u,d,v0)
    np.testing.assert_allclose(cp_ratio,math.exp(y[-1].real),rtol=2e-8)
    angles=[math.atan2(abs(v[0,1]),abs(v[0,0])),math.asin(abs(v[0,2])),
            math.atan2(abs(v[1,2]),abs(v[2,2]))]
    return (np.r_[y[3:6].real,angles],sv[::-1],left[:,::-1],PREV["gamma"](v),cp_ratio)


def align(reference, frame):
    dots = np.diag(reference.conj().T@frame)
    return frame*(dots.conj()/np.abs(dots))


def match(raw,leptons,angle,k):
    short,long,_,delta = CHARGED["get_matrices"](raw,angle)
    size,scale,_ = PREV["branches"](raw,leptons,delta)[2]
    mass = scale*(short+size*np.exp(-1j*k*math.pi/6)*long).T
    basis,sv,_ = np.linalg.svd(mass)
    return math.log(sv[1]/leptons[1]),basis[:,::-1],size


def closest_root(raw,leptons,guess,k):
    bound = np.linspace(2,178,89)
    vals=[match(raw,leptons,x,k)[0] for x in bound]
    pairs=[(a,b) for a,b,fa,fb in zip(bound,bound[1:],vals,vals[1:]) if fa*fb<=0]
    if not pairs:
        return None
    low,high=min(pairs,key=lambda pair:abs(sum(pair)/2-guess))
    f0=match(raw,leptons,low,k)[0]
    for _ in range(34):
        mid=(low+high)/2
        if f0*match(raw,leptons,mid,k)[0]>0:
            low=mid
        else:
            high=mid
    return (low+high)/2


def solve(sines,k,seed,tan):
    frame=np.eye(3,dtype=complex)
    gamma=seed
    for cycle in range(9):
        raw,leptons,left,upper,ratio=run(sines,gamma,frame,tan)
        root=closest_root(raw,leptons,seed,k)
        if root is None:
            return None
        residual,target,size=match(raw,leptons,upper,k)
        target=align(left,target)
        mismatch=np.linalg.norm(left-target)
        if abs(residual)<1e-9 and mismatch<1e-7:
            return dict(gamma_MZ_deg=gamma,gamma_cutoff_deg=upper,abs_t=size,
                        max_left_frame_error=mismatch,middle_mass_log_error=residual,
                        normalized_CP_transport=ratio)
        gamma+=root-upper
        frame=(target@left.conj().T)@frame
    raise AssertionError((gamma,residual,mismatch))


def central(mod_us,mod_cb,ratio):
    ub=mod_cb*ratio
    c13=math.sqrt(1-ub*ub)
    return np.array([mod_us/c13,ub,mod_cb/c13])


def joint_moduli(data, sem, inclusive):
    """Smallest common marginal half-width admitting both |Vcb| inputs and unit row 1."""
    ud, sigma_ud = data["ud_beta"]
    p, sigma_p = data["K_product"]
    form, sigma_f = data["K_form"]
    ratio, *errors = data["ub_over_cb"]
    sigma_r = math.hypot(*errors)
    cb0, *cb_errors = sem["BelleII_Vcb"]
    sigma_cb = math.hypot(*cb_errors)
    cb1, sigma_cb1 = inclusive
    def endpoint(r):
        lo = max(cb0-r*sigma_cb,cb1-r*sigma_cb1)
        hi = min(cb0+r*sigma_cb,cb1+r*sigma_cb1)
        us = (p+r*sigma_p)/(form-r*sigma_f)
        return hi>=lo, (ud+r*sigma_ud)**2+us**2+((ratio+r*sigma_r)*hi)**2, us,hi,ratio+r*sigma_r
    low,high=0.,3.
    assert endpoint(high)[1]>1
    for _ in range(60):
        mid=(low+high)/2
        ok, total, *_ = endpoint(mid)
        if ok and total>=1:
            high=mid
        else:
            low=mid
    ok, value, us, cb, ratio_out=endpoint(high)
    assert ok and abs(value-1)<2e-15
    return high,us,cb,ratio_out


def compute():
    data=read("../1381/inputs.json")
    sem=read("../1445/inputs.json")
    inclusive=read("../1640/inputs.json")["inclusive_holdout"]["Vcb"]
    ratio=data["ub_over_cb"][0]
    direct_us=data["K_product"][0]/data["K_form"][0]
    radius, us_joint, cb_joint, ratio_joint = joint_moduli(data,sem,inclusive)
    published=read("inputs.json")["direct_semileptonic_moduli"]
    cb_main=published["Vcb"][0]
    cases=[("exclusive",direct_us,sem["BelleII_Vcb"][0],ratio),
           ("inclusive",direct_us,inclusive[0],ratio),
           ("CP_even_joint_minimum_width",us_joint,cb_joint,ratio_joint),
           ("main_absolute_moduli",published["Vus"][0],cb_main,published["Vub"][0]/cb_main)]
    raw={}
    for tag,us,cb,rel in cases:
        ss=central(us,cb,rel)
        d=math.sqrt(1-(cb*rel)**2-us**2)
        roots=[solve(ss,k,60.,TAN) for k in (1,2)]
        assert all(r is not None for r in roots)
        if tag=="main_absolute_moduli":
            other=read("results.json")["charged_endpoints"]["direct"]
            for k,r in enumerate(roots,start=1):
                np.testing.assert_allclose(r["gamma_MZ_deg"],other[str(k)]["gamma_MZ_deg"],atol=2e-5,rtol=0)
        raw[tag]=dict(sines=ss.tolist(),derived_ud=d,
                      nuclear_ud_marginal_displacement=(d-data["ud_beta"][0])/data["ud_beta"][1],
                      roots=roots)
    return dict(tan_beta=TAN,kaon_modulus=direct_us,absolute_semileptonic_ratio=ratio,
                minimum_common_CP_even_marginal_width=radius,
                tests=raw)


def stable(val):
    if isinstance(val,dict):
        return {k:stable(v) for k,v in val.items()}
    if isinstance(val,list):
        return [stable(v) for v in val]
    return round(float(val),8) if val is not None else None


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    result=stable(compute())
    dst=HERE/"independent_results.json"
    if cli.parse_args().check:
        assert result==read("independent_results.json")
        print("PASS: independent semileptonic angles, mass-frame RG and commutator integral")
    else:
        with dst.open("x",encoding="utf8") as f:
            f.write(json.dumps(result,indent=2)+"\n")
