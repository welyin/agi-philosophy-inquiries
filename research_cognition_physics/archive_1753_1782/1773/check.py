"""Use one nuclear modulus, two B phases and the charged rare K rate to invert CKM."""
import argparse
import cmath
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
def read(file):
    return json.loads((HERE/file).read_bytes())


CFG = read("inputs.json")
K = read(CFG["kaon_short_distance"])
D = read(CFG["direct_flavour"])
B = read(CFG["B_phases"])["data"]
KP = K["kplus_e11"][0]*(1+K["delta_EM"])
KL = K["kL_e11"][0]
REF = K["lambda_reference"]
TOP = K["Xtop"][0]/REF**4
CHARM = K["Pcharm"][0]
EPS = runpy.run_path(str(HERE/CFG["epsilon_analysis"]))


def root(fun, lo=0., hi=.01):
    assert fun(lo)<0<fun(hi)
    for _ in range(65):
        mid=(lo+hi)/2
        if fun(mid)<0:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def projected(rho, theta, top=TOP, charm=CHARM):
    ph=rho*cmath.exp(-1j*theta)
    zc=-1/(1+ph)
    zt=ph*zc
    return top*zt+charm*zc.real, zt, zc


def from_moduli(mod, gamma_deg):
    us,ub,cb=(mod[k] for k in ("Vus","Vub","Vcb"))
    c13=math.sqrt(1-ub**2)
    s12,s23=us/c13,cb/c13
    c12,c23=math.sqrt(1-s12*s12),math.sqrt(1-s23*s23)
    gamma=math.radians(gamma_deg)
    delta=gamma+math.asin(c12*s23*ub*math.sin(gamma)/(s12*c23))
    e=cmath.exp(1j*delta)
    return np.array([[c12*c13,s12*c13,ub/e],
        [-s12*c23-c12*s23*ub*e,c12*c23-s12*s23*ub*e,s23*c13],
        [s12*s23-c12*c23*ub*e,-c12*s23-s12*c23*ub*e,c23*c13]])


def readings(v):
    ud,us,ub=v[0]
    cd,cs,cb=v[1]
    td,ts,tb=v[2]
    lam_u=us.conjugate()*ud
    lam_c,lam_t=cs.conjugate()*cd,ts.conjugate()*td
    zt,zc=lam_t/lam_u,lam_c/lam_u
    beta=cmath.phase(-cd*cb.conjugate()/(td*tb.conjugate()))
    beta_s=cmath.phase(-ts*tb.conjugate()/(cs*cb.conjugate()))
    gamma=cmath.phase(-ud*ub.conjugate()/(cd*cb.conjugate()))
    z=abs(ud)*(TOP*zt+CHARM*zc.real)
    # Evaluate the published formula with the actual |Vus| separately as a convention test.
    gauge=lam_u.conjugate()/abs(lam_u)
    physical_lam=abs(us)
    scaled_top=K["Xtop"][0]*lam_t*gauge/physical_lam**5
    scaled_charm=(lam_c*gauge).real/physical_lam*CHARM*(REF/physical_lam)**4
    scale=(physical_lam/REF)**8
    np.testing.assert_allclose([KP*abs(z)**2, KL*z.imag**2],
         [KP*scale*abs(scaled_top+scaled_charm)**2, KL*scale*scaled_top.imag**2],rtol=3e-12)
    np.testing.assert_allclose(zt/zc,abs(zt/zc)*cmath.exp(-1j*(beta+beta_s)),atol=2e-13)
    lt=lam_t*gauge
    weights=np.array([-(lt**2).imag/2, abs(lam_u)*lt.imag])
    epsilon=float(weights@EPS["COEFFICIENTS"])
    epsilon_width=math.sqrt(weights@EPS["COV"]@weights)
    return dict(moduli=np.abs(v).tolist(), beta_deg=math.degrees(beta),
                beta_s_deg=math.degrees(beta_s), gamma_deg=math.degrees(gamma),
                J_abs=abs((lam_u.conjugate()*lam_t).imag),rho=abs(zt/zc),
                B_plus_e11=KP*abs(z)**2, B_long_e11=KL*z.imag**2,
                eps_K=epsilon, eps_K_coefficient_width_at_fixed_CKM=epsilon_width,
                eps_K_pull_at_fixed_CKM=(epsilon-EPS["EPS_EXP"])/math.hypot(epsilon_width,EPS["EPS_ERROR"]))


def ckm_completions(zt, zc, ud, beta_s):
    p=zt.conjugate()*zc
    l=p.imag/math.tan(beta_s)-p.real
    assert l>0
    dt=1-ud**2
    t,c=abs(zt)**2,abs(zc)**2
    # x=|Vtd|^2. All physical roots are tested.
    poly=[-1/ud**2,dt/ud**2-t+c,t*dt-t*c/l]
    out=[]
    for x in np.roots(poly):
        if abs(x.imag)>1e-11 or not 0<x.real<dt:
            continue
        x=float(x.real)
        us2=l*x*(dt-x)/(ud**2*t*c)
        if not 0 < us2 <= dt:
            continue
        ds=np.array([ud,math.sqrt(dt-x),math.sqrt(x)])
        ss=ud*math.sqrt(us2)*np.array([1,zc,zt]).conj()/ds
        bs=np.cross(ds,ss).conj()
        v=np.column_stack([ds,ss,bs])
        np.testing.assert_allclose(v.conj().T@v,np.eye(3),atol=2e-12)
        row=readings(v)
        np.testing.assert_allclose(math.radians(row["beta_s_deg"]), beta_s,atol=2e-12)
        out.append(row)
    return out


def infer(branch,ud,S,phi,charged):
    b0=.5*math.asin(S)
    beta=b0 if branch=="acute" else math.pi/2-b0
    bs=-phi/2
    theta=beta+bs
    rho=root(lambda r:KP*ud**2*abs(projected(r,theta)[0])**2-charged)
    _,zt,zc=projected(rho,theta)
    out=ckm_completions(zt,zc,ud,bs)
    assert len(out)==1
    row=out[0]
    np.testing.assert_allclose([row["B_plus_e11"],row["beta_deg"]],[charged,math.degrees(beta)],
                               rtol=2e-10)
    return row


def box(branch,n,phase_budget=0.):
    ud,du=D["holdout_nuclear"]
    S,dS=B["golden_S"]; phi,dphi=B["phi_s_ccs"]
    charged,lo,hi=CFG["kaon_charged_e11"]
    aa=[.5*math.asin(S-n*dS),.5*math.asin(S+n*dS)]
    if branch!="acute":
        aa=[math.pi/2-aa[1], math.pi/2-aa[0]]
    th_lo=aa[0]-(phi+n*dphi)/2-math.radians(phase_budget)
    th_hi=aa[1]-(phi-n*dphi)/2+math.radians(phase_budget)
    assert 0<th_lo<th_hi<math.pi/2
    pl,pu=[CHARM+i*n*K["Pcharm"][1] for i in (-1,1)]
    al,au=[TOP+i*n*K["Xtop"][1]/REF**4 for i in (-1,1)]
    kp_min,kp_max=[KP+i*n*K["kplus_e11"][1]*(1+K["delta_EM"]) for i in (-1,1)]
    rad_min=math.sqrt((charged-n*lo)/kp_max)/(ud+n*du)
    rad_max=math.sqrt((charged+n*hi)/kp_min)/(ud-n*du)
    def envelope(r,low):
        t,p=(al,pl) if low else (au,pu)
        rangle,iangle=(th_hi,th_lo) if low else (th_lo,th_hi)
        return abs(complex(projected(r,rangle,t,p)[0].real,projected(r,iangle,t,p)[0].imag))
    rl=root(lambda r: envelope(r,False)-rad_min)
    ru=root(lambda r: envelope(r,True)-rad_max)
    loamp=projected(rl,th_lo,al,pl)[0].imag*(ud-n*du)
    hiamp=projected(ru,th_hi,au,pu)[0].imag*(ud+n*du)
    return [(KL-n*K["kL_e11"][1])*loamp**2,(KL+n*K["kL_e11"][1])*hiamp**2]


def calculate():
    rows={}
    for key in ("acute","reflected"):
        row=infer(key,D["holdout_nuclear"][0],B["golden_S"][0],B["phi_s_ccs"][0],
                  CFG["kaon_charged_e11"][0])
        row["one_width_outer_KL_e11"]=box(key,1)
        row["two_width_outer_KL_e11"]=box(key,2)
        row["one_width_plus_phase_outer_KL_e11"]=box(key,1,CFG["additional_net_B_phase_budget_degrees"])
        rows[key]=row
    frozen={}
    direct=read(CFG["matter_points"])
    for angle, field in [(62.5,"refits"),(64.,"exploratory_64_deg")]:
        frozen[str(angle)]={name:readings(from_moduli(row["moduli"],angle)) for name,row in direct[field].items()}
    other={key:readings(from_moduli(r["weak"],r["gamma"])) for key,r in
           read(CFG["alternate_matter_points"]).items() if isinstance(r,dict) and "gamma" in r}
    # |J| from KL vs maximum of the ud-cb-ub-cd quartet in the direct B box.
    n=2
    ud,du=D["holdout_nuclear"]
    ub=D["inclusive_ub"][0]+n*math.hypot(*D["inclusive_ub"][1:])
    cb=D["inclusive_cb"][0]+n*D["inclusive_cb"][1]
    us_min_sq=1-(ud+n*du)**2-ub**2
    light=rows["reflected"]["two_width_outer_KL_e11"][0]
    j_min=(ud-n*du)*us_min_sq/(TOP+n*K["Xtop"][1]/REF**4)*math.sqrt(light/(KL+n*K["kL_e11"][1]))
    j_max=(ud+n*du)*ub*cb*math.sqrt(1-(ud-n*du)**2)
    assert j_min>j_max
    return dict(inverse=rows, frozen_matter=frozen, prior_K_calibrated_candidates=other,
                reflected_direct_B_exclusion_in_two_width_box=dict(
                    J_required_at_least=j_min,J_allowed_at_most=j_max,ub_upper=ub,cb_upper=cb))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    output=json.loads(json.dumps(calculate()), parse_float=lambda v:float(f"{float(v):.10g}"))
    path=HERE/"results.json"
    if args.check:
        assert json.loads(path.read_bytes())==output
        print("1773 OK: exact CKM completions, physical kaon normalization and phase sensitivity")
    else:
        path.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",encoding="utf8")
