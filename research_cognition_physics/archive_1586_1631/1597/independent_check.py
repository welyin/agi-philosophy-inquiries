"""Electrical alpha and signed finite radiative sources for the 1597 review."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def signed_gap(ae_box, amu_box, me, mm, other_e=0., other_mu=0.):
    lo = (ae_box[0]-other_e)/(2*me**2) - (amu_box[1]+other_mu)/(2*mm**2)
    hi = (ae_box[1]+other_e)/(2*me**2) - (amu_box[0]-other_mu)/(2*mm**2)
    return max(0.,lo,-hi)


def mass_floor(gap, case, cfg, mix):
    x = cfg["s12"][1]
    y = cfg["s13"][2]
    z = mix[2]
    rmin = cfg["atm_abs"][1]
    beta2 = (1-y)*(1-(1-x)*cfg["solar"][2]/rmin)
    a = math.sqrt(y/((1-y)*(1-z)))
    b = math.sqrt(z/(1-z))
    s12max = cfg["s12"][2]
    generic = (1-y)*math.sqrt(rmin)*((1-s12max)*math.sqrt(1-cfg["solar"][2]/rmin)-s12max)
    if not gap:
        return generic
    h = (case[0]+case[1])/gap
    # Increasing row-norm polynomial, solved by bisection.
    lo, hi = 0., 1.
    for _ in range(70):
        mid = (lo+hi)/2
        if mid**2+h**2+(a*mid+b*h)**2 < beta2:
            lo = mid
        else:
            hi = mid
    return max(generic,lo*math.sqrt(rmin))


def texture(angles, mix):
    s12, s13, s23 = np.sqrt([angles["s12"][0], angles["s13"][0], mix[0]])
    c12, c13, c23 = np.sqrt(1-np.array([s12,s13,s23])**2)
    m2 = math.sqrt(angles["atm_abs"][0])
    m1 = math.sqrt(m2**2-angles["solar"][0])
    def mixing(delta):
        t = np.exp(1j*delta)
        return np.array([[c12*c13,s12*c13,s13/t],
                         [-s12*c23-c12*s23*s13*t,c12*c23-s12*s23*s13*t,s23*c13],
                         [s12*s23-c12*c23*s13*t,-c12*s23-s12*c23*s13*t,c23*c13]])
    def balance(delta):
        U = mixing(delta)
        return abs(m1*U[0,0]*U[1,0]) - abs(m2*U[0,1]*U[1,1])
    lo,hi = 0.,math.pi
    assert balance(lo)*balance(hi)<0
    for _ in range(70):
        mid = (lo+hi)/2
        if balance(lo)*balance(mid)>0:
            lo=mid
        else:
            hi=mid
    delta = -(lo+hi)/2
    U = mixing(delta)
    phase = -m1*U[0,0]*U[1,0]/(m2*U[0,1]*U[1,1])
    assert abs(abs(phase)-1)<1e-12
    N = (U*np.array([m1,m2*phase,0.])) @ U.T
    n=U[:,2]
    left=np.vstack([np.eye(2),-n[:2].conj()/n[2].conj()])
    assert np.max(np.abs(left @ np.diag(np.diag(N)[:2]) @ left.T-N))<1e-15
    return N,left,dict(dirac_deg=math.degrees(delta)%360,
                       majorana_column_squared_phase_deg=math.degrees(np.angle(phase)),
                       mee_meV=float(1000*abs(N[0,0])))


def finite_spectrum(target_k, source_mass, neutral_coefficient):
    # CP-even/odd neutral scalars share the same eta and complex-singlet
    # mixing. Weak-charge assignment Q(E)=-1 gives the negative kernel.
    v=246.22
    lam5=1e-6
    t,w=np.polynomial.legendre.leggauss(72)
    t,w=(t+1)/2,w/2
    sin2=.3
    s=math.sqrt((1-math.sqrt(1-sin2**2))/2)
    c=math.sqrt(1-s*s)
    O=np.array([[c,s],[-s,c]])
    def spectrum(scale):
        E=1300.*scale
        MN=1700.*scale
        A=O @ np.diag((scale*np.array([1100.,1500.]))**2) @ O.T
        Is=[]; Js=[]; fs=[]
        poles=[]
        for sig in (-1,1):
            vals,V=np.linalg.eigh(A+np.diag([sig*lam5*v**2/2,0.]))
            poles.append(np.sqrt(vals).tolist())
            d0=(1-t)*E**2+t*vals[0]
            d1=(1-t)*E**2+t*vals[1]
            Is.append(np.dot(w,np.log1p(t*(vals[1]-vals[0])/d0))/(vals[1]-vals[0]))
            Js.append(-np.dot(w,t*(1-t)**2/(d0*d1)))
            x=vals/MN**2
            fs.append(np.dot(V[0]**2,x*np.log(x)/(x-1)))
        k=sum(Js)/sum(Is)
        charged=E*abs(A[0,1])*sum(Is)/(32*math.pi**2)
        neutral=MN*(fs[1]-fs[0])/(32*math.pi**2)
        return k,charged,neutral,E,MN,poles
    lo,hi=.03,15.
    assert spectrum(lo)[0]<target_k<spectrum(hi)[0]
    for _ in range(65):
        mid=(lo+hi)/2
        if spectrum(mid)[0]<target_k:
            lo=mid
        else:
            hi=mid
    k,charged,neutral,E,MN,poles=spectrum((lo+hi)/2)
    assert abs(k/target_k-1)<1e-12 and neutral>0
    prod=source_mass/charged
    ney=math.sqrt(abs(neutral_coefficient)*1e-9/neutral)
    return dict(k_GeV_minus2=k,E_GeV=E,N_GeV=MN,two_neutral_scalar_pairs_GeV=poles,
                charged_yukawa_product=prod,neutral_yukawa_modulus=ney,
                neutral_yukawa_phase_deg=math.degrees(np.angle(neutral_coefficient))/2,
                lambda5=lam5)


def compute():
    em=read(HERE.parent/"1591"/"inputs.json")
    me,mm=em["me_GeV"],em["mmu_GeV"]
    angles=read(HERE.parent/"1370"/"inputs.json")["IO"]
    mix=read(HERE.parent/"1580"/"inputs.json")["theta23_sin2"]["IO"]
    cases=read(HERE.parent/"1596"/"results.json")["cases"]
    lfv=cases["finite_LFV"]["offdiag_allowance_GeV_minus2"]
    # WP25 eq. 7.6, 7.10, 7.27. Use individual alpha rectangles independently
    # of the published combined uncertainty, and never combine those branches.
    ue=137.035999163
    uc=137.035999045
    dmu=np.array([-88.,164.])*1e-11
    slope=(8.4e-6,8.52e-6)
    nat=2/(4e-7*math.pi*299792458.)
    uq,sq=25812.80831*nat,.00062*nat
    result={"NIST_inverse_alpha": [uq,sq]}
    def lfv_norm(br, small, large, lifetime):
        return math.sqrt(4*em["hbar_GeV_s"]*br*ue/
                         (lifetime*large**3*(1-(small/large)**2)**3))
    me_to_mu=lfv_norm(1.5e-13,me,mm,em["mu_lifetime_s"])
    e_to_tau=lfv_norm(3.3e-8,me,1.775,2.88e-13)
    low13=angles["s13"][1]
    bound=2*me*math.sqrt((1-low13)/low13)*(
        math.sqrt(mix[2])*me_to_mu+math.sqrt(1-mix[1])*e_to_tau)
    assert (uq-ue)*slope[0]>bound/(math.sqrt(.8)*.8)
    result["null_electron_column_max_abs_ae"]=[bound,bound/(math.sqrt(.8)*.8)]
    ae_cases={}
    for lab,u,sig in [("NIST",uq,sq),("Cs",uc,27e-9)]:
        for n in (1,2):
            delta=[(u-ue)+j*n*(sig+15e-9) for j in (-1,1)]
            vals=[x*y for x,y in itertools.product(delta,slope)]
            ae_cases[lab+str(n)]=[min(vals),max(vals)]
    ae_cases["WP25_published_Cs_two_width"]=[-1.52e-12,-.48e-12]
    ae_cases["historical_1591_Cs_two_width"]=[-1.739e-12,-.3024e-12]
    for lab,box in ae_cases.items():
        # WP25 QED input transfer for the wider electrical region bounded
        # conservatively by |delta u|*1e-5 in the nonhadronic reference.
        pad=1.4e-10 if lab.startswith("NIST") else 3e-12
        mu=dmu+np.array([-pad,pad])
        g=signed_gap(box,mu,me,mm)
        gdef=signed_gap(box,mu,me,mm,2e-14,2e-10)
        result[lab]=dict(ae_box=box,g_abs_min=g,meV_min=1000*mass_floor(g,[lfv,0],angles,mix),
                         with_model_defect_g_abs_min=gdef,with_model_defect_meV_min=1000*mass_floor(gdef,[lfv,1.e-8],angles,mix))
    # Two free dipoles can hold both central measurements as alpha varies.
    # The common new source z* must solve the two residual equations.
    de,dm=8.46e-6,8.55e-6
    amu0=38e-11
    ustar=amu0/(de*mm**2/me**2-dm)    # offsets in inverse alpha from trap
    common=de*ustar/(2*me**2)
    assert np.isclose(common,(amu0+dm*ustar)/(2*mm**2),rtol=1e-12,atol=0)
    assert abs(ue+ustar-uq)<2*sq
    result["two_free_anomalies_equal_normalized_response"]={"inverse_alpha_offset_from_trap":ustar,"z_GeV_minus2":common}
    N,L,detail=texture(angles,mix)
    # A signed-Cs leading matching with both physical kernels negative.
    ae=-1.e-12
    k=np.array([ae/(2*me**2),-1.e-8])
    m=np.array([me,mm,1.77686])
    rad=L @ np.diag(m[:2]) @ np.eye(2,3)
    T=np.diag(m)-rad
    C=L @ np.diag(k*m[:2]) @ np.eye(2,3)
    X=C @ np.diag(1/m)
    assert np.linalg.matrix_rank(T,tol=1e-13)==1 and np.max(abs(X@N-N@X.T))<1.e-19
    pref=em["alpha"]*m[2]**3/(4*em["hbar_GeV_s"]/2.903e-13)
    br=[pref*(1-(m[i]/m[2])**2)**3*(abs(C[i,2])**2+abs(C[2,i])**2) for i in (0,1)]
    assert br[0]<3.3e-8 and br[1]<4.4e-8 and C[0,1]==C[1,0]==0
    detail.update({"delta_ae_amu":(2*m[:2]**2*k).tolist(),"B_tauegamma_taumugamma":br,
                   "negative_kernel_spectra":[finite_spectrum(ki,mi,N[i,i]) for i,(ki,mi) in enumerate(zip(k,m[:2]))]})
    result["signed_common_source_example"]=detail
    return result


def same(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(same(v,b[k]) for k,v in a.items())
    if isinstance(a,list):
        return len(a)==len(b) and all(same(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=2.e-6,abs_tol=2.e-19)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=compute()
    path=HERE/"independent_results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    elif not same(result,read(path)):
        raise SystemExit("independent values differ")
    print(json.dumps(result))
