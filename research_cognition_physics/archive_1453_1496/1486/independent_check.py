"""Independent SVD and tau-row CP audit, using measured angles after charged rotation."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent


def load(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def window(row):
    return row[0]-2*row[1],row[0]+2*row[-1]


def tau_prob(R):
    S=6*(R+1/R)+16
    r=8/(S+math.sqrt(S*S-1936))
    return r,(2+11*r*r-4*r-5*r/R)/(6*r*(R-1/R))


def tau_range(Rlo,Rhi):
    """Local interval bound: r decreases on the positive small root, R-1/R increases."""
    boundaries=np.linspace(Rlo,Rhi,2049)
    low,high=1.,0.
    for a,b in zip(boundaries,boundaries[1:]):
        u,v=tau_prob(b)[0],tau_prob(a)[0]
        nlo=2+11*u*u-4*v-5*v/a
        nhi=2+11*v*v-4*u-5*u/b
        assert nlo>0
        low=min(low,nlo/(6*v*(b-1/b)))
        high=max(high,nhi/(6*u*(a-1/a)))
    return low-1e-9,high+1e-9


def phase_box(Rrange,solar,reactor):
    t0,t1=tau_range(*Rrange)
    e0,e1=reactor
    u0,u1=solar
    v0,v1=1-t1/(1-e1),1-t0/(1-e0)
    h0=u0*v0+(1-u1)*(1-v1)*e0
    h1=u1*v1+(1-u0)*(1-v0)*e1
    d0=2*math.sqrt(u0*(1-u1)*v0*(1-v1)*e0)
    d1=2*math.sqrt(u1*(1-u0)*v1*(1-v0)*e1)
    values=[(h-1/6)/d for h in (h0,h1) for d in (d0,d1)]
    phase=[math.acos(max(values))/math.pi,math.acos(min(values))/math.pi]
    # To keep these same theta ranges but permit an uncovered positive CP phase.
    x_below=h1-d0*math.cos(.30*math.pi)
    x_above=h0-d0*math.cos(.62*math.pi)
    cost=min(1/6-x_below,x_above-1/6)
    assert .30<phase[0] and phase[1]<.62 and cost>0
    return {"tau_heavy_weight_outer":[t0,t1],"sin2_theta23_outer":[v0,v1],
            "cos_delta_outer":[min(values),max(values)],"positive_delta_in_pi_outer":phase,
            "min_tau1_probability_change_if_same_angle_windows":cost}


def witness(q,s,R,clebsch,phi=0.,eta=2*math.pi/3):
    # Recover the physical two-by-two matrix directly from the right Gram.
    gram=np.diag([q*q,1.])
    Q=np.array([[math.sqrt(1-s*s),s],[-s,math.sqrt(1-s*s)]])
    G=Q @ gram @ Q.T
    b=math.sqrt(G[0,0])
    D=np.array([[0.,q/b],[b,G[0,1]/b]])
    assert np.allclose(D.T@D,G)
    ca,cb,cc=clebsch
    E=D*np.array([[0.,ca],[cb,cc]])
    L,m,_=np.linalg.svd(E)
    se=abs(L[0,0]); ce=math.sqrt(1-se**2)
    a=np.array([0.,1.,1.]); bvec=np.array([1.,3.,1.])
    r=tau_prob(R)[0]
    M=np.outer(a,a)+np.exp(1j*eta)*r*np.outer(bvec,bvec)
    _,U=np.linalg.eigh(M.conj().T @ M)
    charge=np.array([[ce,-se*np.exp(1j*phi),0],[se*np.exp(-1j*phi),ce,0],[0,0,1.]])
    U=charge @ U
    p=abs(U)**2
    j=(U[0,0]*U[1,1]*U[0,1].conj()*U[1,0].conj()).imag
    s13=p[0,2]; s12=p[0,1]/(1-s13); s23=p[1,2]/(1-s13)
    scale=(1-s13)*math.sqrt(s13*s12*(1-s12)*s23*(1-s23))
    c=(s12*s23+(1-s12)*(1-s23)*s13-p[2,0])/(2*math.sqrt(s12*(1-s12)*s23*(1-s23)*s13))
    assert np.isclose(p[2,0],1/6) and np.isclose(p[2,2],tau_prob(R)[1])
    assert np.isclose((j/scale)**2+c*c,1)
    return {"charged_ratio":m[0]/m[1],"theta_e_degrees":math.degrees(math.asin(se)),
            "reactor_amplitude":4*s13*(1-s13),"solar":s12,"atmospheric":s23,
            "J":j,"PDG_delta_over_pi":math.atan2(j/scale,c)/math.pi}


def compute():
    spec=load(HERE/"inputs.json")
    masses=load(HERE/spec["quark_lepton"])["lattice"]
    osc=load(HERE/spec["oscillations"])
    s=next(row[1] for row in load(HERE/spec["charm"])["edges"] if row[0]=="cd")
    q=masses["md"][0]/masses["ms"][0]
    dm=window(osc["KamLAND_dm21_eV2"])
    atm=window(osc["T2K_dm32_NO_eV2"])
    minos=window(osc["MINOS_effective_atm_eV2"])
    R0=math.sqrt(1+osc["T2K_dm32_NO_eV2"][0]/osc["KamLAND_dm21_eV2"][0])
    solar=window(osc["JUNO_sin2_theta12_NO"])
    reactor=[(1-math.sqrt(1-x))/2 for x in window(osc["DayaBay_sin2_2theta13"])]
    base={k:witness(q,s,R0,z) for k,z in spec["clebsch"].items()}
    mirror=witness(q,s,R0,spec["clebsch"]["small_left"],eta=-2*math.pi/3)
    assert np.isclose(base["small_left"]["PDG_delta_over_pi"],-mirror["PDG_delta_over_pi"])
    return {"direct_spectral_checks":base,"small_left_conjugate":mirror,
            "T2K_KamLAND":phase_box([math.sqrt(1+atm[0]/dm[1]),math.sqrt(1+atm[1]/dm[0])],solar,reactor),
            "MINOS_KamLAND":phase_box([math.sqrt(minos[0]/dm[1]),math.sqrt(1+minos[1]/dm[0])],solar,reactor)}


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    return np.allclose(a,b,atol=1e-12,rtol=1e-10)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    options=parser.add_mutually_exclusive_group()
    options.add_argument("--write",action="store_true")
    options.add_argument("--check",action="store_true")
    value=compute()
    target=HERE/"independent_results.json"
    if parser.parse_args().write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(value,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        assert equal(value,load(target))
    print(json.dumps(value,indent=2))
