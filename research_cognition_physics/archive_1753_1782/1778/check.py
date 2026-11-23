"""Scalar spectral moment and finite-window tests with an unpaired chiral subspace."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
def read(name):
    return json.loads((HERE/name).read_bytes())


CFG=read("inputs.json")
OLD=read(CFG["registered_light_point"])["matching"]
V=OLD["vev_GeV"]
NC=CFG["colours"]


def moments():
    up=np.array(OLD["high_yu"])*V/math.sqrt(2)
    down=np.array([OLD["gate_d"]["light_GeV"],OLD["daughter_yd"][1]*V/math.sqrt(2),
                   OLD["gate_b"]["light_GeV"]])
    out={}
    for floor in CFG["hypothetical_new_mass_floors_GeV"]:
        low=np.r_[up,floor,down,floor]
        y=math.sqrt(2)*low/V
        t2,t4=NC*sum(y**2),NC*sum(y**4)
        out[str(floor)]=dict(T2_min=t2,T4_min=t4,abs_box_term_in_dlambda_dlogmu_min=2*t4/(16*math.pi**2))
    previous=read(CFG["old_fourth_thresholds"])["required_projected_Y_if_only_three_below"]
    assert abs(math.sqrt(2)*1700/V-previous["1700"]["minimum_Y"])<3.e-5
    return dict(light_quark_masses_at_threshold_GeV=[up.tolist(),down.tolist()],floor_rows=out)


def window(floor,A,D,E,Q):
    if E<=D or 2*E>=Q:
        return 0.,0.
    b=max(0.,floor-A*D/math.sqrt(E*E-D*D))
    kinematic=(1-(2*E/Q)**2)**1.5
    return b, NC*2*(b/V)**2*kinematic


def scalar_rate(m,B,Q,known):
    """Dimensionless tree spectral ratio divided by one massless unit-Yukawa Dirac channel."""
    y=math.sqrt(2)*B/V
    x=(m/Q)**2
    a,b=x[:,None],x[None,:]
    mass=m[:,None]+m[None,:]
    threshold=mass<Q
    tagged=np.ones(y.shape,bool)
    tagged[:known,:known]=False
    kernel=np.sqrt(np.maximum(0.,(1-a-b)**2-4*a*b))
    rate=kernel*((1-a-b)*(abs(y)**2+abs(y.T)**2)/2-2*np.sqrt(a*b)*np.real(y*y.T))
    assert np.min(rate[threshold],initial=0)>-1e-8
    return NC*float(sum(rate[threshold&tagged]))


def counter(k):
    cfg=CFG["blind_texture"]
    m=cfg["extra_light_GeV"]
    b=m*cfg["b_over_m"]
    H=m*(1+2*k*k)/(2*k*k)
    B=np.array([[m-b*k,b,0.],[b,m-H-b/k,0.],[0.,0.,m-H]])
    M=B.copy()
    M[1,2]=M[2,1]=H
    ev,U=np.linalg.eigh(M)
    idx=np.argsort(abs(ev))
    ev,U=ev[idx],U[:,idx]
    masses=np.abs(ev)
    bm=U.T@B@(U*np.sign(ev))
    v=np.array([1.,k,k])
    np.testing.assert_allclose(M@v,m*v,atol=2.e-9)
    assert abs(v@B@v)<1.e-10 and abs(bm[0,0])<2.e-8
    np.testing.assert_allclose(masses[0],m,rtol=1.e-10)
    Q=cfg["test_virtuality_GeV"]
    assert masses[1]>Q
    rate=scalar_rate(masses,bm,Q,0)
    assert rate<1.e-19
    return dict(k=k,bare_common_doublet_mass_GeV=H,masses_GeV=masses.tolist(),
                protected_budget_GeV=np.linalg.norm(B[:,0]),full_budget_GeV=np.linalg.norm(B,2),
                analytic_light_Higgs_coupling=0,analytic_extra_tree_scalar_rate_at_Q=0)


def random_checks():
    rng=np.random.default_rng(CFG["random_seed"])
    def complex_noise(shape):
        return rng.normal(size=shape)+1j*rng.normal(size=shape)
    positive=0
    active_rank_two=0
    for count in range(CFG["random_cases"]):
        p=4+count%2
        # Columns ordered (ordinary right singlets with mirrors, right doublet).
        n=p+2
        B=np.zeros((n,n),complex)
        B[:-1,:-1]=150*np.eye(n-1)+25*complex_noise((n-1,n-1))
        B[-1,-1]=20*complex_noise((1,))[0]
        M=B.copy()
        M[-2,-1]=rng.uniform(150,1500)
        M[-1,-2]=rng.uniform(150,3000)
        left,m,right_h=np.linalg.svd(M)
        m=m[::-1]
        left=left[:,::-1]
        right=right_h.conj().T[:,::-1]
        bm=left.conj().T@B@right
        D=np.linalg.norm(B[:p,:],2)
        A=np.linalg.norm(B,2)
        d=np.linalg.svd(B[:p,:],compute_uv=False)[::-1]
        assert np.min(d-m[:p])>-1.e-8
        s=np.linalg.svd(B,compute_uv=False)
        assert sum(s**4)+1.e-3>=sum(m[:p]**4)
        E=rng.uniform(2,9)*D
        Q=3.1*E
        b,rate=window(m[3],A,D,E,Q)
        keep=m<=E
        extra=(np.arange(n)>=3)&keep
        part=bm[np.ix_(extra,keep)]
        singular=np.linalg.svd(part,compute_uv=False)
        assert (len(singular)>=p-3) and singular[p-4]+1.e-8>=b
        assert scalar_rate(m,bm,Q,3)+1.e-9>=(p-3)*rate
        positive+=int(b>0)
        active_rank_two+=int(b>0 and p==5)
    return dict(cases=CFG["random_cases"],nonzero_bounds=positive,two_extra_directions_with_positive_bound=active_rank_two)


def main():
    cfg=CFG["finite_example"]
    A=V*cfg["full_max_Yukawa"]/math.sqrt(2)
    D=V*cfg["protected_max_Yukawa"]/math.sqrt(2)
    rows=[]
    for Q,E in cfg["virtuality_and_mass_window_GeV"]:
        b,rate=window(cfg["extra_mass_floor_GeV"],A,D,E,Q)
        rows.append(dict(Q_GeV=Q,E_GeV=E,extra_singular_strength_GeV=b,
                         relative_rate_one_quark_sector=rate,relative_rate_both_sectors=2*rate))
    m=300.
    q=3000.
    np.testing.assert_allclose(scalar_rate(np.array([m]),np.array([[m]]),q,0),
                               NC*2*(m/V)**2*(1-4*(m/q)**2)**1.5,rtol=1.e-14)
    return dict(moments=moments(),bounded_full_Y=rows,
                unbounded_full_Y_counterexamples=[counter(k) for k in CFG["blind_texture"]["k"]],
                random_tests=random_checks())


def tidy(x):
    if isinstance(x,dict):
        return {str(k):tidy(v) for k,v in x.items()}
    if isinstance(x,list):
        return [tidy(v) for v in x]
    return float(f"{x:.10g}") if isinstance(x,(float,np.floating)) else x


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    verify=parser.parse_args().check
    data=tidy(main())
    path=HERE/"results.json"
    if verify:
        assert data==read("results.json")
        print("1778 verified: full Yukawa moments, finite scalar spectrum, and blind chiral threshold")
    else:
        path.write_text(json.dumps(data,indent=2)+"\n",encoding="utf8")
