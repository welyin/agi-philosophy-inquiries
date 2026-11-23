"""Gauge-invariant doublet chains with finite mass thresholds in a single matter model."""
from argparse import ArgumentParser
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
TWOPI=2*math.pi


def fetch(path):
    return json.loads((HERE/path).read_text(encoding="utf8"))


def scalar_spectrum(n,epsilon):
    lengths=2*n+1
    chain=np.zeros((2*n,lengths))
    for branch in range(2):
        parent=0
        for step in range(n):
            child=1+branch*n+step
            chain[child-1,parent]=-epsilon
            chain[child-1,child]=1
            parent=child
    profile=np.array([1]+[epsilon**i for i in range(1,n+1)]*2)
    norm=np.linalg.norm(profile)
    assert np.linalg.norm(chain@profile)<1e-14
    heavy=np.linalg.svd(chain,compute_uv=False)
    assert len(heavy)==2*n and min(heavy)>0.5
    assert np.isclose(np.prod(heavy),norm)
    return heavy,norm


def ckm(angles):
    x,y,z,delta=angles
    a,b,c=np.sin([x,y,z])
    u,v,w=np.cos([x,y,z])
    p=np.exp(1j*delta)
    return np.array([[u*v,a*v,b/p],[-a*w-u*c*b*p,u*w-a*c*b*p,c*v],
                     [a*c-u*w*b*p,-u*c-a*w*b*p,w*v]])


def physical_match(rows,charges,eps,norm):
    vals=rows[:3]
    theoretical=[]
    left=[]
    for sector,indices,frame in zip(vals,("uc","dc","ec"),(np.eye(3),ckm(rows[3]),np.eye(3))):
        overlaps=eps**np.asarray(charges[indices])/norm
        target=frame@np.diag(sector)
        core=target/overlaps
        measured=core*overlaps
        assert np.allclose(measured,target,rtol=1e-12,atol=1e-15)
        # Charged masses are pairwise distinct. Verify the flavor left-projector frame.
        h=measured@measured.conj().T
        assert np.allclose(frame.conj().T@h@frame,np.diag(np.asarray(sector)**2),atol=1e-15)
        left.append(frame)
        theoretical.extend(np.linalg.svd(core,compute_uv=False))
    assert np.allclose(left[0].conj().T@left[1],ckm(rows[3]))
    return [float(min(theoretical)),float(max(theoretical))]


def gauge_data(charges,n):
    # one-loop increments in the positive screening sign, with physical SM hypercharge.
    heavy_f=np.zeros(3)
    abelian=np.zeros((2,2))
    reps=[("Q",3,2,1/6,1),("uc",3,1,-2/3,1),("dc",3,1,1/3,1),
          ("L",1,2,-1/2,0),("ec",1,1,1,0)]
    def add(multiplicity, y, f):
        v=np.array([y,f])
        abelian[:]+=multiplicity*np.outer(v,v)
    for name,color,weak,y,tc in reps:
        for x in charges[name]:
            add(2*color*weak/3,y,x)
            for k in range(1,x+1):
                add(2*color*weak/3,y,k-1)
                add(2*color*weak/3,-y,-k)
            heavy_f += x*np.array([4*color*weak*y*y/3,
                                   2*color*(weak==2)/3, 2*tc*weak/3])
    for k in range(-n,n+1):
        add(2/3,1/2,k)
    add(1/3,0,-1)
    b_sm = np.array([41/6,-19/6,-7])
    b_scalar = np.array([1/6,1/6,0])
    total=b_sm+heavy_f+2*n*b_scalar
    assert np.isclose(total[0],abelian[0,0])
    return b_sm,heavy_f,b_scalar,abelian


def inv_alpha(gauge):
    g1,g2,g3=gauge
    return 4*math.pi/np.array([3*g1*g1/5,g2*g2,g3*g3])


def run_gauge(couplings,ref,f,stop,hs,eps,beta):
    sm,ferm,scal,abelian=beta
    heavy=hs*f/eps
    outer=max(max(heavy),f)
    def at(energy):
        new_f=max(0.,math.log(energy/f))
        new_s=sum(max(0.,math.log(energy/m)) for m in heavy)
        return couplings-(sm*math.log(energy/ref)+ferm*new_f+scal*new_s)/TWOPI
    assert ref<=f<=stop and min(heavy)>f
    return at(stop),at(outer),outer


def run():
    cfg=fetch("inputs.json")
    charge=fetch(cfg["charges"])["quark_charge_control"]
    low=fetch(cfg["gauge_boundary"])
    high=fetch(cfg["high"])["high_scale"]
    mid=cfg["at_1e12"]
    rows=[[mid["up"],mid["down"],mid["charged_leptons"],mid["angles"]],
          [[a[0] for a in high["up"]],[a[0] for a in high["down"]],
           [a[0] for a in fetch(cfg["high_leptons"])["scales"]["1e16"][3:]],
           [a[0] for a in high["angles_theta12_theta13_theta23_delta"]]]]
    beta=gauge_data(charge,cfg["D_chains"])
    output={"Abelian_b_above_all":np.round(beta[-1],10).tolist(),"epsilon_examples":[]}
    low_g=[low[name][0] for name in ("g1_GUT","g2","g3")]
    for eps in cfg["epsilons"]:
        hs,norm=scalar_spectrum(cfg["D_chains"],eps)
        numbers = {"epsilon":eps,"null_profile_norm":norm,
                   "heavy_min_max_over_Lambda":[float(min(hs)),float(max(hs))],
                   "cores_at_1e12_and_1e16":[physical_match(r,charge,eps,norm) for r in rows]}
        threshold={}
        # Complete spectrum here is below the upper endpoint near each boundary.
        for i,label in enumerate(["hypercharge","weak","color"]):
            def passes(f):
                v=run_gauge(inv_alpha(low_g),low["mu_ref_GeV"],f,cfg["reach_GeV"],hs,eps,beta)[0]
                return v[i]>=1/cfg["alpha_control"]
            lo,hi=1e4,1e15
            if passes(lo):
                threshold[label]=None
                continue
            assert passes(hi)
            for _ in range(70):
                m=math.sqrt(lo*hi)
                if passes(m):
                    hi=m
                else:
                    lo=m
            threshold[label]=hi
        numbers["min_f_at_MZ_boundary_for_alpha_below_0p30_at_1e16"]=threshold
        examples=[]
        for f,stop,g in [(1e12,1e16,mid["gauge_GUT1_2_3"]),
                          (1e16,1e18,cfg["at_1e16_gauge_GUT1_2_3"])]:
            end,ref,m0=run_gauge(inv_alpha(g),f,f,stop,hs,eps,beta)
            k=np.diag([ref[0],4*math.pi/cfg["horizontal_boundary_g"]**2])
            k=k-beta[-1]*math.log(stop/m0)/TWOPI
            assert math.isclose(k[0,0],end[0],rel_tol=1e-12,abs_tol=1e-12)
            examples.append({"f_GeV":f,"stop_GeV":stop,"all_mass_matching_GeV":float(m0),
                             "formal_inverse_alpha_at_end":end.tolist(),
                             "minimum_eigenvalue_abelian_inverse_kinetic_at_end":float(min(np.linalg.eigvalsh(k))),
                             "finite_control_at_end":bool(min(np.linalg.eigvalsh(k))>=1/cfg["alpha_control"]
                                                       and min(end[1:])>=1/cfg["alpha_control"])})
        numbers["finite_examples"]=examples
        assert not examples[0]["finite_control_at_end"] and examples[1]["finite_control_at_end"]
        output["epsilon_examples"].append(numbers)
    return output


if __name__=="__main__":
    p=ArgumentParser()
    group=p.add_mutually_exclusive_group()
    group.add_argument("--write-new",action="store_true")
    group.add_argument("--check",action="store_true")
    args=p.parse_args()
    data=run()
    if args.write_new:
        with (HERE/"results.json").open("x",encoding="utf8") as f:
            json.dump(data,f,indent=2)
            f.write("\n")
    else:
        assert data==fetch("results.json")
    print("1452 joint one-loop gauge and tree mass checks passed")
