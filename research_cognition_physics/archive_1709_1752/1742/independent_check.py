"""Separate phase envelope, two-flavour/finite-mass bound and Wilson transport."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent
def read(path):
    return json.loads((HERE/path).read_bytes())
PI=math.pi


def flavor():
    src=runpy.run_path(str(HERE/"../1739/check.py"))
    x=src["CFG"]["nearby_measured_inputs_in_marginal_widths"]["absolute_and_oscillations"]
    _,pt,_=src["charged_endpoint"](2,.125,x["weak"],x["up"],x["g3"],x["top"],x["charm"])
    # In 1731.endpoint E has the T index on the first slot. E.T supplies the weak-doublet frame.
    def lines(Y):
        return np.linalg.eigh(Y@Y.conj().T)[1]
    e=lines(pt["E"])
    d=lines(pt["D"])[:,0]
    vector=np.array([2.,1.,1.])*d
    norm=np.linalg.norm(vector)
    pieces=abs(vector[:,None]*e.conj())
    lower=max(0,2*np.max(pieces[:,0])-sum(pieces[:,0]))
    tau_max=sum(pieces[:,2])
    assert abs(lower-1.92872203)<2e-6
    assert abs(norm**2-1-3*abs(d[0])**2)<1e-12
    # Unrestricted right-handed e,mu,tau rotation can send x into tau.
    free_frame=np.linalg.svd(vector[:,None],full_matrices=True)[0][:,[1,2,0]]
    hide=abs(free_frame.conj().T@vector)
    assert max(hide[:2])<1e-14
    return dict(lengths=pieces[:,0], norm=norm, min_e_single_diagonal=lower,
                max_tau_single_diagonal=tau_max,unmixed_emu_squared_floor=norm**2-tau_max**2,
                fixed_absolute_components=abs(e.conj().T@vector), unconstrained_absolute_components=hide)


def running():
    gauge=read("../1369/inputs.json")
    base=read("../1698/results.json")["shared_witness"]
    registered=read("../1739/gauge_results.json")["R1698_two_new_colours_and_soft_scalars_at_g3_minus_0p67"]
    mass=read("../1739/results.json")["rows"]["absolute_and_oscillations"]["threshold_GeV"]
    # Fix weak mass, refit crossing separately by interpolation between a1 and a2.
    L=np.array([.6,1.,0.]); D=np.array([.4,0.,1.])
    event=[]
    def add(ms,db):
        event.extend((m,db) for m in ms)
    add(base["colored_fermions_GeV"]+[mass]*2,(2/3)*D)
    add(base["colored_scalars_GeV"]+[math.hypot(mass,base["soft_GeV"])]*2,D/3)
    add(base["weak_fermions_GeV"],(2/3)*L)
    add(base["weak_scalars_GeV"],L/3)
    add([registered["three_equal_complete_weak_poles"]]*3,L)
    sm=np.array([4.1,-19/6,-7.])
    susy=np.array([6.6,1.,-3.])
    mz=gauge["mu_ref_GeV"]
    low=np.array([gauge[k][0] for k in ("g1_GUT","g2","g3")])
    low[2]-=.67*gauge["g3"][1]
    inv=4*PI/low**2
    quad_x,quad_w=np.polynomial.legendre.leggauss(24)
    grid=sorted(set([mz,base["soft_GeV"]]+[m for m,_ in event]))
    # Final interval has constant group beta. Solve a1=a2 in that interval.
    logs=np.zeros(3); quad=0.; pieces={}
    prev=mz
    for i,start in enumerate(grid):
        assert start==prev
        if start==base["soft_GeV"]:
            inv-=np.array([0.,2.,3.])/(12*PI)
        is_susy=start>=base["soft_GeV"]
        beta=(susy if is_susy else sm).copy()
        for m,step in event:
            if m<=start:
                beta+=step
        span=math.log(grid[i+1]/start) if i+1<len(grid) else 2*PI*(inv[0]-inv[1])/(beta[0]-beta[1])
        gm=np.array([23/30,3/2,4/3]) if is_susy else np.array([23/20,9/4,2])
        end=inv-beta*span/(2*PI)
        for j,b in enumerate(beta):
            logs[j]+=gm[j]*(math.log(inv[j]/end[j])/b if abs(b)>1e-12 else span/(2*PI*inv[j]))
        samples=inv[None,:]-(span*(quad_x+1)/2)[:,None]*beta[None,:]/(2*PI)
        quad+=span/2*(quad_w@(1/samples)@gm)/(2*PI)
        inv=end
        prev=grid[i+1] if i+1<len(grid) else start*math.exp(span)
        if start==mz or prev==base["soft_GeV"] or i+1==len(grid):
            pieces[str(prev)]=math.exp(sum(logs))
    assert abs(quad-sum(logs))<1e-12 and max(inv)-min(inv)<2e-8
    assert abs(prev/registered["U"]-1)<3e-9
    return dict(U=prev,alpha5_inverse=np.mean(inv), short=math.exp(quad), by_group=np.exp(logs),
                cumulative=pieces,quad_log_difference=quad-sum(logs))


def limits(family,transport):
    mp,mpi=.938272,.134977
    sec=365.25*86400
    hbar=6.582119569e-25
    u,alpha=transport["U"],transport["alpha5_inverse"]
    L=transport["short"]*1.25
    g2=4*PI/alpha
    base_W=[.159,.161]
    form1=[.169,.170]
    rows={}
    for mode,min_life,ml,n in (("electron",2.4e34,.000510999,0),("muon",1.6e34,.1056584,1)):
        E=(mp*mp+ml*ml-mpi*mpi)/(2*mp)
        momentum=math.sqrt(E*E-ml*ml)
        x=ml/E
        output={}
        for name,W,W1 in (("central",base_W[n],form1[n]),("rectangle",base_W[n]-sum(([.015,.020,.025],[.015,.020,.026])[n]),
                                  form1[n]+sum(([.014,.018,.029],[.014,.018,.028])[n]))):
            # Isospin relation for both (parity-related) form factors.
            W,W1=W/math.sqrt(2),W1/math.sqrt(2)
            b=ml/mp*W1/W
            schur=(1-b*b)**2*(1-x*x)/(1+b*b+2*x*b)
            Q=np.array([[1.,-b],[-b,1.]])
            G=Q@np.array([[1.,-x],[-x,1.]])@Q
            small=G[0,0]-G[0,1]**2/G[1,1]
            assert abs(small-schur)<1e-12 and b<1
            kin=E*momentum/(8*PI*mp)
            pre=kin*W**2*schur
            maxf=math.sqrt(hbar/(sec*min_life*pre))*u*u/(L*g2)
            lifetime=hbar/(sec*pre*(family["min_e_single_diagonal"]*L*g2/u**2)**2)
            output[name]=dict(W0_GeV2=W,W1_GeV2_upper=W1,helicity_worst_factor=schur,prefactor=pre,
                    allowed_abs_f_at_crossing=maxf, tagged_e_lifetime_years_if_electron_channel=lifetime if n==0 else None,
                    conditional_min_MX_over_U_if_electron_channel=math.sqrt(family["min_e_single_diagonal"]/maxf) if n==0 else None)
        rows[mode]=output
    both={}
    for name in ("central","rectangle"):
        cap=sum(rows[tag][name]["allowed_abs_f_at_crossing"]**2 for tag in ("electron","muon"))
        norm2=family["norm"]**2
        min_tau=max(0.,1-cap/norm2)
        # If e and mu may be rotated arbitrarily with the old tau row left intact,
        # at least one two-sided working limit is exceeded for all allowed phases.
        severity=family["unmixed_emu_squared_floor"]/cap
        both[name]=dict(f_e_squared_plus_f_mu_squared_allowed=cap, tau_fraction_required_at_U=min_tau,
                        at_least_one_rate_over_its_limit_with_tau_fixed=severity,
                        necessary_frozen_factor_MX_over_U_with_tau_fixed=severity**.25)
    return dict(modes=rows, simultaneous=both)


def run():
    f=flavor()
    r=running()
    data=limits(f,r)
    assert data["simultaneous"]["rectangle"]["at_least_one_rate_over_its_limit_with_tau_fixed"]>1
    # Varying a long matching multiplies a lifetime by L^-2 and the inferred mass by sqrt(L).
    h=data["modes"]["electron"]["rectangle"]["conditional_min_MX_over_U_if_electron_channel"]
    hf=data["modes"]["electron"]["central"]["conditional_min_MX_over_U_if_electron_channel"]
    reference=data["modes"]["electron"]["central"]["W0_GeV2"]
    # Leading me=0 cross-check against the older pion matrix element (already for pi0).
    rate0=.938272*(1-(.134977/.938272)**2)**2/(32*PI)
    corr=(rate0*reference**2/data["modes"]["electron"]["central"]["prefactor"])**.25
    out=dict(flavor=f,transport=r,bounds=data,
             A_long_max_for_crossing_tagged_at_rectangle=1.25/h**2,
             A_long_max_for_crossing_with_tau_fixed_at_rectangle=1.25/math.sqrt(data["simultaneous"]["rectangle"]["at_least_one_rate_over_its_limit_with_tau_fixed"]),
             alternative_older_W0_leading_electron_mass_ratio_central_and_low=
                 [hf*math.sqrt(val/reference)*corr for val in (.131,.131-.004-.013)],
             mass_ratio_using_alternative_long_factor_1p247=h*math.sqrt(1.247/1.25))
    def tidy(x):
        if isinstance(x,dict): return {k:tidy(v) for k,v in x.items()}
        if isinstance(x,(np.ndarray,list,tuple)): return [tidy(v) for v in x]
        if isinstance(x,(float,np.floating)): return float(f"{x:.11g}")
        return x
    return tidy(out)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.check:
        assert result==read("independent_results.json")
        print("1742 independent finite-mass e/mu, flavor norm and gauge evolution passed")
    else:
        (HERE/"independent_results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf8")

