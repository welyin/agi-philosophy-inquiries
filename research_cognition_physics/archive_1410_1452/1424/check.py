"""Flavour-weighted disappearance frequencies and a finite matter-path first-dip check."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE=Path(__file__).resolve().parent
PHASE=1.26693268


def row_parameters(cfg, old, ordering):
    q=old[ordering]
    solar=cfg["solar_IO_fallback"] if ordering=="IO" else {
        "s12_squared":q["solar_sine_squared"],"dm21_eV2":q["solar_gap"]}
    return [solar["s12_squared"], q["s23_squared"] if "s23_squared" in q else q["atmospheric_sine_squared"],
            old["DayaBay_amplitude"],solar["dm21_eV2"],old["quark_sines_and_phase"][-1]]


def unpack(rows, width, sign):
    out=[]
    for row, side in zip(rows, sign):
        err=row[2] if len(row)==3 and side>0 else row[1]
        out.append(row[0]+side*width*err)
    return out


def difference(values, cos_sign):
    x,y,amplitude,dm,phase=values
    z=(1-math.sqrt(1-amplitude))/2
    u=math.sqrt(z*y/(1-y))
    return dm*(1-2*x-2*math.sqrt(x*(1-x))*u*cos_sign*math.cos(phase))/(1+u*u)


def branch_intervals(rows, width):
    vals=[unpack(rows,width,s) for s in itertools.product([-1,1],repeat=5)]
    xlo,ylo,alo,dmlo,phlo=np.min(vals,axis=0)
    xhi,yhi,ahi,dmhi,phhi=np.max(vals,axis=0)
    umax=math.sqrt((1-math.sqrt(1-ahi))/2*yhi/(1-yhi))
    # Ensures monotonicity in u for cos<0 and in x for both signs.
    assert 2*math.sqrt(xlo*(1-xlo))*math.cos(phhi)*(1-umax**2)>2*umax*(1-2*xlo)
    assert umax*math.cos(phlo)*(1-2*xlo)/math.sqrt(xlo*(1-xlo))<2
    intervals=[]
    for cos_sign in [1,-1]:
        a=[difference(v,cos_sign) for v in vals]
        intervals.append([min(a),max(a)])
    gap=intervals[1][0]-intervals[0][1]
    return intervals,gap


def matrix(values, cos_sign):
    x,y,a,dm,alpha=values
    z=(1-math.sqrt(1-a))/2
    delta=-alpha if cos_sign==1 else math.pi+alpha
    u12,u13,u23=np.eye(3,dtype=complex),np.eye(3,dtype=complex),np.eye(3,dtype=complex)
    for obj,i,j,ss,phase in [(u12,0,1,x,0),(u13,0,2,z,delta),(u23,1,2,y,0)]:
        obj[i,i]=obj[j,j]=math.sqrt(1-ss)
        obj[i,j]=math.sqrt(ss)*np.exp(-1j*phase)
        obj[j,i]=-obj[i,j].conj()
    return u23@u13@u12


def dip_frequency(H_vac, potential, path, charges):
    low,high=path["dip_window_GeV"]
    length=path["baseline_km"]
    def avg(energy):
        p=[]
        for charge in charges:
            h=H_vac.copy() if charge>0 else H_vac.conj().copy()
            h[0,0]+=charge*2e9*energy*potential
            eigen,vec=np.linalg.eigh(h)
            amp=np.dot(abs(vec[1,:])**2,np.exp(-2j*PHASE*length*eigen/energy))
            p.append(abs(amp)**2)
        return sum(p)/len(p)
    mesh=np.linspace(low,high,81)
    seed=np.argmin([avg(e) for e in mesh])
    assert 0<seed<len(mesh)-1
    left,right=mesh[seed-1],mesh[seed+1]
    for _ in range(52):
        e1,e2=(2*left+right)/3,(left+2*right)/3
        if avg(e1)<avg(e2):
            right=e2
        else:
            left=e1
    return math.pi*(left+right)/(4*PHASE*length)


def compute():
    cfg=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    old=json.loads((HERE/cfg["angles"]).read_text(encoding="utf-8"))
    ee,ee_err=cfg["daya_bay_ee_eV2"]
    path=cfg["matter_benchmark"]
    out={}
    for ordering,osign in [("NO",1),("IO",-1)]:
        rows=row_parameters(cfg,old,ordering)
        central=[p[0] for p in rows]
        boxes=[]
        for width in cfg["windows"]:
            intervals,gap=branch_intervals(rows,width)
            signed=intervals if osign==1 else [[-r,-l] for l,r in intervals]
            boxes.append({"width":width,"predicted_ee_minus_mumu_eV2":signed,"branch_gap_eV2":gap})
        lo,hi=1.0,2.0
        for _ in range(50):
            mid=(lo+hi)/2
            if branch_intervals(rows,mid)[1]>0:
                lo=mid
            else:
                hi=mid
        dips={}
        for p in [1,-1]:
            u=matrix(central,p)
            dm=central[3]
            m=np.array([0.,dm,osign*ee+central[0]*dm])
            h=(u*m)@u.conj().T
            f=abs(m[2]-abs(u[1,1])**2/(1-abs(u[1,2])**2)*dm)
            assert np.isclose(ee-f,osign*difference(central,p),atol=1e-16,rtol=1e-12)
            v=path["V_per_density_eV"]*path["rho_times_Ye_g_cm3"]
            for name,potential,charges in [
                ("vacuum",0,[1]),("matter_nu",v,[1]),("matter_antinu",v,[-1]),("matched_average",v,[1,-1])]:
                measured=dip_frequency(h,potential,path,charges)
                dips[f"cos_{p}_{name}_minus_moment_eV2"]=measured-f
        out[ordering]={"boxes":boxes,"box_width_closing_gap":lo,"first_dip_offsets":dips}
    ranges={}
    for key in ["MINOS_beam_two_flavour_eV2","MINOS_combined_two_flavour_eV2"]:
        m,dl,dh=cfg[key]
        ranges[key]={"central_difference":ee-m,"one_width":[ee-ee_err-(m+dh),ee+ee_err-(m-dl)],
                     "two_width":[ee-2*ee_err-(m+2*dh),ee+2*ee_err-(m-2*dl)]}
        assert ranges[key]["two_width"][0]<-4e-5 and ranges[key]["two_width"][1]>4e-5
    lim=cfg["finite_phase21_ceiling"]
    for r in np.linspace(0,1,25):
        for phase21 in np.linspace(.0003,lim,51):
            phasor=r*np.exp(2j*(1-r)*phase21)+(1-r)*np.exp(-2j*r*phase21)
            assert abs(np.angle(phasor))/2 <= phase21**3/(12*math.cos(2*lim))+1e-16
            assert np.isclose(abs(phasor)**2,1-4*r*(1-r)*math.sin(phase21)**2)
    phase_error=max(row_parameters(cfg,old,o)[3][0]+2*row_parameters(cfg,old,o)[3][1] for o in out)
    out["local_frequency_error_per_channel_bound_eV2"]=phase_error*lim**2/(12*math.cos(2*lim))
    out["heterogeneous_published_frequency_differences_eV2"]=ranges
    return out


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=json.dumps(compute(),indent=2,allow_nan=False)+"\n"
    dest=HERE/"results.json"
    if args.write:
        with dest.open("x",encoding="utf-8") as f: f.write(result)
    elif json.loads(result)!=json.loads(dest.read_text(encoding="utf-8")):
        raise SystemExit("1424 replay mismatch")
    print("1424 disappearance frequency intervals, finite-window phasor and matter dip calculation verified")
