"""Absolute heavy-to-light decay strengths for the six calibrated colour assignments."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


OLD = runpy.run_path(str(HERE/"../1766/independent_check.py"))
HEAVY = read("../1744/inputs.json")["pole_mass_GeV"][0]
GF = read("../1404/inputs.json")["const_GeV"]["GF"]
VEV = (math.sqrt(2)*GF)**-.5
MW, MZ, MH, MT = 80.379, 91.1876, 125.25, 172.5  # Fixed two-body phase-space convention, GeV.
HBARC_GEV_UM = 1.973269804e-10


def tensors(row, key):
    # Saved moduli are after transport; restore the three numerical shooting coordinates.
    target = np.array([row["moduli"][v] for v in ("Vus","Vub","Vcb")])
    trial = target.copy()
    for _ in range(4):
        made = OLD["eval_point"](np.log(trial),key,.25)
        diff = np.array([made["moduli"][v] for v in ("Vus","Vub","Vcb")])-target
        if np.max(abs(diff))<1e-12:
            break
        trial -= diff
    assert max(abs(diff))<1e-12
    bg,F,_ = OLD["load_mother"](*trial,.25)
    y=bg["yW"].copy()
    y[2]=made["tau_gain"]*F.conj().T@y[2]
    y[3]=F.conj().T@y[3]
    transport=np.eye(3,dtype=complex)
    for start,_,dt,grid in reversed(bg["segments"]):
        y,transport=bg["model"]["advance"](y,transport,-dt,grid[::-1],bg["pt"]["flow"])
        if start==HEAVY:
            up,yu,_=np.linalg.svd(y[0,:,:3])
            down,yd,_=np.linalg.svd(y[1,:,:3])
            V=up[:,::-1].conj().T@down[:,::-1]
            return yu[::-1],yd[::-1],V
    raise ValueError("missing colour threshold")


def vector_shape(x,y):
    # x,y are daughter fermion / boson squared masses in units of the parent squared mass.
    phase=math.sqrt((1-x-y)**2-4*x*y)
    return phase*((1-x)**2+y*(1+x)-2*y*y)


def higgs_shape(x,y):
    phase=math.sqrt((1-x-y)**2-4*x*y)
    return phase*((1+x)*(1+x-y)+4*x)


def response(a,f,ckm_column,up_masses):
    t=math.sqrt(1-f*f)
    mass=np.array([[a*f,a*t],[0.,HEAVY]])
    # Get the eigenvalue excess and mixing from the same 2x2 matrix without subtracting TeV poles.
    gap=math.hypot(HEAVY**2-a*a,2*HEAVY*a*t)
    excess=2*(HEAVY*a*t)**2/(gap+HEAVY**2-a*a)
    big=math.sqrt(HEAVY**2+excess)
    small=HEAVY*a*f/big
    eps=excess/gap
    light2=(small/big)**2
    common=big**3/(32*np.pi*VEV**2)
    W = 2*common*eps*abs(ckm_column)**2*np.array([
        vector_shape((u/big)**2,(MW/big)**2) for u in up_masses])
    Z=common*eps*(1-eps)*vector_shape(light2,(MZ/big)**2)
    h=common*eps*(1-eps)*higgs_shape(light2,(MH/big)**2)
    total=W.sum()+Z+h
    # Independent rotation of the whole Higgs source; one Higgs-dependent row only.
    left,s,right=np.linalg.svd(mass)
    np.testing.assert_allclose(s,[big,small],rtol=1e-12)
    Yuk=left.T@np.array([[a*f,a*t],[0.,0.]])@right.T/VEV
    np.testing.assert_allclose(np.abs([Yuk[1,0],Yuk[0,1]]),
                               np.sqrt(eps*(1-eps))*np.array([big,small])/VEV,rtol=4e-9,atol=2e-18)
    reduced=h/(common*higgs_shape(light2,(MH/big)**2))*(1-light2)**2/light2
    np.testing.assert_allclose(reduced, f**-2-1, rtol=5e-13)
    return dict(light_GeV=small,parent_a_GeV=a,heavy_excess_GeV=excess/(big+HEAVY),overlap=f,
                left_doublet_probability=eps,partial_W_u_c_t_GeV=W,partial_Z_GeV=Z,partial_h_GeV=h,
                total_three_modes_GeV=total,proper_length_three_modes_um=HBARC_GEV_UM/total,
                BR_W_Z_h=[sum(W)/total,Z/total,h/total],h_strength_over_light_mass=reduced,
                overlap_from_h=(1+reduced)**-.5)


def unresolved(data,indices,V,up_masses):
    """Evolution of an equal incoherent population with the common W states kept coherently."""
    a,b=data.values()
    widths=np.array([a["total_three_modes_GeV"],b["total_three_modes_GeV"]])
    shifts=np.array([a["heavy_excess_GeV"],b["heavy_excess_GeV"]])
    eps=np.array([a["left_doublet_probability"],b["left_doublet_probability"]])
    shape=[vector_shape((u/HEAVY)**2,(MW/HEAVY)**2) for u in up_masses]
    # Common-energy width matrix for overlapping poles. The variation of kinematics
    # between their individual evaluation energies is of order the relative mass shift.
    w_coeff=HEAVY**3/(16*np.pi*VEV**2)*math.sqrt(np.prod(eps))
    g=w_coeff*np.vdot(V[:,indices[0]], np.asarray(shape)*V[:,indices[1]])
    mixed=abs(g)**2/np.prod(widths)
    norm_change=mixed/(1+(2*(shifts[0]-shifts[1])/sum(widths))**2)
    diag=.5/(widths*(1-norm_change))
    off=-g*sum(diag)/(sum(widths)+2j*(shifts[0]-shifts[1]))
    occupation=np.array([[diag[0],off],[off.conjugate(),diag[1]]])
    G=np.array([[widths[0],g],[g.conjugate(),widths[1]]])
    H=np.diag(shifts)
    integrated=(G@occupation+occupation@G)/2+1j*(H@occupation-occupation@H)
    np.testing.assert_allclose(integrated, np.eye(2)/2, atol=2e-12)
    hprob=diag*np.array([a["partial_h_GeV"],b["partial_h_GeV"]])
    return dict(common_W_fractional_offdiagonal=abs(g)/math.sqrt(np.prod(widths)),
                common_W_normalization_change=norm_change/(1-norm_change),
                tagged_fractions=hprob/sum(hprob))


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):
        return [clean(v) for v in x]
    return float(f"{x:.10g}")


def calculate():
    archive=read("../1766/independent_results.json")
    results={}
    for source in ("refits","exploratory_64_deg"):
        entries={}
        for key,row in archive[source].items():
            yu,yd,V=tensors(row,key)
            c=abs(V[0,0])
            powers={"c":c,"h":np.sqrt((1+c)/2),"1":1.}
            up_masses=yu*VEV/np.sqrt(2)
            up_masses[-1]=MT
            data={}
            tags=("d","s","b")
            for i,symbol in enumerate(key):
                if symbol=="1":
                    continue
                data[tags[i]]=response(yd[i]*VEV/np.sqrt(2),powers[symbol],V[:,i],up_masses)
            assert len(data)==2
            a,b=data.values()
            split=abs(a["heavy_excess_GeV"]-b["heavy_excess_GeV"])
            ratio=sum(v["partial_h_GeV"]/(v["total_three_modes_GeV"]*2) for v in data.values())
            norms={tag:v["BR_W_Z_h"][2]/(2*ratio) for tag,v in data.items()}
            mixing=unresolved(data,[tags.index(tag) for tag in data],V,up_masses)
            np.testing.assert_allclose(mixing["tagged_fractions"],list(norms.values()),rtol=2e-14)
            bound=max(v["total_three_modes_GeV"] for v in data.values())
            assert bound<.2 and split<bound
            vals=sorted(v["h_strength_over_light_mass"] for v in data.values())
            np.testing.assert_allclose(vals[1]/vals[0], ((1+c)/c)**2,rtol=1e-12)
            entries[key]=dict(c=c,heavy_modes=data,conditional_flavour_given_h_in_equal_QCD_sample=norms,
                              ratio_of_two_reduced_h_strengths=vals[1]/vals[0],
                              overlapping_widths=mixing,
                              mass_split_GeV=split, mass_split_over_larger_width=split/bound)
        results[source]=entries
    return clean(dict(common_mass_GeV=HEAVY,vev_GeV=VEV, threshold_pair_partonic_GeV=2*HEAVY,
                      phase_space_mW_mZ_mh_mt_GeV=[MW,MZ,MH,MT],by_calibration=results))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    args=p.parse_args()
    result=calculate()
    if args.check:
        assert result==read("independent_results.json")
        print("1768 verified six colour assignments at both calibrations, three partial widths and inverse Higgs strengths")
    else:
        (HERE/"independent_results.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
