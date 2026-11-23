"""Three-family cyclic selector built from charged spectral lines and a named neutral reference."""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
OLD=runpy.run_path(str(HERE.parent/"1744"/"check.py"))
OMEGA=np.exp(2j*math.pi/3)
FOURIER=OMEGA**np.outer(np.arange(3),np.arange(3))/math.sqrt(3)


def cycle(ce,ref):
    a=ce.conj().T@ref
    assert min(abs(a))/np.linalg.norm(a) > 1e-8
    w=ce*(a/abs(a))
    S=w@np.roll(np.eye(3),1,axis=0)@w.conj().T
    q=w@FOURIER
    np.testing.assert_allclose(q.conj().T@S@q,np.diag(OMEGA**(-np.arange(3))),atol=2e-15)
    return q


def pull(x,window):
    m,d,u=window
    return (x-m)/(d if x<=m else u)


def full(b,f,offsets):
    charged=b["model"]["JOINT"]["LEPTONS"]
    common=1.
    for _ in range(3):
        bs=dict(b)
        bs["yW"]=b["yW"].copy()
        bs["yW"][2]*=common
        row=OLD["daughter"](bs,f)
        actual=np.array(row["MZ_leptons"])
        if abs(actual[2]/charged[2,0]-1)<2e-12:
            break
        common*=charged[2,0]/actual[2]
    o=row["neutrinos"]
    ap=[pull(x,w) for x,w in zip(o["angle_readings"],offsets)]
    freq=b["nu"]["FREQ"]
    dp=[pull(o[k],freq[w]) for k,w in (("dm31_eV2","juno_NO_dm31"),("dm32_eV2","accelerator_NO_dm32"))]
    cp=(actual-charged[:,0])/charged[:,1]
    limits=freq["juno_NO_dm21"]
    allowed=np.array([limits[0]-2*limits[1],limits[0]+2*limits[2]])
    for k,w in (("dm31_eV2","juno_NO_dm31"),("dm32_eV2","accelerator_NO_dm32")):
        win=freq[w]
        lo,hi=np.array([win[0]-2*win[1],win[0]+2*win[2]])*limits[0]/o[k]
        allowed=[max(allowed[0],lo),min(allowed[1],hi)]
    return dict(charged_common=common,charged_pulls=cp.tolist(),angles=o["angle_readings"].tolist(),
                angle_pulls=ap,frequency_pulls=dp,dm31=o["dm31_eV2"],dm32=o["dm32_eV2"],
                allowed_solar_two_width_intersection=list(allowed) if allowed[0]<=allowed[1] else None,
                fc90=o["SK_T2K_FC68_90"][1],delta_deg=o["delta_deg"],mbb_meV=o["mbb_meV"],
                quark_four_ratio_radius=row["absolute_four_ratio_radius"],
                charged_clock_angles_pass=max(abs(np.r_[cp,ap,dp]))<=2 and bool(o["SK_T2K_FC68_90"][1]))


def solve(step=.125):
    b=OLD["mother"](step)
    nu=b["nu"]
    spectra=np.array([b["cW"],b["fW"],1.])
    ce,sv,_=np.linalg.svd(b["yW"][2,:,:3]); ce=ce[:,::-1]; sv=sv[::-1]
    raw=b["N"]
    u=np.linalg.svd(raw,full_matrices=True)[0]
    refs=dict(lower_column=raw[:,0],upper_column=raw[:,1],dominant_dirac=u[:,0],
              subleading_dirac=u[:,1],dark_line=u[:,2])
    pt=b["pt"]
    tr=np.linalg.inv(b["S"]@pt["low_E"])
    mass=tr.T@(pt["C"]+nu["FROZEN"][0]*np.outer(pt["b"],pt["b"]))@tr
    eigen=np.linalg.svd(mass)[2].conj().T
    refs["majorana_heavy"]=eigen[:,0]
    refs["majorana_light_massive"]=eigen[:,1]
    masses=np.asarray(pt["flow"]["fields"]["neutral_poles_GeV"])[[0,2]]
    kinetic=(raw/masses)@(raw/masses).conj().T
    for i,name in enumerate(("e","mu","tau")):
        refs["kinetic_probe_"+name]=kinetic@ce[:,i]
    sol=json.loads((HERE.parent/"1484"/"inputs.json").read_bytes())
    atm=json.loads((HERE.parent/"1659"/"inputs.json").read_bytes())["atmospheric_NOvA_upper_octant"]
    offsets=[sol["DayaBay_sin2_2theta13"],sol["JUNO_sin2_theta12_NO"],[atm["central"],atm["low"],atm["high"]]]
    # Exact spectrum follows from the two traces and determinant; all are phase/order independent.
    normalized=sv/sv[-1]
    A=np.diag(normalized**2)
    vals=spectra**2
    trace=np.trace(A)*np.mean(vals)
    elementary=np.linalg.det(A)*np.prod(vals)*sum(normalized**-2)*np.mean(1/vals)
    determinant=np.linalg.det(A)*np.prod(vals)
    sq=np.sort(np.roots([1,-trace,elementary,-determinant]))
    charged_ratio=np.sqrt(sq)/normalized
    cases={}
    invariant_abs_gram=None
    for name,ref in refs.items():
        q=cycle(ce,ref)
        for perm in itertools.permutations(range(3)):
            key=name+"_"+"".join(map(str,perm))
            f=(q*spectra[list(perm)])@q.conj().T
            new=np.linalg.svd(f@b["yW"][2,:,:3],compute_uv=False)[::-1]/sv
            np.testing.assert_allclose(new,charged_ratio,rtol=5e-12)
            left=np.linalg.svd(f@b["yW"][2,:,:3])[0][:,::-1]
            inv=np.linalg.inv(f)
            heavy=abs(left.conj().T@(inv@inv-np.eye(3))@left)
            if invariant_abs_gram is None:
                invariant_abs_gram=heavy
            else:
                np.testing.assert_allclose(invariant_abs_gram,heavy,atol=1e-13)
            cases[key]=full(b,f,offsets)
    return dict(scale_GeV=b["Mweak"],legs=spectra.tolist(),all_cyclic_charged_gains=charged_ratio.tolist(),
                tau_scaled_threshold_relative_shifts=(charged_ratio/charged_ratio[-1]-1).tolist(),
                charged_daughter_normalized_heavy_abs_gram=invariant_abs_gram.tolist(),
                cases=cases)


def simple(obj):
    if isinstance(obj,dict): return {k:simple(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple,np.ndarray)): return [simple(t) for t in obj]
    if isinstance(obj,(np.floating,float)): return float(f"{obj:.10g}")
    if isinstance(obj,np.bool_): return bool(obj)
    return obj


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    arg=ap.parse_args()
    out=HERE/"independent_results.json"
    ans=simple(solve())
    if arg.check:
        assert ans==json.loads(out.read_bytes())
        print("1748 cyclic flavour frames, exact charged polynomial and low-energy windows verified")
    else:
        with out.open("x",encoding="utf8") as f:
            json.dump(ans,f,indent=2,allow_nan=False)
            f.write("\n")
