"""Independent response to a named opposite variation of the two charged sources."""
import argparse
import json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def read(name):
    return json.loads((ROOT/name).read_text(encoding="utf-8-sig"))

def run():
    cfg=read("inputs.json")
    frozen=read("results.json")
    h=read(cfg["high_quarks"])["high_scale"]
    low=read(cfg["low_spectra"])["yukawa"]
    lepton=np.array(read(cfg["high_leptons"])["scales"]["1e16"][3:])
    floors={}
    for label,ell in [("MZ",np.array(low["charged_leptons"])),("1e16",lepton)]:
        maxnorm=1+sum(((c+2*e)/(ell[-1,0]-2*ell[-1,1]))**2 for c,e in ell[:2])
        g=frozen[label]["sum_normalized_gain_floor"]
        floors[label]=g*g/maxnorm-1
    a,b,c,delta=np.array(h["angles_theta12_theta13_theta23_delta"])[:,0]
    x,y,z=np.sin([a,b,c]); X,Y,Z=np.cos([a,b,c])
    r12=np.array([[X,x,0],[-x,X,0],[0,0,1]])
    r23=np.array([[1,0,0],[0,Z,z],[0,-z,Z]])
    d=y*np.exp(-1j*delta)
    mix=r23@np.array([[Y,0,d],[0,1,0],[-d.conjugate(),0,Y]])@r12
    p1,p2,p3,p4=cfg["charged_only_completion_relative_phases"]
    up=np.diag(np.array(h["up"])[:,0]*np.exp(1j*np.array([p1,p2,0])))
    down=(mix*(np.array(h["down"])[:,0]*np.exp(1j*np.array([p3,p4,0]))))@mix.T
    gain,other=[complex(*v) for v in frozen["charged_only_completion"]["coeff_c_d_at_phi0"]]
    E=gain*down+other*up
    Q=gain*down-other*up
    l,ss,rh=np.linalg.svd(E)
    q=(l.conj().T@Q@rh.conj().T)[::-1,::-1]
    assert np.allclose(ss[::-1],lepton[:,0],rtol=1e-6,atol=1e-10)
    widths=abs(q)**2/np.linalg.norm(E)**2
    assert np.isclose(np.sum(widths),np.linalg.norm(Q)**2/np.linalg.norm(E)**2)
    assert widths.sum()>floors["1e16"]
    return {"strict_working_floor_massless_inclusive_relative_width":
            {k:round(v,8) for k,v in floors.items()},
            "actual_charged_witness_inclusive_relative_width":round(float(widths.sum()),8),
            "flavour_violating_two_orders_over_reference_inclusive":
            {key:round(float(widths[i,j]+widths[j,i]),8)
             for key,i,j in [("em",0,1),("et",0,2),("mt",1,2)]}}

if __name__=="__main__":
    p=argparse.ArgumentParser()
    f=p.add_mutually_exclusive_group()
    f.add_argument("--write-new",action="store_true");f.add_argument("--check",action="store_true")
    args=p.parse_args()
    data=run()
    if args.write_new:
        with (ROOT/"independent_results.json").open("x",encoding="utf8") as stream:
            json.dump(data,stream,indent=2);stream.write("\n")
    elif args.check:
        assert data==read("independent_results.json")
        print("1454 charged completion and relative-source width checked independently")
    else:
        print(json.dumps(data,indent=2))
