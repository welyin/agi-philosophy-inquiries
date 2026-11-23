"""Cross-check by diagonalizing the two dephasing projections on traceless Hermitians."""
import argparse
import json
from math import sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_text(encoding="utf-8"))


def generators():
    result=[np.diag([1.,-1.,0.])/sqrt(2), np.diag([1.,1.,-2.])/sqrt(6)]
    for i in range(3):
        for j in range(i):
            a=np.zeros((3,3),complex)
            a[i,j]=1
            result.extend([(a+a.T)/sqrt(2),1j*(a-a.T)/sqrt(2)])
    return np.array(result)


def mixing(angles):
    result=np.eye(3,dtype=complex)
    for i,j,angle,delta in [(1,2,angles[2],0),(0,2,angles[1],angles[3]),(0,1,angles[0],0)]:
        c,s=np.cos(angle),np.sin(angle)
        rot=np.eye(3,dtype=complex)
        rot[i,i]=rot[j,j]=c
        rot[i,j]=s*np.exp(-1j*delta)
        rot[j,i]=-rot[i,j].conjugate()
        result=result@rot
    return result


def run():
    lo=read("../1386/inputs.json")
    hi=read("../1397/inputs.json")["high_scale"]
    author=read("results.json")["same_scale"]
    frames=generators()
    metrics={}
    for scale,angles,masses in [
        ("MZ",lo["quark_sines_and_phase"],lo["yukawa"]),
        ("GUT_SM",hi["angles_theta12_theta13_theta23_delta"],hi)]:
        a=np.array(angles)[:,0]
        if scale=="MZ":
            a[:3]=np.arcsin(a[:3])
        v=mixing(a)
        design=[]
        for frame in (frames,v.conj().T@frames@v):
            off=frame.copy()
            off[:,np.arange(3),np.arange(3)]=0
            design.append(off.reshape(8,-1))
        gram=sum((o.conj()@o.T).real for o in design)
        exact=float(np.linalg.eigvalsh(gram)[0])
        assert np.isclose(exact,author[scale]["projection_gap_central"])
        chi=0
        for sector in ["up","down"]:
            center,err=np.array(masses[sector]).T
            chi=max(chi,max((center[:-1]+2*err[:-1])/(center[1:]-2*err[1:])))
        contrast=(1-chi**2)/(1+chi**2)
        floor=author[scale]["at_least_one_reconstructed_pair_magnitude_floor"]*contrast
        metrics[scale]={"central_direct_Hermitian_gap":round(exact,12),
                        "two_width_max_light_over_heavy":round(chi,12),
                        "spin_summed_normalized_amplitude_floor":round(floor,12)}
    rng=np.random.default_rng(431)
    for i in range(75):
        m,n=sorted(rng.uniform(.01,1,size=2))
        left,right=rng.normal(size=2)+1j*rng.normal(size=2)
        for test in [right,-(2*m*n)/(m*m+n*n)*left]:
            x,y=n*left+m*test,m*left+n*test
            power=abs(x)**2+abs(y)**2
            lower=(n*n-m*m)**2/(n*n+m*m)*abs(left)**2
            assert power+1e-12>=lower
            if test!=right:
                assert abs(power-lower)<1e-12
    return metrics


if __name__=="__main__":
    flags=argparse.ArgumentParser()
    opt=flags.add_mutually_exclusive_group()
    opt.add_argument("--write-new",action="store_true")
    opt.add_argument("--check",action="store_true")
    if flags.parse_args().write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf-8") as out:
            json.dump(run(),out,indent=2)
            out.write("\n")
    else:
        assert run()==read("independent_results.json")
        print("1431 independent projection gap and spin-summed response checked")
