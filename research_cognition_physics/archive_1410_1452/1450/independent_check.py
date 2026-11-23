"""Independent zero-det-index contraction and sharp domain-number-one witnesses."""
import argparse
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent

def read(name):
    return json.loads((HERE/name).read_text(encoding="utf-8-sig"))

def ckm(triples):
    a,b,c,phase=triples
    sina,sinb,sinc=np.sin([a,b,c])
    cosa,cosb,cosc=np.cos([a,b,c])
    r12=np.array([[cosa,sina,0],[-sina,cosa,0],[0,0,1]])
    r23=np.array([[1,0,0],[0,cosc,sinc],[0,-sinc,cosc]])
    v=sinb*np.exp(-1j*phase)
    r13=np.array([[cosb,0,v],[0,1,0],[-v.conjugate(),0,cosb]])
    return r23@r13@r12

def check_contraction():
    rng=np.random.default_rng(1450)
    for charges in [[-2,0,3],[-1,-1,2],[3,0,3]]:
        q=np.array(charges)
        diff=q[:,None]-q[None,:]
        for eps in (.08,.22,.7):
            core=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
            core[diff<0]=0
            core+=4*np.eye(3)
            gram=eps**abs(diff)
            assert min(np.linalg.eigvalsh(gram)) > -1e-12
            mass=gram*core
            assert np.allclose(np.linalg.inv(mass),gram*np.linalg.inv(core))
            for mat in [core,np.linalg.inv(core)]:
                assert np.linalg.norm(gram*mat,ord=2)<=np.linalg.norm(mat,ord=2)+1e-12

def run():
    config=read("inputs.json")
    mz=read(config["low_masses"])
    high=read(config["high_quarks"])["high_scale"]
    low_ckm=np.array(mz["quark_sines_and_phase"])[:,0]
    vlow=ckm([*np.arcsin(low_ckm[:3]),low_ckm[3]])
    vhigh=ckm(np.array(high["angles_theta12_theta13_theta23_delta"])[:,0])
    lepton_high=read(config["high_charged_leptons"])["scales"]["1e16"][3:]
    groups={"MZ":([mz["yukawa"][k] for k in ["up","down","charged_leptons"]],vlow),
            "1e16":([high["up"],high["down"],lepton_high],vhigh)}
    results={}
    for label,(sector,rotation) in groups.items():
        rows=[np.array(block) for block in sector]
        widest=[row[0,0]+2*row[0,1] for row in rows[:2]]
        target=1/max(widest)
        spectra=[row[:,0].copy() for row in rows]
        spectra[1][0]=widest[1]
        eps=.22
        charges=[np.array([1,0,0]),np.zeros(3,dtype=int),np.array([2,0,0])]
        smallest=2.
        largest=0.
        for yy,x,frame in zip(spectra,charges,[np.eye(3),rotation,np.eye(3)]):
            core=frame@np.diag(yy/eps**x)
            masses=core@np.diag(eps**x)
            assert np.allclose(masses,frame@np.diag(yy))
            sig=np.linalg.svd(core,compute_uv=False)
            smallest=min(smallest,sig[-1])
            largest=max(largest,sig[0])
        assert abs(smallest*target-1)<1e-12 and largest<=target
        assert sum(sum(q) for q in charges[:2])==1
        results[label]={"zero_up_necessary_B":round(1/widest[0],5),
                        "zero_down_necessary_B":round(1/widest[1],5),
                        "sharp_B_D_equal_one":round(target,5),
                        "example_core_singular_max":round(float(largest),8)}
    check_contraction()
    return results

if __name__=="__main__":
    parse=argparse.ArgumentParser()
    flags=parse.add_mutually_exclusive_group()
    flags.add_argument("--write-new",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parse.parse_args()
    out=run()
    if args.write_new:
        with (HERE/"independent_results.json").open("x",encoding="utf8") as f:
            json.dump(out,f,indent=2)
            f.write("\n")
    elif args.check:
        assert out==read("independent_results.json")
        print("1450 sharp mass-resolved domain-number threshold verified")
    else:
        print(json.dumps(out,indent=2))
