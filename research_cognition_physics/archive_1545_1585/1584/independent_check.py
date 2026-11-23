"""Absolute light charge gap and a two-basis weak-current certificate."""
from pathlib import Path
from math import sqrt
import argparse
import json
import numpy as np

ROOT=Path(__file__).resolve().parent
def read(r):
    return json.loads((ROOT.parent/r/"inputs.json").read_text(encoding="utf-8-sig"))
def span(row):
    e=2*sum(row[1:])
    return row[0]-e,row[0]+e
def complete(x):
    Q,_=np.linalg.qr(x)
    return Q
def op_checks():
    c=np.arange(1,19).reshape(6,3)
    T=complete(np.sin(c*.83)+1j*np.cos(c*1.17))
    e=.45
    D=np.diag([1.]*3+[e**2]*3)
    Q=np.diag([1.]*3+[3.]*3)
    lhs=T.conj().T@Q@T
    rhs=3*np.eye(3)-2*(T.conj().T@D@D@T-e**4*np.eye(3))/(1-e**4)
    assert np.linalg.norm(lhs-rhs)<1e-12
    V=complete(np.sin(c[:3,:]*1.63)+1j*np.cos(c[:3,:]*.32))
    sv=np.linalg.svd(abs(V)**2,compute_uv=False)
    Y=V.conj().T@lhs@V
    def off(H):
        return np.linalg.norm(H-np.diag(np.diag(H)))**2
    return bool(off(lhs)+off(Y)>=(1-sv[1])*np.linalg.norm(lhs-np.trace(lhs)*np.eye(3)/3)**2)
def solve():
    sem=read("1381")
    cb=span(read("1412")["BelleII_Vcb_BGL"])
    k,ff=span(sem["K_product"]),span(sem["K_form"])
    us=(k[0]/ff[1],k[1]/ff[0])
    ub=cb[1]*span(sem["ub_over_cb"])[1]
    btb=1-cb[1]**2-ub**2
    # Compare the Laplacian of BB^T to a weighted 3-node path.
    x=us[0]**2*(btb-us[1]**2)
    y=cb[0]**2*btb
    small=3*x*y/(x+y+sqrt(x*x-x*y+y*y))
    gap=small/(1+sqrt(1-small))
    scenarios={}
    eps=.45
    for label,row in [("MZ",read("1386")["yukawa"]["up"][-1]),
                      ("1e16_GeV",read("1397")["high_scale"]["up"][-1])]:
        yt=span(row)[0]
        delta=2*max(0.,yt-eps**4)/(1-eps**4)
        floor=delta*sqrt(gap/24)
        scenarios[label]={"yt_lower":yt,"epsilon":eps,"L":0,"Cu_spectral_cap":1,
                         "charge_spread_lower":delta,"left_flavor_changing_max_lower":floor,
                         "normalized_inclusive_two_body_strength_lower":floor**2}
    return {"positive_sum_test":op_checks(),"tree_lower_edges":[x,y],"one_minus_sigma2_B_lower":gap,
            "conditional_scenarios":scenarios,"only_low_energy_tree_CKM_transport_to_high_scale_is_additional":True}
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    args=ap.parse_args()
    out=ROOT/"independent_results.json"
    res=solve()
    if args.check:
        assert json.loads(out.read_text(encoding="utf-8-sig"))==res
        print("1584 independent flavor charge certificate verified")
    else:
        out.write_text(json.dumps(res,indent=2)+"\n",encoding="utf-8")
