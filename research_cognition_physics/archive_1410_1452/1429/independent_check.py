"""1429. Alternative mediator: reciprocal lepton gains on the two off-diagonal ports."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def compute():
    data=json.loads((HERE/"../1367/inputs.json").read_text(encoding="utf-8"))
    values=data["points"]["GUT_SM"]
    def solve(m):
        d,s,e,mu=m
        D,E=s/d,mu/e
        pd,pe=D+1/D,E+1/E
        assert pe>9*pd
        k=E/D
        eta=.5*math.log(k)
        # Independent trace/determinant reconstruction of the ordinary h-dependent ceiling.
        old_pe=2+9*(pd-2)
        old_re=(old_pe+math.sqrt(old_pe**2-4))/2
        h=pe-9*(pd-2)
        k_bal=(h+math.sqrt(h*h-4))/2
        return dict(pe_minus_9pd=pe-9*pd,minimum_gain_ratio=k,
                    minimum_log_amplitude_displacement=eta,upper_gain_at_optimum=math.exp(eta),
                    minimum_gain_ratio_if_quark_links_balanced=k_bal,
                    electron_floor_ordinary=mu/old_re)
    results={}
    for n in (0,2):
        endpoints=[[(v-n*w)*u,(v+n*w)*u] for v,w,u in
                   zip(values["central"],values["half_width"],data["units"])]
        samples=[solve(m) for m in itertools.product(*endpoints)]
        results[str(n)]={k:[min(row[k] for row in samples),max(row[k] for row in samples)] for k in samples[0]}
    # Reach both central optimal endpoints using matrix singular values, without importing author code.
    d,s,e,mu=[v*u for v,u in zip(values["central"],data["units"])]
    D,E=s/d,mu/e
    for q,k in [(math.log(D),results["0"]["minimum_gain_ratio"][0]),
                (0.,results["0"]["minimum_gain_ratio_if_quark_links_balanced"][0])]:
        h=max(D+1/D-2*math.cosh(q),0)
        a,b,c=math.exp(q/2),math.exp(-q/2),math.sqrt(h)
        quark=np.array([[0.,a],[b,c]])
        lepton=np.array([[0.,a*math.sqrt(k)],[b/math.sqrt(k),-3*c]]).T
        for M,ratio in ((quark,D),(lepton,E)):
            u,v=np.linalg.svd(M,compute_uv=False)
            assert math.isclose(u/v,ratio,rel_tol=2e-12)
    return results


def near(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(near(a[k],b[k]) for k in a)
    if isinstance(a,list):
        return len(a)==len(b) and all(near(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=2e-11,abs_tol=1e-15)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    opts=p.add_mutually_exclusive_group()
    opts.add_argument("--write",action="store_true")
    opts.add_argument("--check",action="store_true")
    args=p.parse_args()
    values=compute()
    path=HERE/"independent_results.json"
    if args.write:
        with path.open("x",encoding="utf-8") as stream:
            json.dump(values,stream,indent=2)
            stream.write("\n")
    else:
        assert near(values,json.loads(path.read_text(encoding="utf-8")))
    print(json.dumps(values))
