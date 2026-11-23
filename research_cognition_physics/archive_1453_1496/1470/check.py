"""Exact accepted-amplitude bounds and one raw Poisson count without background subtraction."""
import argparse
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def poisson_lower_tail(n, mean):
    term=math.exp(-mean)
    total=term
    for k in range(1,n+1):
        term*=mean/k
        total+=term
    return total


def poisson_upper(n,confidence):
    lo,hi=0.,2.*n+50.
    for _ in range(80):
        mid=(lo+hi)/2
        if poisson_lower_tail(n,mid)>1-confidence:
            lo=mid
        else:
            hi=mid
    return (lo+hi)/2


def lower(p,q,a,b):
    # a=norm(opposite), b=norm(same), with a>=b.
    return ((p+q)*a-abs(p-q)*b)/2


def evaluate():
    inp=json.loads((HERE/"inputs.json").read_text(encoding="utf-8"))
    count=inp["VBF_WH_same_SRminus"]
    obs=inp["associated_VHbb"]
    sw=inp["response_surpluses"]["WH"]
    sz=inp["response_surpluses"]["ZH"]
    plus,ep=count["unit_same_signal"]
    minus,em=count["unit_opposite_signal"]
    U=poisson_upper(count["raw"],inp["confidence_raw_only"])
    rows={}
    for n in [0,1,2]:
        a=math.sqrt(minus-n*em)
        b=math.sqrt(plus+n*ep)
        p=math.sqrt(obs["WH"][0]-n*obs["WH"][1]-sw)
        q=math.sqrt(obs["ZH"][0]-n*obs["ZH"][1]-sz)
        amp=lower(p,q,a,b)
        assert p>=q and amp>math.sqrt(U)
        # A one-dimensional real amplitude pair realizes the lower bound exactly.
        W=(a-b)/2
        Z=(-a-b)/2
        assert math.isclose(abs(p*W-q*Z),amp)
        assert math.isclose((W-Z)**2,minus-n*em) and math.isclose((W+Z)**2,plus+n*ep)
        qcrit=max(0.,(2*math.sqrt(U)-(a-b)*p)/(a+b))
        assert qcrit<=p
        deltaZcrit=obs["ZH"][0]-n*obs["ZH"][1]-qcrit**2
        rows[str(n)]={"min_effective_amplitudes":[p,q],"min_opposite_signal":amp**2,
                      "poisson_lower_tail_at_minimum":poisson_lower_tail(count["raw"],amp**2),
                      "deltaZ_critical_at_selected_deltaW":deltaZcrit,
                      "rate_transfer_required_to_keep_gap":U/amp**2,
                      "unseen_coherent_amplitude_needed":amp-math.sqrt(U),
                      "relative_amplitude_needed":1-math.sqrt(U)/amp,
                      "minimum_destructive_overlap":(a*a-b*b)/(a*a+b*b)}
    a,b=math.sqrt(minus),math.sqrt(plus)
    return {"raw_signal_plus_background_upper":U,
            "same_opposite_yield_ratio":plus/minus,
            "relative_inner_product_twice_over_normsum":(plus-minus)/(plus+minus),
            "unknown_individual_template_norm_range":[(a-b)/2,(a+b)/2],
            "boxes":rows}


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write-new",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=evaluate()
    path=HERE/"results.json"
    if args.write_new:
        with path.open("x",encoding="utf-8") as out:
            json.dump(result,out,indent=2)
            out.write("\n")
    else:
        assert result==json.loads(path.read_text(encoding="utf-8"))
        print("1470 author OK")
