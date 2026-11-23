"""When does a fixed number of antisymmetric exchange lines count participating families?"""
import argparse
import itertools
import json
import math
from pathlib import Path
import runpy
import numpy as np

HERE=Path(__file__).resolve().parent
CFG=json.loads((HERE/"inputs.json").read_bytes())
PORT=runpy.run_path(str(HERE/CFG["poles_source"]))


def response(h,i,k):
    return np.array([math.prod(h[i]-h[j] for j in subset) for subset in itertools.combinations(range(len(h)),k)])


def enclosures(h,i,k,kappa):
    vals=response(h,i,k)
    # For the unnormalized numerator the signs of the extrema determine all possible signs.
    neg=-sum(x for x in vals if x<0)
    pos=sum(x for x in vals if x>0)
    return pos-kappa*neg,kappa*pos-neg


def combinatorial_check():
    count=0
    for n in range(2,10):
        h=np.arange(n,dtype=float)**2+np.arange(n)*.5
        for k in range(1,n+1):
            for i in range(1,n-1):
                vals=response(h,i,k)
                nonzero=vals[vals!=0]
                assert bool(len(nonzero))==(k<=n-1)
                if k<=n-2:
                    assert min(nonzero)<0<max(nonzero)
                if k==n-1:
                    assert len(nonzero)==1 and np.sign(nonzero[0])==(-1)**(n-i-1)
                # Uniform exterior trace equals the coefficient of u^(n-k) in p(h_i+u).
                direct=np.polynomial.polynomial.polyfromroots(h-h[i])[n-k]
                np.testing.assert_allclose(sum(vals),direct,rtol=2e-13,atol=5e-8)
                count+=1
    return count


def family_tests(h):
    h=h/h[-1]
    a,b=np.diff(h)
    old=np.array([response(h,i,2).mean() for i in range(3)])
    assert list(np.sign(old))==[1,-1,1]
    heavy={}
    for spacing in CFG["counterexamples"]["fourth_heavier_from_h2_in_units_h3_minus_h2"]:
        extra=h[1]+spacing*b
        spectrum=np.r_[h,extra]
        vals=response(spectrum,1,2)
        threshold=b*(extra-h[1])/(a*(b+extra-h[1]))
        computed=sum(vals[vals>0])/-sum(vals[vals<0])
        np.testing.assert_allclose(threshold,computed,rtol=3e-14)
        for scale in (.999,1.001):
            lo,hi=enclosures(spectrum,1,2,scale*threshold)
            assert hi>0 and (lo<0)==(scale>1)
        heavy[str(spacing)]=dict(fourth_eigenvalue=extra,uniform_response=float(vals.mean()),
                                 anisotropy_at_zero=threshold)
    extra=(h[0]+h[1])/2
    hidden=np.r_[h,extra]
    seen=np.array([response(hidden,i,2).mean() for i in range(3)])
    assert list(np.sign(seen))==[1,-1,1]
    for d in (4,6,9):
        spectrum=list(h)+[h[-1]+offset for offset in range(1,d-2)]
        vals=response(spectrum,1,2)
        crit=sum(vals[vals>0])/-sum(vals[vals<0])
        assert crit>(d-3)*b/(2*a)
    return dict(first_three_squared_masses=h,three_family_uniform_responses=old,
                minimal_condition_number_for_any_fourth_heavier_than_third=b/(2*a),
                minimal_condition_number_if_next_multiple_is_six=3*b/(2*a),
                fourth_uniform_flip_if_crosses=h[1]+a*b/(b-a),
                extra_eigenvalue_inside_light_gap=extra,seen_uniform_responses_with_inside_extra=seen,
                heavy_tests=heavy)


def solve():
    bg=PORT["OLD"]["mother"](CFG["step"])
    charged=bg["yW"][2,:,:3]
    down=bg["yW"][1,:,:3]
    sectors={}
    for label,src in (("down",down),("charged_leptons",charged)):
        h=np.linalg.svd(src,compute_uv=False)[::-1]**2
        sectors[label]=family_tests(h)
    return dict(matching_scale_GeV=bg["Mweak"],combinatorial_cases=combinatorial_check(),sectors=sectors)


def native(x):
    if isinstance(x,dict):
        return {k:native(v) for k,v in x.items()}
    if isinstance(x,(list,np.ndarray)):
        return [native(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return float(f"{x:.12g}")
    return x


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--check",action="store_true")
    cli=p.parse_args()
    result=native(solve())
    path=HERE/"results.json"
    if cli.check:
        assert result==json.loads(path.read_bytes())
        print("1753 preparation-independent signs, bounded fourth-family response and dilution checks verified")
    else:
        with path.open("x",encoding="utf8") as out:
            json.dump(result,out,indent=2,allow_nan=False)
            out.write("\n")
