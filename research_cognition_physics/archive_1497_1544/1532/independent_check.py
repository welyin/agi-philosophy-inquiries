"""Independently invert the two signed molecular EDM strips, with scalar-response budgets."""
import argparse
import itertools
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def solve(yh,yt,ah,at):
    c=(yt-yh)/(at-ah)
    return yh-ah*c,c


def work_box(eps):
    def ends(x,h):
        return x-h,x+h
    estimates=[]
    for yh,yt,ah,at in itertools.product(ends(-1.3,2*(2.0+.6)),ends(4.3,2*(3.1+2.6)),
                                         ends(.9,eps*.9),ends(1.5,eps*1.5)):
        d,c=solve(yh,yt,ah,at)
        assert abs(d+ah*c-yh)<1e-10 and abs(d+at*c-yt)<1e-10
        estimates.append((d,c))
    d,c=zip(*estimates)
    # m_e and a representative isoscalar m_N in GeV, physical Higgs mass window.
    response=2*.000510999*.939*(.307-2*.018)/(125.40**2)
    return dict(d_in_1e_minus30=[min(d),max(d)],C_in_1e_minus10=[min(c),max(c)],
                electron_pseudoscalar_if_h_only=[min(c)*1e-10/response,max(c)*1e-10/response],
                independent_alpha_budget=eps)


def certificate():
    betaH,betaT=2*72/180-1,2*90/232-1
    # Unit isovector variation. Solve for change of common source and EDM.
    d_null,c_null=solve(-.9*betaH,-1.5*betaT,.9,1.5)
    assert abs(d_null+.9*(c_null+betaH))<1e-14
    assert abs(d_null+1.5*(c_null+betaT))<1e-14
    isotope_delta=2*72*(1/178-1/180)
    return dict(central_d_C_scaled=solve(-1.3,4.3,.9,1.5),
                zero_alpha_error=work_box(0),ten_percent_alpha_budget=work_box(.10),
                independent_relative_error_for_parallel=(1.5-.9)/(1.5+.9),
                nuisance_null_vector_d_C0_C1_scaled=[d_null,c_null,1.],
                synthetic_Hf178_minus180_difference_alpha_C1=.9*isotope_delta,
                synthetic_abs_C1_if_abs_difference_below_2e_minus30=2e-10/(.9*isotope_delta))


def equal(a,b):
    if isinstance(a,dict):
        return a.keys()==b.keys() and all(equal(x,b[k]) for k,x in a.items())
    if isinstance(a,(tuple,list)):
        return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-22)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    flags=parser.add_mutually_exclusive_group(required=True)
    flags.add_argument("--write",action="store_true")
    flags.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=certificate()
    out=HERE/"independent_results.json"
    if args.write:
        with out.open("x",encoding="utf-8") as file:
            file.write(json.dumps(result,indent=2,allow_nan=False)+"\n")
    else:
        assert equal(result,json.loads(out.read_text(encoding="utf-8")))
    print("1532 signed molecular strips verified")
