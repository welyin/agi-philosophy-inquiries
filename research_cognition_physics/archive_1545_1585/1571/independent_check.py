"""Exact weights of published orbifold generators and a conditional weak-triplet threshold."""
import argparse
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODELS = {
    "SU6_paper_eq11": [F(1,3)]*3 + [-F(1,3)]*3,
    "SU7_paper_eq43": [F(1,4)]*3 + [-F(1,4)]*2 + [F(1,4), -F(1,2)],
    "SU7_paper_eq48": [F(3,10)]*3 + [-F(2,5)]*2 + [F(1,10), -F(1,5)],
    "SU7_integer_direct_product": [F(x,6) for x in (-2,-2,-2,1,1,1,3)],
}


def winding(values):
    return math.lcm(*[v.denominator for v in values])


def model(y):
    n = len(y)
    assert sum(y) == 0 and len(set(y[:3])) == len(set(y[3:5])) == 1
    color_weak = [F(1,3),F(1,3),-F(2,3),F(1,2),-F(1,2)] + [F(0)]*(n-5)
    full = [a+b for a,b in zip(y,color_weak)]
    k = 2 * sum(x*x for x in y)
    z = winding(full)
    assert all((z*x).denominator == 1 for x in full)
    assert z*z*(F(1,3)+F(1,4)+k/4) == sum((z*x)**2 for x in full)/2
    assert (sum((z*x)**2 for x in full)/2).denominator == 1
    assert all(any((i*x).denominator != 1 for x in full) for i in range(1,z))
    pure = winding(y)
    # Root at slot 3 to 7 for paper SU7: the actual adjoint carries a nonintegral 6Y.
    noninteg = [str(abs(6*(a-b))) for i,a in enumerate(y) for b in y[i+1:]
                if (6*(a-b)).denominator != 1]
    return {"k":str(k), "fundamental_6Y":[str(6*x) for x in y],
            "nonintegral_adjoint_6Y":sorted(set(noninteg)),
            "pure_Y_integer_period":pure, "brane_Q_phase_turns_at_identity":str(F(pure,6)%1),
            "ordinary_diagonal_winding":z, "wound_half_norm":str(sum((z*x)**2 for x in full)/2)}


def threshold(k, count):
    inp=json.loads((HERE/"../1568/inputs.json").read_text(encoding="utf-8-sig"))
    A, a3 = inp["A_em_MSbar"][0], 1/inp["alpha_s_lattice"][0]
    by,bw,bc=inp["SM_beta_Y_2_3"]
    s_eff = .23098
    scheme = .00030  # example within the registered .001 transport, not an inferred correction
    s = s_eff - scheme
    # Directly solve three independent inverse couplings for aU, log(MX/MZ), log(MX/MT).
    assert k > 0 and count > 0
    t = (A*(1-s)-k*a3)/(by-k*bc)
    aU = a3 - bc*t
    trip_time = 3*(A*s-aU-bw*t)/count
    assert 0 < trip_time < t
    MX=inp["MZ_GeV"]*math.exp(2*math.pi*t)
    MT=MX*math.exp(-2*math.pi*trip_time)
    actual = [k*aU+by*t, aU+bw*t+count*trip_time/3, aU+bc*t]
    desired = [A*(1-s),A*s,a3]
    assert max(abs(x-y) for x,y in zip(actual,desired)) < 3e-14
    return {"effective_angle":s_eff, "effective_minus_MS_example":scheme, "MSbar_angle":s,
            "M_triplet_GeV":MT, "M_unification_GeV":MX, "inverse_alpha_U":aU,
            "real_triplet_count":count, "inverse_weak_shift":count*trip_time/3, "reconstructed_inverse_couplings":actual}


def run():
    data={key:model(y) for key,y in MODELS.items()}
    assert [v["ordinary_diagonal_winding"] for v in data.values()] == [6,12,30,6]
    assert [v["pure_Y_integer_period"] for v in data.values()] == [3,4,10,6]
    return {"weight_tests":data,"one_split_real_triplet":threshold(4/3,1), "canonical_SU5_three_adjoints":threshold(5/3,3)}


def rounded(x):
    if isinstance(x,dict):
        return {k:rounded(v) for k,v in x.items()}
    if isinstance(x,list):
        return [rounded(v) for v in x]
    return float(format(x,".13g")) if isinstance(x,float) else x


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    select=parser.add_mutually_exclusive_group(required=True)
    select.add_argument("--check", action="store_true")
    select.add_argument("--write-new", action="store_true")
    args=parser.parse_args()
    result=rounded(run())
    loc=HERE/"independent_results.json"
    if args.write_new:
        with loc.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    elif json.loads(loc.read_text(encoding="utf-8-sig")) != result:
        raise SystemExit("1571 independent result changed")
    print("1571: four compact embeddings and a weak-triplet threshold verified")
