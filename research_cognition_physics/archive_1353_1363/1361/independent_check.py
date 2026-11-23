"""Check the raised and lowered responses directly on six published moment intervals."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path


# Units are hundredths, independently transcribed from the source (table III).
Q2 = [126, 200, 316, 501, 794, 1259]
MEANS = [239, 249, 255, 278, 282, 280]
STAT = [8, 8, 6, 6, 7, 13]
SYST = [14, 10, 10, 19, 19, 18]


def witness(n, c):
    """Choose one common H and construct explicit responses and residuals."""
    h = F(3, 10)
    moments, factors = [], []
    for q, m, stat, syst in zip(Q2, MEANS, STAT, SYST):
        q2 = F(q, 100)
        # Even give all bins one common positive systematic displacement.
        datum = F(m + syst, 100)
        lo = max(datum - F(stat, 100), n*c - h/q2 - F(2, 100))
        hi = min(datum + F(stat, 100), n - h/q2 + F(2, 100))
        if lo > hi:
            return None
        moment = (lo+hi)/2
        coefficient = min(F(1), max(c, (moment+h/q2)/n))
        assert abs(moment - (n*coefficient - h/q2)) <= F(2, 100)
        moments.append(float(moment))
        factors.append(float(coefficient))
    return {"moments": moments, "responses": factors}


def run():
    bounds = [F(m+u+v+2, 100)+F(30,q) for q,m,u,v in zip(Q2,MEANS,STAT,SYST)]
    upper = min(bounds)
    assert Q2[bounds.index(upper)] == 316
    crit = {}
    for count in (3, 4, 6, 9):
        floor = upper/count
        assert witness(count,floor) is not None
        # For H <= .30, even every datum at its upper edge cannot pass this floor.
        q2 = F(316,100)
        assert count*(floor+F(1,10**7))-F(3,10)/q2-F(2,100) > F(271,100)
        crit[str(count)] = float(floor)
    c = F(32,100)
    compatible = [n for n in range(3,15) if witness(n,c) is not None]
    assert compatible == [3,4,5,6,7,8]
    return {"critical_floors":crit,"common_positive_systematic_witness_for_four":witness(4,upper/4),
            "at_response_point_32":compatible,
            "with_integral_charge_and_odd_count":[n for n in compatible if n % 6 == 3]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-new", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = run()
    path = Path(__file__).with_name("independent_results.json")
    if args.write_new:
        with path.open("x", encoding="utf8") as stream:
            json.dump(result,stream,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        assert result == json.loads(path.read_bytes())
    print(json.dumps(result))
