"""Reconstruct two-axis conventions, independent holdout and scaling degeneracy."""
import argparse
import itertools
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run():
    gf, mz, v, a = 1.166379e-5, 91.1876, -.03816, -.50111
    alpha_ext = 128.946
    unit = math.sqrt(math.sqrt(2)*gf)*mz
    left, right = unit*(v+a), unit*(v-a)
    nu = right-left
    triple = right-2*left
    g = math.sqrt(2*triple*nu)
    e = math.sqrt(triple*right)
    sine, cosine = e/g, triple/g
    assert abs(cosine*cosine+sine*sine-1) < 1.e-14
    gp = e/cosine
    assert abs(g*(-.5)*cosine + gp*(-.5)*(-sine)-left) < 2e-15
    assert abs(g*.5*cosine + gp*(-.5)*(-sine)-nu) < 2e-15
    assert abs(gp*sine-right) < 2e-15
    ranges = []
    for k in (1,2):
        inv = []
        for radius, vec, axial in itertools.product(
                [mz-k*.0021,mz+k*.0021], [v-k*.00047,v+k*.00047], [a-k*.00035,a+k*.00035]):
            inv.append(4*math.pi/(radius**2*math.sqrt(2)*gf*(3*axial**2-2*axial*vec-vec**2)))
        ranges.append((min(inv),max(inv),(alpha_ext-k*.015)/max(inv)))
    chiral = [nu*x/(mz*mz*math.sqrt(2)*gf) for x in (left,right)]
    predicted_v = sum(chiral)/2
    predicted_a = (chiral[0]-chiral[1])/2
    reported = [-.035,-.503]
    assert all(abs(p-x)<.017 for p,x in zip([predicted_v,predicted_a],reported))
    factor = 1.1
    original = [triple*left,triple*right]
    alternative = [(triple/factor)*(factor*x) for x in (left,right)]
    assert max(abs(x-y) for x,y in zip(original,alternative)) < 1.e-15
    return dict(canonical_angles=[sine,cosine],triple=triple,nu_mu=nu,alpha_boxes=ranges,
                charm_ii_VA=[predicted_v,predicted_a],charm_ii_chiral_moduli=list(map(abs,chiral)),
                observed_VA=reported,reported_marginal_widths=[.017,.017],
                erasing_one_parameter_template=dict(scale=factor,triple=triple/factor,
                    products=alternative,pole_electron_width_ratio=factor**2))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=json.loads(json.dumps(run()))
    target=ROOT/"independent_results.json"
    if args.write:
        with target.open("x",encoding="utf-8") as f:
            json.dump(result,f,indent=2)
            f.write("\n")
    else:
        assert result==json.loads(target.read_text(encoding="utf-8-sig"))
        print("1475 independent check passed")
