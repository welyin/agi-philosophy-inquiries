"""713 entry: inherited680 congruence applied to the actual710 QQQL tensor."""
import argparse
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_charge_changing_vertex as original
TARGET=HERE/'local_interaction_entry_results.json'


def mul(a,b):
    out={}
    for x,c in a.items():
        for y,d in b.items():
            if x&y:continue
            crossings=sum((y&((1<<i)-1)).bit_count() for i in range(32) if x>>i&1)
            out[x|y]=out.get(x|y,0)+(-1)**crossings*c*d
    return {k:v for k,v in out.items() if v}


def run():
    p={sum(1<<i for i in ids):Fraction(c) for ids,c in original.COEFF.items()}
    powers=[{0:Fraction(1)}]
    for _ in range(5):powers.append(mul(powers[-1],p))
    assert powers[4] and not powers[5]
    plus={};minus={}
    for k,term in enumerate(powers[:5]):
        for mask,c in term.items():
            plus[mask]=plus.get(mask,0)+c/factorial(k)
            minus[mask]=minus.get(mask,0)+(-1)**k*c/factorial(k)
    assert mul(plus,minus)==mul(minus,plus)=={0:1}
    for k,term in enumerate(powers[:5]):assert all(mask.bit_count()==4*k for mask in term)
    names=('research_note_680.md','research_note_689.md','research_note_690.md',
           'research_note_691.md','research_note_699.md','research_note_710.md','research_note_712.md')
    return dict(entry_round=713,new_formal_round=False,original_raw_tensor_terms=len(p),
        nonzero_terms_by_power=[len(t) for t in powers],
        top_power_coefficients={str(k):str(v) for k,v in powers[4].items()},
        exponential_terms=len(plus),exact_inverse_product_is_one=True,
        nilpotency_order=5,coefficient_normalization_absorbable_into_coupling=True,
        scope='Actual holomorphic QQQL wedge example; general even positive-degree local interaction follows inherited680 polynomial congruence. Fixed candidate/Hb only; no new all-theory no-go.',
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert json.loads(TARGET.read_text('utf8'))==r
    print(json.dumps(r,ensure_ascii=False,indent=2))
