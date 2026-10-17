"""Exact finite checks for a conditional causal source identity, not a QFT simulation."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

D = Path(__file__).resolve().parent
S = D.parents[1]


def calc():
    base = json.loads((D/'source_baseline.json').read_text(encoding='utf-8'))
    for rel, h in base['source_sha256'].items():
        assert hashlib.sha256((S/rel).read_bytes()).hexdigest() == h, rel
    for name, h in base['external_source_sha256'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == h, name
    a, xc = F(1), F(1, 2)
    gaps = []
    for tp in [F(0), F(1, 3), F(1), F(3), F(10)]:
        for fraction in [F(0), F(1, 3), F(1)]:
            to = fraction*tp
            for y in [a, F(5, 4), F(2)]:
                for xo in [F(0), F(1, 4), xc]:
                    gap = (tp+y)-xo-(tp-to)
                    assert gap == y-xo+to >= a-xc
                    gaps.append(gap)
    # Exact quadratic minimal stress, including mass: indices in an orthonormal frame.
    metric = [-1, 1, 1, 1]
    m2 = F(2)

    def tensor(q, d):
        contraction = sum(metric[i]*d[i]*d[i] for i in range(4))+m2*q*q
        return [[d[i]*d[j]-(F(metric[i], 2)*contraction if i == j else 0)
                 for j in range(4)] for i in range(4)]

    def cross(q, d, f, df):
        contraction = sum(metric[i]*d[i]*df[i] for i in range(4))+m2*q*f
        return [[d[i]*df[j]+df[i]*d[j]-(metric[i]*contraction if i == j else 0)
                 for j in range(4)] for i in range(4)]

    cases = 0
    for q, d, f, df, eta in [
        (F(2,3), [F(1),F(2),F(-1),F(3,2)], F(3,5), [F(2),F(-1),F(0),F(1)], F(2,7)),
        (F(0), [F(0)]*4, F(3,5), [F(2),F(-1),F(0),F(1)], F(2,7)),
        (F(2,3), [F(1),F(2),F(-1),F(3,2)], F(0), [F(0)]*4, F(2,7)),
    ]:
        full = tensor(q+eta*f, [d[i]+eta*df[i] for i in range(4)])
        ref, pulse, mixed = tensor(q,d), tensor(f,df), cross(q,d,f,df)
        for i in range(4):
            for j in range(4):
                assert full[i][j]-ref[i][j] == eta*mixed[i][j]+eta**2*pulse[i][j]
                if q == 0 and not any(d) or f == 0 and not any(df):
                    assert mixed[i][j] == 0
        cases += 1
    # Commuting oscillator translations: exact coherent-state mean energy increment.
    wa, ub = F(1), F(1)
    actual = ((wa+ub)**2-ub**2)/2
    bare = wa**2/2
    assert actual == F(3,2) and bare == F(1,2) and actual-bare == wa*ub
    # A finite diagnostic of nonselective versus conditional mean; not a field apparatus.
    bell = np.array([1.,0.,0.,1.])/np.sqrt(2)
    X = np.array([[0.,1.],[1.,0.]])
    branches = []
    for sign in [1,-1]:
        effect = (np.eye(2)+sign*X)/2
        state = np.kron(np.eye(2),effect)@bell
        probability = float(state@state)
        mean = float(state@(np.kron(X,np.eye(2))@state)/probability)
        assert abs(probability-.5) < 1e-14 and abs(mean-sign) < 1e-14
        branches.append([probability,mean])
    assert abs(sum(p*m for p,m in branches)) < 1e-14
    return {
        'status':'PASS_WITH_SCOPE',
        'exact':{'causal_margin':str(a-xc), 'minimum_checked_gap':str(min(gaps)),
                 'support_diagnostics':len(gaps), 'stress_identity_cases':cases,
                 'commuting_shift_actual_increment':str(actual),
                 'commuting_shift_bare_increment':str(bare),
                 'commuting_shift_cross_increment':str(actual-bare)},
        'diagnostics':{'conditional_X_probability_and_mean':branches},
        'interpretation':{'general_causal_proof':'README equations (3)-(6)',
                          'new_QFT_device':False,'finite_G_error_certificate':False,
                          'same_noise_or_branch_sources_claimed':False,
                          'new_formal_round':False},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--write',action='store_true')
    args = ap.parse_args()
    result = calc()
    out = D/'results.json'
    if args.write:
        with out.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        mode = 'created_once'
    else:
        assert json.loads(out.read_text(encoding='utf-8')) == result
        mode = 'read_only_match'
    print(json.dumps({'mode':mode,**result},ensure_ascii=False))


if __name__ == '__main__':
    main()
