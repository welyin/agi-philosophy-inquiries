"""Separate finite-outcome check of proof (11); no operational norm claim."""
from pathlib import Path
import argparse
import json
import numpy as np
import composition_stability as core


def instrument(rng):
    v = rng.normal(size=(6, 2)) + 1j*rng.normal(size=(6, 2))
    v, _ = np.linalg.qr(v)
    branches = []
    for i in range(3):
        k = v[2*i:2*i+2]
        z = k.T.reshape(-1)
        branches.append(np.outer(z, z.conj()))
    total = sum(branches)
    assert np.max(abs(np.einsum('iaja->ij', total.reshape(2,2,2,2))-np.eye(2))) < 1e-12
    return branches


def calculate():
    rng = np.random.default_rng(105711)
    words = core.R_TYPES + core.P_TYPES + core.Q_TYPES
    basis = np.array([core.pauli_word(w) for w in words[1:]])
    max_normalization_error = max_embedding_error = max_tv_ratio = 0.0
    min_complement_eigenvalue = 1.0
    for _ in range(24):
        h = np.einsum('a,aij->ij', rng.normal(size=87), basis)
        h *= 0.9/np.max(abs(np.linalg.eigvalsh(h)))
        w = (np.eye(16)+h)/16
        ordered = min((core.twirl(w,(1,)), core.twirl(w,(3,))),
                      key=lambda v: core.trace_distance(w,v))
        distance = core.trace_distance(w, ordered)
        for _ in range(4):
            aa, bb = instrument(rng), instrument(rng)
            effects = [np.kron(a,b)/4 for a in aa for b in bb]
            complement = np.eye(16)-sum(effects)
            min_complement_eigenvalue = min(min_complement_eigenvalue,
                                             float(np.linalg.eigvalsh(complement).min()))
            p = np.array([core.pairing(4*w, 4*e) for e in effects])
            q = np.array([core.pairing(4*ordered, 4*e) for e in effects])
            assert min(p.min(),q.min()) >= -1e-12
            max_normalization_error = max(max_normalization_error,abs(p.sum()-1),abs(q.sum()-1))
            pp = np.array([core.pairing(w,e) for e in effects]+[core.pairing(w,complement)])
            qq = np.array([core.pairing(ordered,e) for e in effects]+[core.pairing(ordered,complement)])
            max_embedding_error = max(max_embedding_error, float(np.max(abs(pp[:-1]-p/16))),
                                       abs(pp[-1]-15/16), abs(qq[-1]-15/16))
            tv = float(np.sum(abs(p-q))/2)
            embedded_tv = float(np.sum(abs(pp-qq))/2)
            assert abs(tv/16-embedded_tv) < 1e-12
            assert embedded_tv <= distance+1e-12
            assert tv <= 16*distance+1e-12
            max_tv_ratio = max(max_tv_ratio, tv/(16*distance))
    assert min_complement_eigenvalue > -1e-12
    assert max_normalization_error < 1e-12 and max_embedding_error < 1e-12
    return {'all_checks_passed':True,'processes':24,'instrument_pairs':96,
            'outcomes_per_pair':9,'povm_embedding_scale':16,
            'max_probability_normalization_error':max_normalization_error,
            'max_embedding_error':max_embedding_error,
            'min_complement_eigenvalue':min_complement_eigenvalue,
            'max_ratio_TV_over_16T':max_tv_ratio,
            'scope':'single-copy finite local CP instrument classical outcomes only',
            'proof_of_universal_bound_is_analytic':True}


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    result=calculate()
    path=Path(__file__).with_name('record_bound_results.json')
    if args.save:
        with path.open('x',encoding='utf-8') as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write('\n')
    else:
        core.compare(json.loads(path.read_text(encoding='utf-8')),result)
    print(json.dumps(result,ensure_ascii=False))
