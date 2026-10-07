"""Working 790: first Ward insertion of a quadratic coordinate vector field.

Reuses the finite exact BV algebra, not a continuous functional Laplacian.
Tests the commutator of Wick contraction with the full free BV differential.
"""
from fractions import Fraction as Q
from pathlib import Path
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('bv777_generator', HERE.parent/'777/joint_auxiliary_bv_reduction.py')
bv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bv)
v = {name: bv.var(name) for name in bv.NAMES}


def run():
    # b is a flat physical field w in THIS diagnostic; rho=b* is retained.
    # No nonminimal b,h,a sector is included in this finite free model.
    q, r, w = [v[name] for name in ('q', 'r', 'b')]
    s2 = bv.add(bv.scale(bv.mul(r, r), Q(1, 2)), bv.mul(v['p'], v['c']))
    s0 = lambda f: bv.bracket(s2, f)
    vector = [bv.add(bv.mul(q, r), bv.scale(bv.mul(w, w), Q(1, 2))),
              bv.add(bv.mul(q, r), bv.mul(r, w), bv.scale(bv.mul(q, q), Q(1, 2))),
              bv.add(bv.mul(w, r), bv.scale(bv.mul(r, r), Q(1, 2)))]
    generator = bv.add(*(bv.mul(v[anti], f) for anti, f in zip(('p', 't', 'rho'), vector)))
    names = ('q', 'r', 'b')
    h = [[Q(3, 5), Q(2, 9), Q(1, 6)],
         [Q(2, 9), Q(2, 7), Q(4, 11)],
         [Q(1, 6), Q(4, 11), Q(1, 3)]]
    project = lambda f: bv.substitute(f, {anti: {} for _, anti in bv.PAIRS})
    rows, anomalies, means = [], [], []
    for sigma in (Q(5, 4), Q(7, 6)):
        covariance = [[1+sigma, Q(0), sigma/3], [Q(0), Q(0), Q(0)],
                      [sigma/3, Q(0), sigma]]
        c = [[covariance[i][j]-h[i][j] for j in range(3)] for i in range(3)]

        def gamma(f, kernel=c):
            return bv.add(*(bv.scale(bv.diff(bv.diff(f, names[i]), names[j]), kernel[i][j]/2)
                            for i in range(3) for j in range(3)))

        # The full differential, including s(r*)=r, is essential here.
        anomaly = bv.add(s0(gamma(generator)), bv.scale(gamma(s0(generator)), -1))
        k = [gamma(f) for f in vector]
        dk = bv.mul(r, k[1])
        rho = bv.add(project(gamma(s0(generator))), bv.scale(dk, -1))
        assert project(anomaly) == bv.scale(rho, -1)
        assert anomaly == project(anomaly) and rho
        assert not bv.add(s0(gamma(generator, covariance)),
                          bv.scale(gamma(s0(generator), covariance), -1))
        interaction = bv.scale(bv.prod(w, w, w), Q(5, 42))
        ell = bv.scale(w, Q(7, 13))
        old_raw = project(gamma(interaction))
        new_raw = project(gamma(bv.add(interaction, bv.scale(s0(generator), -1))))
        new_ell = bv.add(ell, rho)
        assert bv.add(new_raw, new_ell, dk) == bv.add(old_raw, ell)
        assert bv.add(new_raw, ell, dk) != bv.add(old_raw, ell)
        anomalies.append(anomaly)
        means.append(k)
        rows.append(dict(sigma=str(sigma), contracted_generator=bv.display(gamma(generator)),
                         local_first_ward_anomaly=bv.display(anomaly),
                         source_contact=bv.display(rho), source_identity_residual='0'))
    assert anomalies[0] == anomalies[1] and means[0] != means[1]
    return dict(round=790, status='working_probe_not_signed_off', cases=rows,
                equation='A1(G3)|antifields=0 = -rho',
                source_transport_sign_verified=True, state_independent_anomaly=True,
                state_dependent_mean=True,
                original_full_continuum_path_source_matching_proven=False,
                auxiliary_gauge_endpoint_matching_proven=False)


if __name__ == '__main__':
    result = run()
    HERE.joinpath('generator_source_probe_results.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
