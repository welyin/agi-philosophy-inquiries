"""830: original finite-graph energy-moment obstruction to an exact encoder.

Sparse CAR tests retain the original 32 species/spin modes at every node.
They check coefficient identities, not bosonic time evolution. The full
unbounded Hamiltonian and its domain are treated analytically in the note.
"""
from pathlib import Path
import argparse, json, sys
import numpy as np
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
TARGET = HERE / 'native_encoding_energy_moments_results.json'
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(HERE.parent / '829'))
from research_layout import Layout, ResearchRuntime
import majorana_code_source_bridge as code


def polynomial_checks():
    gamma, _, _, _, compress, _, _ = code.code_data()
    # Three commuting Hermitian bilinears; their product is a logical word
    # already identified in 829. No boson coefficients are commuted here.
    terms = [code.mul((0, 0, 1), code.mul(gamma[i], gamma[j]))
             for i, j in ((0, 1), (4, 5), (12, 14))]
    assert all(code.adj(a) == a for a in terms)
    assert all(code.commute(a, b) for a in terms for b in terms)
    assert all(compress(a)[0] == 'zero' for a in terms)
    square = [code.mul(a, b) for a in terms for b in terms]
    assert sum(compress(a) == ('scalar', 0) for a in square) == 3
    assert sum(compress(a)[0] == 'zero' for a in square) == 6
    cube = [code.product([a, b, c]) for a in terms for b in terms for c in terms]
    logical = [a for a in cube if compress(a)[0] == 'logical']
    assert len(logical) == 6 and len(set(logical)) == 1
    assert all(compress(a)[0] in ('logical', 'zero') for a in cube)
    return dict(compressed_H=0, compressed_H_squared=3,
                compressed_H_cubed_logical_coefficient=6,
                square_of_compressed_H_is_not_compressed_H_squared=True,
                higher_energy_moments_not_claimed_logically_blind=True,
                example_is_algebraic_not_the_original_H=True)


def packet(order, y_nu):
    # Same normal strict-Gauss scalar packet as 717/743/746: R=0,
    # |x_H| < .7, |x5| < 1, H^5 volume d^5x/sqrt(1+|x|^2/6).
    u, wu = np.polynomial.legendre.leggauss(order)
    radial = .7 * (u + 1) / 2
    wr = .7 * wu / 2
    r, s = np.meshgrid(radial, u, indexing='ij')
    weight = wr[:, None] * wu[None, :] * r**3 / np.sqrt(1 + (r*r+s*s)/6)
    weight *= np.exp(-2/(1-(r/.7)**2)-2/(1-s*s))
    weight /= weight.sum()
    return dict(mean_x5=float(np.sum(weight*s)),
                mean_Higgs_radius_squared=float(np.sum(weight*r*r)),
                original_Dirac_lower_bound_on_H2_difference=float(
                    2*abs(y_nu)**2*np.sum(weight*r*r)))


def native_checks(old):
    rng = np.random.default_rng(830)
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.diag([1., -1.])
    pair = (1 << 30) | (1 << 31)
    rows = []
    for nodes in (1, 2, 5):
        n = 32 * nodes
        edges = [] if nodes == 1 else ([(0, 1)] if nodes == 2 else
                                     [(v, (v+1) % nodes) for v in range(nodes)])
        for trial in range(4):
            x = rng.normal(size=(nodes, 5)) * .6
            h = np.zeros((n, n), complex)
            d = np.zeros_like(h)
            for v in range(nodes):
                sl = slice(32*v, 32*(v+1))
                h[sl, sl], d[sl, sl] = old.mass_x(x[v])
            hopping = np.zeros_like(h)
            for v, w in edges:
                C = old.matter.gauge.group_exp(rng.normal(size=8)*.3, 3)
                W = old.matter.gauge.group_exp(rng.normal(size=3)*.3, 2)
                rep = old.matter.representation(C, W, np.exp(1j*rng.normal()))
                spin = sum(z*p for z, p in zip(rng.normal(size=3), (sx, sy, sz)))
                edge = (.3+.17j)*rep@np.kron(np.eye(16), spin)
                sv, sw = slice(32*v, 32*(v+1)), slice(32*w, 32*(w+1))
                hopping[sv, sw] = edge
                hopping[sw, sv] = edge.conj().T
            zero = np.zeros_like(h)
            apply = lambda state, hh, dd: old.car.quadratic(state, hh, dd)
            dot = old.dot
            norm2 = lambda state: float(dot(state, state).real)
            vacuum, occupied = {0: 1.}, {pair: 1.}
            m0, m1 = apply(vacuum, zero, d), apply(occupied, zero, d)
            d1, c1 = apply(occupied, h, zero), apply(occupied, hopping, zero)
            b0, b1 = apply(vacuum, h+hopping, d), apply(occupied, h+hopping, d)
            dirac_expected = 2*abs(old.matter.Y['nu'])**2*float(x[0, :4]@x[0, :4])
            hopping_expected = float(np.linalg.norm(hopping[:, 30:32])**2)
            delta = norm2(b1)-norm2(b0)
            errors = [abs(norm2(m1)-norm2(m0)),
                      abs(norm2(d1)-dirac_expected),
                      abs(norm2(c1)-hopping_expected),
                      abs(dot(d1, c1)), abs(dot(m1, d1)), abs(dot(m1, c1)),
                      abs(dot(occupied, b1)), abs(dot(vacuum, b0)),
                      abs(delta-dirac_expected-hopping_expected)]
            # H_b acts on the scalar wavefunction, with the same fermion
            # basis state. These orthogonalities remove all diagonal cross
            # terms with H_b; H_b is not set to zero in the analytic identity.
            assert max(errors) < 2e-12
            assert norm2(apply(vacuum, h+hopping, zero)) == 0
            assert dirac_expected > 0 and delta > 0
            rows.append(dict(nodes=nodes, CAR_modes=n, links=len(edges), trial=trial,
                             original_Dirac_contribution=dirac_expected,
                             original_hopping_contribution=hopping_expected,
                             H2_diagonal_difference_at_configuration=delta,
                             max_identity_residual=float(max(errors))))
    fine, coarse = packet(80, old.matter.Y['nu']), packet(48, old.matter.Y['nu'])
    change = max(abs(fine[k]-coarse[k]) for k in fine)
    assert change < 1e-8 and abs(fine['mean_x5']) < 1e-15
    assert fine['original_Dirac_lower_bound_on_H2_difference'] > 0
    return dict(full_original_coefficient_rows=rows,
                original_Ynu_abs_squared=float(abs(old.matter.Y['nu'])**2),
                original_Ys_abs_squared=float(abs(old.matter.Y['s'])**2),
                old_normal_Gauss_packet=fine, quadrature_change_48_to_80=change,
                boson_H_and_all_interactions_retained_analytically=True,
                no_boson_time_evolution_or_continuum_simulation=True,
                one_and_two_node_checks_are_coefficient_checks_not_full_code_carriers=True,
                five_node_case_has_ten_original_sterile_code_modes=True)


def run():
    with ResearchRuntime(Layout()).installed():
        import joint_record_mass_feedback as old
        native = native_checks(old)
    return dict(round=830, all_checks_passed=True, fresh_test_groups=1,
                energy_moment_algebra=polynomial_checks(), native=native,
                scope='For the old normal-Gauss vacuum/pair input and a fixed original finite graph, exact closed energy-conserving encoding into the 829 degree-four-blind code is impossible: input H^2 moments differ, output moments cannot. Full native boson dynamics is retained in the analytic proof. This does not exclude approximate or relational encoding, other input preparations, or the original nonstationary continuum model.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write:
        assert not TARGET.exists(), 'Do not overwrite saved science.'
        TARGET.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
    else:
        assert result == json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(dict(round=830, all_checks_passed=True,
                         polynomial=result['energy_moment_algebra'],
                         normal_packet=result['native']['old_normal_Gauss_packet'],
                         maximum_sparse_residual=max(r['max_identity_residual'] for r in result['native']['full_original_coefficient_rows'])),
                     ensure_ascii=False, indent=2))
