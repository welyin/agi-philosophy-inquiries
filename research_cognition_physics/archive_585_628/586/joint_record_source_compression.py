"""586: retain original matter instruments when compressing classical records.

Full-model statements are analytic. Numerical local differential diagnostics
are not a complete Gauss lattice propagation or a quantum gravity simulation.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_quantum_measure_records as prior_packet
import joint_local_lapse_source as lapse

HERE = Path(__file__).resolve().parent
TARGET = HERE / "joint_record_source_compression_results.json"
HBAR, WEIGHT = .7, .8
M = original.M


def instruments(s):
    """Two original 577 reads, grouped by parity; direct parity alternative."""
    t = np.sin(s)/2
    dt = np.cos(s)/2
    e = np.stack(((1+t)/2, (1-t)/2), axis=-1)
    de = np.stack((dt/2, -dt/2), axis=-1)
    f, df = [], []
    for a, b in ((0, 0), (0, 1), (1, 0), (1, 1)):
        f.append(e[..., a]*e[..., b])
        df.append(de[..., a]*e[..., b]+e[..., a]*de[..., b])
    f, df = np.stack(f, axis=-1), np.stack(df, axis=-1)
    c = np.stack((f[..., 0]+f[..., 3], f[..., 1]+f[..., 2]), axis=-1)
    dc = np.stack((df[..., 0]+df[..., 3], df[..., 1]+df[..., 2]), axis=-1)
    return np.sqrt(f), df/(2*np.sqrt(f)), np.sqrt(c), dc/(2*np.sqrt(c))


def delta_a(s):
    return np.cos(s)**2/(8*(1+np.sin(s)**2/4))


def metric_record_hessian(phi):
    f, s = original.F(phi), phi[-1]
    a = 1-s*s/(6*M)
    return -f*f*(a*a*np.sin(s)+a*s*np.cos(s)/(3*M))


def parity_check():
    s = np.linspace(-3.4, 3.4, 241)
    lf, dlf, lc, dlc = instruments(s)
    fine = np.sum(dlf*dlf, axis=-1)
    coarse = np.sum(dlc*dlc, axis=-1)
    residual = float(np.max(abs(fine-coarse-delta_a(s))))
    probability_error = float(np.max(abs(lc[:, 0]**2-(lf[:, 0]**2+lf[:, 3]**2))))
    assert residual < 2e-16 and probability_error < 3e-16
    # The odd branch really is identical; the even branch need not be.
    s1, s2 = .25, .55
    x, _, u, _ = instruments(np.array(s1))
    y, _, v, _ = instruments(np.array(s2))
    even_kernel_gap = float(u[0]*v[0]-x[0]*y[0]-x[3]*y[3])
    odd_kernel_gap = float(u[1]*v[1]-x[1]*y[1]-x[2]*y[2])
    assert even_kernel_gap > .004 and abs(odd_kernel_gap) < 1e-15
    return dict(pointwise_energy_coefficient_error=residual,
                parity_probability_error=probability_error,
                even_kernel_gap=even_kernel_gap, odd_kernel_gap=odd_kernel_gap,
                original_readout="E_r=1/2+r*sin(s)/4",
                zero_wait_composition=True)


def fisher_check():
    """Matrix chain rule, independently constructed multivariable effects."""
    rng = np.random.default_rng(586)
    errors, minimum = [], 1.
    for _ in range(40):
        x = rng.normal(size=3)/3
        coeff = rng.normal(size=(6, 3))
        z = coeff@x
        p = np.exp(z-z.max()); p /= p.sum()
        dp = p[:, None]*(coeff-p@coeff)
        fine = (dp.T/p)@dp
        cp = np.array([p[:3].sum(), p[3:].sum()])
        dcp = np.array([dp[:3].sum(axis=0), dp[3:].sum(axis=0)])
        coarse = (dcp.T/cp)@dcp
        conditional = np.zeros((3, 3))
        for a, group in enumerate((slice(0, 3), slice(3, 6))):
            q = p[group]/cp[a]
            dq = dp[group]/cp[a]-p[group, None]*dcp[a]/cp[a]**2
            conditional += cp[a]*(dq.T/q)@dq
        errors.append(float(np.max(abs(fine-coarse-conditional))))
        minimum = min(minimum, float(np.linalg.eigvalsh(fine-coarse).min()))
    assert max(errors) < 2e-15 and minimum > -1e-14
    # Positive control: within-group ratios independent of configuration.
    a = np.array([.2, .3, .5]); b = np.array([.6, .1, .3])
    def f(s):
        c = .5+.2*np.sin(s)
        return np.sqrt(np.r_[c*a, (1-c)*b])
    x, y = f(.2), f(.9)
    c1, c2 = .5+.2*np.sin(.2), .5+.2*np.sin(.9)
    same = max(abs(x[:3]@y[:3]-np.sqrt(c1*c2)),
               abs(x[3:]@y[3:]-np.sqrt((1-c1)*(1-c2))))
    assert same < 2e-16
    return dict(matrix_chain_rule_error=max(errors),
                minimum_loss_eigenvalue=minimum, constant_conditional_kernel_error=float(same))


def compact_packet(nodes, phase=.23):
    z, qw = np.polynomial.legendre.leggauss(nodes)
    x, y = np.meshgrid(z, z, indexing="ij")
    h = prior_packet.HCENTER+prior_packet.HRADIUS*x
    s = prior_packet.SCENTER+prior_packet.SRADIUS*y
    chi = np.exp(-1/(1-x*x)-1/(1-y*y))
    weights = qw[:, None]*qw[None, :]*prior_packet.HRADIUS*prior_packet.SRADIUS*h**3
    prob = weights*chi**2; prob /= prob.sum()
    coords = np.stack((h, s), axis=-1)
    f = M-(h*h+s*s)/6
    g = f[..., None, None]*(np.eye(2)-coords[..., :, None]*coords[..., None, :]/(6*M))
    # Actual curved wavefunction chi/sqrt(mu), with nonzero probability current.
    grad = np.stack((-2*x/(prior_packet.HRADIUS*(1-x*x)**2),
                     -2*y/(prior_packet.SRADIUS*(1-y*y)**2)), axis=-1)
    grad = grad-coords/(2*f[..., None])
    grad = grad.astype(complex)+1j*phase*np.stack((s, h), axis=-1)
    def kinetic(l, dl):
        vector = l[..., None]*grad
        vector[..., 1] += dl
        density = np.einsum("...i,...ij,...j->...", vector.conj(), g, vector).real
        return float(np.sum(prob*density))*HBAR**2/(2*WEIGHT)
    lf, dlf, lc, dlc = instruments(s)
    before = kinetic(np.ones_like(s), np.zeros_like(s))
    after_f = sum(kinetic(lf[..., j], dlf[..., j]) for j in range(4))
    after_c = sum(kinetic(lc[..., j], dlc[..., j]) for j in range(2))
    df = HBAR**2/(2*WEIGHT)*float(np.sum(prob*g[..., 1, 1]*np.sum(dlf**2, axis=-1)))
    dc = HBAR**2/(2*WEIGHT)*float(np.sum(prob*g[..., 1, 1]*np.sum(dlc**2, axis=-1)))
    delta = HBAR**2/(2*WEIGHT)*float(np.sum(prob*g[..., 1, 1]*delta_a(s)))
    a = 1-s*s/(6*M)
    hh = -f*f*(a*a*np.sin(s)+a*s*np.cos(s)/(3*M))
    acceleration_gap = HBAR**2/WEIGHT**2*float(np.sum(prob*hh*delta_a(s)))
    assert abs((after_f-before)-df) < 1e-10
    assert abs((after_c-before)-dc) < 1e-10
    assert abs((after_f-after_c)-delta) < 1e-10
    assert delta > 0 and acceleration_gap < 0
    return dict(nodes=nodes, phase=phase, kinetic_before=before,
                kinetic_after_fine=after_f, kinetic_after_direct_parity=after_c,
                energy_increase_fine=df, energy_increase_direct_parity=dc,
                energy_gap=delta, energy_direct_form_error=abs(after_f-after_c-delta),
                mean_b_acceleration_gap=acceleration_gap,
                future_original_probability_t2_coefficient=acceleration_gap/8,
                immediate_even_probability=float(np.sum(prob*lc[..., 0]**2)),
                immediate_original_plus_probability=float(np.sum(prob*(.5+np.sin(s)/4))))


def packet_check():
    rows = [compact_packet(n) for n in (48, 80, 128)]
    assert abs(rows[-1]["energy_gap"]-rows[-2]["energy_gap"]) < 1e-13
    assert abs(rows[-1]["mean_b_acceleration_gap"]-rows[-2]["mean_b_acceleration_gap"]) < 1e-12
    real = compact_packet(80, 0.)
    assert abs(real["energy_gap"]-rows[1]["energy_gap"]) < 1e-15
    return dict(rows=rows, real_source_energy_gap=real["energy_gap"],
                full_Gauss_packet_with_Haar_constant_links=True,
                other_full_H_terms_unchanged_by_pointwise_completeness=True,
                energy_not_a_discretized_spectrum=True)


def differential_check():
    """Independent full five-coordinate nested differential diagnostic.

    Other coordinates are parameters here; cancellations of their terms in the
    FULL graph commutator are proved in the note, not simulated by this test.
    """
    point = np.array([.14, .53, -.17, .1, .42])
    neighbor = np.array([.23, .41, -.09, .14, .31])
    unit = np.eye(5, dtype=int)
    origin = np.zeros(5, dtype=int)
    expected = HBAR**2/WEIGHT**2*delta_a(point[-1])*metric_record_hessian(point)
    rows = []
    for step in (.024, .012, .006):
        @lru_cache(None)
        def data(key):
            q = point+step*np.array(key)
            g = original.inverse(q)
            drift = -original.F(q)*q/(3*M)
            v = WEIGHT*float(original.node_potential(q))+float(original.distance_squared(q, neighbor))/2
            wave = np.exp(-.21*np.dot(q, q)+.13j*q[1]+.09j*q[2]*q[4])
            lf, _, lc, _ = instruments(q[-1])
            value = np.r_[lf, lc]*wave
            return g, drift, v, q[-1], np.r_[value, np.sin(q[-1])*value]

        def h_apply(fn, key):
            key = np.array(key)
            g, drift, v, _, _ = data(tuple(key))
            middle = fn(tuple(key)); result = v*middle
            lap = np.zeros_like(middle)
            for i in range(5):
                plus, minus = fn(tuple(key+unit[i])), fn(tuple(key-unit[i]))
                lap += g[i, i]*(plus-2*middle+minus)/step**2
                lap += drift[i]*(plus-minus)/(2*step)
                for j in range(i):
                    cross = (fn(tuple(key+unit[i]+unit[j]))-fn(tuple(key+unit[i]-unit[j]))
                             -fn(tuple(key-unit[i]+unit[j]))+fn(tuple(key-unit[i]-unit[j])))
                    lap += 2*g[i, j]*cross/(4*step**2)
            return result-HBAR**2*lap/(2*WEIGHT)

        @lru_cache(None)
        def first(key):
            return h_apply(lambda k: data(k)[-1], key)

        def second_input(key):
            first_value = first(key)
            return np.r_[first_value, np.sin(data(key)[3])*first_value[:6]]
        outer = h_apply(second_input, origin)
        d2 = -(outer[6:12]-2*outer[12:]+np.sin(point[-1])*outer[:6])/HBAR**2
        lf, _, lc, _ = instruments(point[-1])
        wave0 = np.exp(-.21*np.dot(point, point)+.13j*point[1]+.09j*point[2]*point[4])
        value = (lf@d2[:4]-lc@d2[4:])/wave0
        rows.append(dict(step=step, direct_nested_commutator_real=float(value.real),
                         imaginary_part=float(value.imag), error=float(abs(value-expected))))
    assert rows[-1]["error"] < rows[0]["error"]/9, rows
    assert rows[-1]["error"] < 3e-5, rows
    return dict(expected_acceleration_gap=float(expected), rows=rows,
                full_five_coordinate_local_operator=True, full_graph_propagation=False)


def geometry_check():
    packet = compact_packet(80)
    eps, psi0 = .8, 1.13
    numerator = packet["energy_gap"]*WEIGHT
    def gap(psi): return numerator/(eps**3*psi**6)
    expected = -6*gap(psi0)/psi0
    derivatives = []
    for step in (.004, .002, .001):
        value = (gap(psi0+step)-gap(psi0-step))/(2*step)
        derivatives.append(dict(step=step, derivative=value, error=abs(value-expected)))
    assert derivatives[-1]["error"] < derivatives[0]["error"]/12
    source = np.zeros((5, 5, 5)); source[0, 0, 0] = gap(psi0)
    smooth = lapse.average(source)
    x = np.arange(5)*2*np.pi/5
    N = 1+.2*np.cos(x[:, None, None])+np.zeros_like(source)
    value = float(np.sum(N*smooth))
    predicted = float(lapse.average(N)[0, 0, 0]*source[0, 0, 0])
    assert abs(value-predicted) < 1e-17
    return dict(node_energy_gap=gap(psi0), node_volume=eps**3*psi0**6,
                spatial_weight_partial_derivative=expected, derivatives=derivatives,
                smeared_lapse_energy_gap=value, lapse_duality_error=abs(value-predicted),
                local_source_sum_error=float(abs(smooth.sum()-source.sum())),
                fixed_state_partial_response_not_a_dynamical_Einstein_solution=True)


def run():
    evidence = dict(original_parity_instruments=parity_check(),
                    multivariable_fisher_chain_rule=fisher_check(),
                    original_compact_Gauss_packet=packet_check(),
                    original_local_nested_commutator=differential_check(),
                    same_source_geometry_response=geometry_check())
    deps = ("joint_curved_quantum_source.py", "joint_quantum_measure_records.py",
            "joint_local_lapse_source.py", "research_note_506.md", "research_note_557.md",
            "research_note_574.md", "research_note_577.md", "research_note_584.md",
            "research_note_585.md", "research_round_585_checks.json",
            "cognition_forward_bridge_review_585.md")
    return dict(round=586, tests_run=len(evidence), failures=0, errors=0, evidence=evidence,
                dependency_hashes={n: hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope="Original fixed finite curved-target Gauss matter. Coarse-record instrument criterion, source energy and same-H subsequent-record witness; ideal instruments and preparation remain inputs. No new autonomous detector, quantum Einstein solution or continuum limit.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results", action="store_true")
    args = parser.parse_args(); result = run()
    if args.write_results:
        with TARGET.open("x", encoding="utf8", newline="\n") as stream:
            stream.write(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
    else:
        assert json.loads(TARGET.read_text(encoding="utf8")) == result
    print(json.dumps({k: result[k] for k in ("round", "tests_run", "failures", "errors")}))
