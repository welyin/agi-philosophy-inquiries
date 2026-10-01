"""557: one dynamic SU(2) link, physical edge records and source backreaction.

The graph, group, electric kinetic term and ideal neutral pointer coupling are
model inputs. The pulse is not autonomous evolution under the source H.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import joint_matter_energy_moment_control as poly

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_gauge_link_reference_results.json'


def left(q):
    w, x, y, z = q
    return np.array([[w, -x, -y, -z], [x, w, -z, y],
                     [y, z, w, -x], [z, -y, x, w]], dtype=float)


J = [left(np.eye(4)[i]) / 2 for i in (1, 2, 3)]


def geometry_checks():
    rng = np.random.default_rng(557)
    worst = dict(invariance=0., gradient_identity=0., finite_difference=0.)
    for _ in range(60):
        x, y = rng.normal(size=(2, 4))
        quats = rng.normal(size=(3, 4))
        quats /= np.linalg.norm(quats, axis=1)[:, None]
        u, g, h = map(left, quats)
        c = x @ u @ y
        f = np.sum((x-u@y)**2) / 2
        f2 = np.sum((g@x-g@u@h.T@h@y)**2) / 2
        worst['invariance'] = max(worst['invariance'], abs(f2-f))
        group = np.array([-x @ j @ u @ y for j in J])
        identity = (np.dot(x,x)*np.dot(y,y)-c*c)/4
        worst['gradient_identity'] = max(worst['gradient_identity'], abs(group@group-identity))
        eps = 1e-5
        for i, j in enumerate(J):
            plus = math.cos(eps/2)*np.eye(4)+2*math.sin(eps/2)*j
            minus = math.cos(eps/2)*np.eye(4)-2*math.sin(eps/2)*j
            fd = (np.sum((x-plus@u@y)**2)-np.sum((x-minus@u@y)**2))/(4*eps)
            worst['finite_difference'] = max(worst['finite_difference'], abs(fd-group[i]))
    assert worst['invariance'] < 1e-12
    assert worst['gradient_identity'] < 1e-12
    assert worst['finite_difference'] < 2e-9
    casimir = -sum(j@j for j in J)
    assert np.max(abs(casimir-.75*np.eye(4))) < 1e-14
    return worst


def heat_checks():
    path = HERE/'round557_drafts/link_heat_kernel_probe.py'
    spec = importlib.util.spec_from_file_location('heat_probe_557', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    result = module.run()
    assert result == json.loads((path.parent/'link_heat_kernel_probe_results.json').read_text('utf8'))
    # Conditional first moments at a fixed source configuration, independent
    # centered matter noise. This is a commuting position POVM, not a Wigner law.
    x = np.array([.4, -.7, 1.1, .2])
    y = np.array([-.2, .8, .1, -.5])
    tq, v, kappa = .3, .7, 1.2
    expected = kappa*v/2*np.sum((x-y)**2)
    for row in result['rows']:
        eta = row['fundamental_attenuation']
        measured = kappa*v/2*(x@x+y@y+8*tq-8*tq-2*eta*(x@y)/math.exp(-3*row['t']/4))
        assert abs(measured-expected) < 1e-12
    return result


def pointer_kernel_checks():
    nodes, weights = np.polynomial.hermite.hermgauss(50)
    weights = weights/math.sqrt(math.pi)
    maxerr = 0.
    rows = []
    for f, nu in ((0., .3), (.8, .6), (3.2, 1.1)):
        # y=F+noise integrates K_y(F)^2 dy. F derivative of log K_y is
        # (y-F)/(2 nu^2), so both the cross and energy terms are independent.
        z = math.sqrt(2)*nu*nodes
        y = f+z
        gradient = z/(2*nu*nu)
        computed = np.array([weights@y, weights@(y*y), weights@gradient,
                             weights@(gradient*gradient)])
        predicted = np.array([f, f*f+nu*nu, 0., 1/(4*nu*nu)])
        maxerr = max(maxerr, float(np.max(abs(computed-predicted))))
        # Off-diagonal kernel verifies a nontrivial state update, not identity.
        gap = .9
        actual = float(weights @ np.exp(gap*z/(2*nu*nu)-gap*gap/(4*nu*nu)))
        target = math.exp(-gap*gap/(8*nu*nu))
        assert abs(actual-target) < 1e-14
        rows.append(dict(F=f, noise_variance=nu*nu, moments=computed.tolist(),
                         off_diagonal_factor=target))
    assert maxerr < 1e-12
    return dict(maximum_moment_residual=maxerr, cases=rows)


def source(k, endpoint, v, kappa, b, nu, omega=1.3):
    """Radial polynomial matter source and constant link; all ten coordinates.

    Haar invariance sets U=I inside integrals because this initial source is
    separately radial at the two ends. Group derivatives remain in the energy.
    The eleventh Gaussian integration variable represents pointer momentum.
    """
    d = 11
    a = 1/(2*v)
    variance = 1/(2*omega)
    q = [poly.var(i,d) for i in range(d)]
    radii = [poly.add(*(poly.mul(q[j],q[j]) for j in range(i,i+4))) for i in (0,5)]
    P = poly.const(1,d)
    if k == 1:
        P = poly.scale(poly.add(poly.const(2,d),poly.scale(radii[endpoint],-omega)),1/math.sqrt(2))
    elif k != 0:
        raise ValueError(k)
    assert abs(poly.norm2(P,variance)-1) < 1e-12
    density = poly.mul(poly.conj(P),P)
    L,C,u = poly.matter()
    ell = float(np.linalg.eigvalsh(L)[0])
    W = {}
    for i, offset in enumerate((0,5)):
        delta = [poly.add(radii[i],poly.const(-u[0],d)),
                 poly.add(poly.mul(q[offset+4],q[offset+4]),poly.const(-u[1],d))]
        W = poly.add(W,poly.scale(poly.quadratic(delta,L),v/4))
    diff = [poly.add(q[i],poly.scale(q[5+i],-1)) for i in range(4)]
    F = poly.scale(poly.add(*(poly.mul(z,z) for z in diff)),kappa*v/2)
    c = poly.add(*(poly.mul(q[i],q[5+i]) for i in range(4)))
    df = [poly.deriv(F,i) for i in range(10)]
    dg = [poly.scale(poly.add(*(poly.scale(poly.mul(q[i],q[5+j]),gen[i,j])
                               for i in range(4) for j in range(4) if gen[i,j])), -kappa*v)
          for gen in J]
    grad_m = poly.add(*(poly.mul(z,z) for z in df))
    grad_u = poly.add(*(poly.mul(z,z) for z in dg))
    rhs_m = poly.scale(F,4*kappa*v)
    rhs_u = poly.scale(poly.add(poly.mul(radii[0],radii[1]),
                               poly.scale(poly.mul(c,c),-1)),(kappa*v)**2/4)
    residual = max([abs(z) for z in poly.add(grad_m,poly.scale(rhs_m,-1),
                            grad_u,poly.scale(rhs_u,-1)).values()]+[0.])
    assert residual < 1e-13
    mean = lambda expr: poly.real(poly.expect(poly.mul(density,expr),variance))
    G = poly.add(poly.scale(grad_m,a),poly.scale(grad_u,b))
    derivatives = [poly.add(poly.deriv(P,i),poly.scale(poly.mul(q[i],P),-omega)) for i in range(10)]
    initial_T = a*sum(poly.norm2(z,variance) for z in derivatives)
    E1 = initial_T+mean(poly.add(W,F))
    first, second, response = mean(F),mean(poly.mul(F,F)),mean(G)
    bound_F2 = 8*kappa*kappa*v*v*u[0]**2+16*kappa*kappa*v*E1/ell
    bound_G = 4*a*kappa*v*E1+b*kappa*kappa*v*v*u[0]**2/2+b*kappa*kappa*v*E1/ell
    assert second <= bound_F2 and response <= bound_G
    # Independent direct integration of derivatives of the post-pulse full
    # source-plus-pointer wavefunction in pointer momentum representation.
    pp = poly.scale(q[10],1/(2*nu*math.sqrt(variance)))
    shifted = [poly.add(z,poly.scale(poly.mul(poly.mul(pp,f),P),-1j))
               for z,f in zip(derivatives,df)]
    after_T = a*sum(poly.norm2(z,variance) for z in shifted)
    after_el = b*sum(poly.norm2(poly.mul(poly.mul(pp,z),P),variance) for z in dg)
    delta = after_T+after_el-initial_T
    expected_delta = response/(4*nu*nu)
    assert abs(delta-expected_delta) < 1e-9*max(1.,expected_delta)
    assert abs((after_T-initial_T)-mean(poly.scale(grad_m,a))/(4*nu*nu)) < 1e-9
    # A radial Gaussian witness admits independent closed source moments.
    analytic_residual = None
    if k == 0:
        analytic = [2*kappa*v/omega, 6*(kappa*v)**2/omega**2,
                    8*a*kappa*kappa*v*v/omega+3*b*(kappa*v)**2/(4*omega**2)]
        analytic_residual = float(np.max(abs(np.array([first,second,response])-analytic)))
        assert analytic_residual < 1e-10
    radius = math.sqrt((bound_F2+nu*nu)/.05)
    return dict(radial_excitation=k,excited_endpoint=endpoint,omega=omega,
                v=v,kappa=kappa,electric_coefficient=b,noise_variance=nu*nu,
                E1=E1,edge_mean=first,edge_second_moment=second,gradient_budget=response,
                source_energy_increase=delta,post_read_source_energy=E1+delta,
                matter_energy_increase=after_T-initial_T,electric_energy_increase=after_el,
                independent_energy_identity_residual=abs(delta-expected_delta),
                polynomial_gradient_residual=float(residual),
                gaussian_closed_form_residual=analytic_residual,
                bounds=dict(edge_second=float(bound_F2),gradient=float(bound_G),
                    source_energy_increase=float(bound_G/(4*nu*nu))),
                actual_record_second_moment=second+nu*nu,
                record_radius_at_failure_0_05=radius,
                witnessed_Markov_bound=(second+nu*nu)/radius**2)


def run():
    geometry = geometry_checks()
    heat = heat_checks()
    pointer = pointer_kernel_checks()
    rows = [source(0,0,1.,.8,.3,.6), source(1,0,.7,1.1,.4,.3),
            source(1,1,1.5,.5,.7,.9)]
    checks = ['local_SU2_covariance_Casimir_and_independent_group_gradients',
              'positive_link_kernel_normalization_attenuation_and_unbiased_edge_probability',
              'constant_link_physical_source_and_exact_square_root_instrument_leakage',
              'same_quartic_E1_bounds_for_edge_reading_and_source_backreaction',
              'neutral_pointer_probability_state_update_and_energy_integrals',
              'full_post_pulse_derivative_energy_including_matter_and_dynamic_link']
    dependencies = ('joint_matter_energy_moment_control.py','joint_singlet_common_mass_rg_results.json',
                    'round557_drafts/link_heat_kernel_probe.py',
                    'round557_drafts/link_heat_kernel_probe_results.json')
    return dict(round=557,tests_run=len(checks),failures=0,errors=0,checks=checks,
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in dependencies},
        geometry=geometry,link_heat_kernel=heat,pointer_kernel=pointer,witness_sources=rows,
        scope=dict(given_two_vertex_graph_and_SU2_are_inputs=True,
                   quartic_matter_inherited_from_round544=True,
                   local_Gauss_physical_sector_preserved_by_direct_invariant_pointer=True,
                   only_edge_configuration_record_uses_E1_not_full_momentum_menu=True,
                   pointer_and_ideal_pulse_are_additional_model_inputs=True,
                   source_energy_increase_is_not_total_apparatus_dissipation=True,
                   no_autonomous_detector_under_original_source_H_proved=True,
                   no_continuum_spacetime_or_gravity_generation=True))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args = parser.parse_args()
    result = run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result == json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=557,tests=result['tests_run'],
        energy_residual=max(r['independent_energy_identity_residual'] for r in result['witness_sources']),
        link_survival_t_half=result['link_heat_kernel']['rows'][0]['singlet_survival_for_constant_link_input'])))
