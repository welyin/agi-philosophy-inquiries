"""1169: joint original-clock / relative-CO prediction in one moving model.

The Newton Hamiltonian and R_CO actual measurement contract are inputs.
The common Gaussian reference is approximate; exact ordered mass fibers remain.
Default recomputes and compares read-only.  --write-once never overwrites.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
M, MU, G_SOURCE, SOFT, SIGMA = 1e8, 1e12, 10., 100., .01
T, CAPTURE = 32*math.pi+math.pi/2, math.pi/4
THRESHOLD, WIDTH = .032, .001
CONTRACT_ERROR, STAT_ERROR, NUM_ERROR = .003, .0005, .000002
SOURCES = {
 'archive_1106_1119/1112/check.py': '6e49d71ec5535fd6999432481938663a5543e2704db993ec7a00b8f475b6cd9f',
 'archive_935_955/research_note_941.md': 'a12e2dd67212a77538f850908efc6057e450d4def00dccd59aaa7cfaea346b7f',
 'archive_935_955/research_note_953.md': 'e510568c8812a86a6829be66d6733195858fd7e72c6686e14f060afce3bc2d88',
 'archive_956_989/research_note_978.md': '0f1481e7e8602a6980855a3ec1a3efe7c18b945d9d03b79b10a6258580a304be',
 'archive_956_989/research_note_985.md': 'd8b193fd13cae40c986381134d14c356ece5bb6dbaaaa5c0bf4ec8b0a5ec3e8a',
 'archive_1161_/1168/proof.md': 'c5674849ab457dd3b6a4bba2d0ce81b03263f94c1b3881303aa76821e46c7711',
 'archive_1161_/1168/acceptance.json': '98cf567ed554b0dfd4b20a066349c3a58808d16478dfaa608b214b577d09ef72',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def small(value, label, tol=3e-9):
    error = float(np.linalg.norm(value))
    assert error < tol, (label, error, tol)
    return error


def kernel(radius):
    s = radius*radius+SOFT*SOFT
    return (1/math.sqrt(s), -radius/s**1.5,
            np.array([(2*radius*radius-SOFT*SOFT)/s**2.5, -s**-1.5, -s**-1.5]))


def fg(kappa, t):
    """Stable power series, with kappa permitted to be negative."""
    x = kappa*t*t
    f = sum(x**n/math.factorial(2*n) for n in range(10))
    g = t*sum(x**n/math.factorial(2*n+1) for n in range(10))
    shift = t*t*sum(x**n/math.factorial(2*n+2) for n in range(9))
    return f, g, shift


def motion(radius, t, mass=M, coupling=G_SOURCE):
    """Complete two-body Gaussian, including source recoil.  No fixed Y."""
    u0, du, hess = kernel(radius)
    total, reduced = mass+MU, mass*MU/(mass+MU)
    initial = np.diag([SIGMA**2, SIGMA**2, 1/(4*SIGMA**2), 1/(4*SIGMA**2)])
    omega = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    means, covs, symplectic = [], [], []
    for axis in range(3):
        kappa = coupling*(1+mass/MU)*hess[axis]
        force = coupling*(1+mass/MU)*(du if axis == 0 else 0.)
        f, g, shift = fg(kappa, t)
        relative = np.array([-f, f, -g/mass, g/MU])
        velocity = np.array([-kappa*g, kappa*g, -f/mass, f/MU])
        com = np.array([mass/total, MU/total, t/total, t/total])
        momentum = np.array([0., 0., 1., 1.])
        matrix = np.array([com-MU/total*relative,
                           com+mass/total*relative,
                           mass/total*momentum-reduced*velocity,
                           MU/total*momentum+reduced*velocity])
        d = force*shift
        mean = np.array([-MU/total*d, mass/total*d, -reduced*force*g, reduced*force*g])
        means.append(mean)
        covs.append(matrix@initial@matrix.T)
        symplectic.append(small(matrix@omega@matrix.T-omega, 'complete symplectic flow', 2e-8))
    means, covs = np.array(means), np.array(covs)
    relvec = np.array([-1., 1., 0., 0.])
    relmean = means[:, 1]-means[:, 0]
    relvar = np.einsum('i,aij,j->a', relvec, covs, relvec)
    energy = float(np.sum((means[:, 2]**2+covs[:, 2, 2])/(2*mass)
                    +(means[:, 3]**2+covs[:, 3, 3])/(2*MU))
                    -coupling*mass*(u0+du*relmean[0]+.5*np.sum(hess*(relmean**2+relvar))))
    return dict(mean=means, covariance=covs, relative_mean=relmean,
                relative_variance=relvar, max_symplectic_error=max(symplectic), energy=energy)


def clock_time(radius):
    """Closed integrals on the common moving reference, not the exact clock."""
    u0, du, hess = kernel(radius)
    kappa = G_SOURCE*(1+M/MU)*hess[0]
    force = G_SOURCE*(1+M/MU)*du
    z = math.sqrt(kappa)*T
    amp = force/kappa
    # Stable expansions avoid cancellation of sinhc(z)-1 at small z.
    s1 = sum(z**(2*n)/math.factorial(2*n+1) for n in range(1, 10))
    s2 = sum((.5*2**(2*n)-2)*z**(2*n)/math.factorial(2*n+1) for n in range(2, 10))
    sv = .5*sum((2*z)**(2*n)/math.factorial(2*n+1) for n in range(1, 10))
    irel = amp*T*s1
    irel2 = amp*amp*T*s2
    ivx2 = (MU/(M+MU))**2*force**2/kappa*T*sv
    tau = T-G_SOURCE*(u0*T+du*irel+.5*hess[0]*irel2)-.5*ivx2
    next1 = z**20/math.factorial(21)
    next2 = (.5*2**20-2)*z**20/math.factorial(21)
    nextv = .5*(2*z)**20/math.factorial(21)
    tail = 2*(G_SOURCE*(abs(du*amp*T)*next1+.5*abs(hess[0]*amp*amp*T)*next2)
              +.5*(MU/(M+MU))**2*force**2/kappa*T*nextv)
    assert tail < 1e-30
    nodes, weights = np.polynomial.legendre.leggauss(24)
    rates = []
    for node in nodes:
        t = T*(node+1)/2
        _, gg, shift = fg(kappa, t)
        eta, px = force*shift, -(M*MU/(M+MU))*force*gg
        rates.append(1-G_SOURCE*(u0+du*eta+.5*hess[0]*eta*eta)-px*px/(2*M*M))
    quadrature = T/2*float(weights@rates)
    return dict(tau=tau, age=CAPTURE+tau, age_mod_2pi=(CAPTURE+tau)%(2*math.pi),
                integrated_relative_mean=irel, integrated_relative_square=irel2,
                integrated_vx_square=ivx2, series_tail_bound=tail,
                independent_quadrature_error=small(tau-quadrature, 'clock phase integral', 2e-11))


def uniform_bounds(radius):
    """Conservative whole-time prepared-state bounds, not infinite H norm."""
    d = .052 if radius == 1000 else .014
    s, sp, p = .0143, 52., (105000. if radius == 1000 else 27000.)
    u0, du, hess = kernel(radius)
    normh = float(max(abs(hess)))
    # Whole interval follows cosh/sinh monotonicity; negative axes dominated
    # by the positive radial majorant.  This endpoint tests that majorant.
    kappa = G_SOURCE*(1+M/MU)*hess[0]
    f, g, _ = fg(kappa, T)
    radial_var_bound = f*f*2*SIGMA**2+g*g*(1/M**2+1/MU**2)/(4*SIGMA**2)
    baseline = motion(radius, T)
    assert abs(baseline['relative_mean'][0]) < d
    assert math.sqrt(radial_var_bound) < s
    # PX linear coefficients; transverse kappa magnitude <= radial kappa.
    reduced = M*MU/(M+MU)
    pvar_bound = 2*SIGMA**2*(reduced*kappa*g)**2 + (
        (M/(M+MU)+reduced*f/M)**2+(M/(M+MU)-reduced*f/MU)**2)/(4*SIGMA**2)
    assert math.sqrt(pvar_bound) < sp
    assert abs(baseline['mean'][0, 2]) < p
    moment3 = (d+105**(1/6)*s)**3
    soft = 4*T*G_SOURCE*(M+3)*moment3/(radius-1)**4
    log_tail_probability = math.log(6)-(1-d)**2/(6*s*s)
    fourth_r = d+15**.25*s
    eighth_r = d+945**.125*s
    l4_tail_polynomial = 1/SOFT+u0+abs(du)*fourth_r+.5*normh*eighth_r**2
    log_tail_error = math.log(T*G_SOURCE*(M+3)*l4_tail_polynomial)+log_tail_probability/4
    assert log_tail_error < math.log(1e-40)
    tail = 1e-40  # explicit positive upward bound; never floating zero
    grad = abs(du)+normh*d
    potential = G_SOURCE*(grad*s+.5*normh*math.sqrt(15)*s*s)
    kinetic = (2*p*sp+math.sqrt(15)*sp*sp)/(2*M*M)
    mass_denominator = 3*p*p/(2*M**3)
    clock = 3*T*(potential+kinetic+mass_denominator)
    total = soft+tail+clock
    return dict(relative_mean_bound=d, per_axis_sd_bound=s, momentum_mean_bound=p,
                momentum_sd_bound=sp, radial_sd_majorant=math.sqrt(radial_var_bound),
                momentum_sd_majorant=math.sqrt(pvar_bound), soft_remainder_bound=soft,
                log_tail_probability=log_tail_probability, log_tail_error=log_tail_error,
                positive_tail_bound=tail, clock_potential_rate_bound=potential,
                clock_kinetic_rate_bound=kinetic, clock_inverse_mass_rate_bound=mass_denominator,
                clock_shape_bound=clock, full_prepared_isometry_error=total)


def taylor_calibration(radius, order):
    """Independent finite Gaussian quadrature, only calibration of analytic bound."""
    mot = motion(radius, T)
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    grid = np.stack(np.meshgrid(nodes, nodes, nodes, indexing='ij'), axis=-1).reshape(-1, 3)
    ww = np.einsum('i,j,k->ijk', weights, weights, weights).reshape(-1)/math.pi**1.5
    xi = (mot['relative_mean']+grid*np.sqrt(2*mot['relative_variance'])).astype(np.longdouble)
    rr, aa = np.longdouble(radius), np.longdouble(SOFT)
    s0 = rr*rr+aa*aa
    u0 = 1/np.sqrt(s0)
    z = (2*rr*xi[:, 0]+np.sum(xi*xi, axis=1))/s0
    diff = u0*np.expm1(-np.log1p(z)/2)
    du = -rr/s0**np.longdouble(1.5)
    hx = (2*rr*rr-aa*aa)/s0**np.longdouble(2.5)
    ht = -1/s0**np.longdouble(1.5)
    remainder = diff-du*xi[:, 0]-.5*(hx*xi[:, 0]**2+ht*np.sum(xi[:, 1:]**2, axis=1))
    return T*G_SOURCE*(M+3)*math.sqrt(float(ww@np.asarray(remainder*remainder, float)))


def cdf(x):
    return .5*math.erfc(-x/math.sqrt(2))


def independent_decimal_prediction(radius):
    """Independent 80-digit series evaluation of the closed reference formula."""
    with localcontext() as context:
        context.prec = 80
        D = Decimal
        pi = D('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998')
        m, mu, rr, aa, sig = D('1e8'), D('1e12'), D(radius), D(100), D('.01')
        tt = D('32.5')*pi
        base = rr*rr+aa*aa
        uu, du = 1/base.sqrt(), -rr/(base*base.sqrt())
        hx = (2*rr*rr-aa*aa)/(base*base*base.sqrt())
        kappa, force = 10*(1+m/mu)*hx, 10*(1+m/mu)*du
        xx = kappa*tt*tt
        f = sum(xx**n/D(math.factorial(2*n)) for n in range(15))
        gg = tt*sum(xx**n/D(math.factorial(2*n+1)) for n in range(15))
        shift = tt*tt*sum(xx**n/D(math.factorial(2*n+2)) for n in range(15))
        displacement = -force*shift
        variance = f*f*2*sig*sig+gg*gg*(1/m**2+1/mu**2)/(4*sig*sig)
        amp = force/kappa
        i1 = amp*tt*sum(xx**n/D(math.factorial(2*n+1)) for n in range(1, 15))
        i2 = amp*amp*tt*sum((D(2)**(2*n-1)-2)*xx**n/D(math.factorial(2*n+1)) for n in range(2, 15))
        iv = (mu/(m+mu))**2*force**2/kappa*tt*sum(D(2)**(2*n-1)*xx**n/D(math.factorial(2*n+1)) for n in range(1, 15))
        tau = tt-10*(uu*tt+du*i1+hx*i2/2)-iv/2
        angle = (pi/4+tau-32*pi)/2
        sine = sum((-1)**n*angle**(2*n+1)/D(math.factorial(2*n+1)) for n in range(25))
        z = (displacement-D('.032'))/(variance+D('.001')**2).sqrt()
        cdf_series = D('.5')+sum((-1)**n*z**(2*n+1)/(D(2)**n*D(math.factorial(n))*D(2*n+1)) for n in range(36))/(2*pi).sqrt()
        cdf_tail = abs(z**73/(D(2)**36*D(math.factorial(36))*73))/(2*pi).sqrt()
        sin_tail = abs(angle**51/D(math.factorial(51)))
        assert cdf_tail < D('1e-40') and sin_tail < D('1e-60')
        return dict(reference_joint=float(cdf_series*sine**6),
                    cdf_next_term_bound=float(cdf_tail), sin_next_term_bound=float(sin_tail),
                    decimal_precision=80)


def old_clock_group(phases):
    spec = importlib.util.spec_from_file_location('round1112_source_1169', BASE/'archive_1106_1119/1112/check.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    material = module.Material(1)
    eigen, vec = np.linalg.eigh(material.hfk)
    small(eigen-np.repeat(np.arange(4), 72), 'full 288D spectrum')
    small(material.hfk-material.gamma@np.kron(material.hc, np.eye(72))@material.gamma.conj().T,
          'original full Gamma dressing')
    ready = material.encode_joint(np.eye(3, dtype=complex)[:, None, :], 0).reshape(288, 3)
    mask = np.zeros((4, 8, 3, 3))
    mask[:, [1, 3, 5, 7], :, :] = 1
    mask = mask.reshape(288)
    output = []
    for age in [CAPTURE, *(p['age'] for p in phases)]:
        evolved = (vec*np.exp(-1j*eigen*age))@(vec.conj().T@ready)
        effect = evolved.conj().T@(mask[:, None]*evolved)
        prob = math.sin(age/2)**6
        output.append(dict(age=age, mr_probability=prob,
            input_isometry_error=small(evolved.conj().T@evolved-np.eye(3), 'all unknown qutrit columns'),
            full_mr_effect_error=small(effect-prob*np.eye(3), 'same physical Mr on all columns')))
    # The exact fibers do not collapse into one scalar time: coefficients of
    # [PX^2,u] and [PY^2,u] in [H_0,H_3], evaluated without cancellation.
    px_comm = -G_SOURCE*3*(2*M+3)/(2*M*(M+3))
    py_comm = -G_SOURCE*3/(2*MU)
    assert px_comm != 0 and py_comm != 0
    return dict(dimension=288, new_positive_mass_baseline=M, old_internal_h_unchanged=True,
                full_effect_checks=output, ordered_fiber_commutator_px_coefficient=px_comm,
                ordered_fiber_commutator_py_coefficient=py_comm)


def calculate():
    for name, expected in SOURCES.items():
        assert digest(BASE/name) == expected, name
    radii = [1000, 2000]
    phases = [clock_time(r) for r in radii]
    old = old_clock_group(phases)
    movement, certificates, predictions = [], [], []
    for radius, phase in zip(radii, phases):
        fiber_rows = []
        for j in range(4):
            mot = motion(radius, T, M+j)
            initial = motion(radius, 0., M+j)
            energy_error = abs(mot['energy']-initial['energy'])
            small(energy_error, 'autonomous quadratic energy including recoil', 3e-8)
            momentum_error = small(mot['mean'][:, 2]+mot['mean'][:, 3], 'source opposite momentum')
            fiber_rows.append(dict(j=j, mean_shortening=-float(mot['relative_mean'][0]),
                relative_sd=math.sqrt(mot['relative_variance'][0]),
                mean_x=float(mot['mean'][0, 0]), source_shift=float(mot['mean'][0, 1]),
                px=float(mot['mean'][0, 2]), py=float(mot['mean'][0, 3]),
                total_momentum_error=momentum_error, quadratic_energy_error=energy_error,
                symplectic_error=mot['max_symplectic_error']))
        movement.append(dict(radius=radius, complete_mass_fibers=fiber_rows, phase=phase))
        bound = uniform_bounds(radius)
        cal9, cal15 = taylor_calibration(radius, 9), taylor_calibration(radius, 15)
        assert cal9 < bound['soft_remainder_bound'] and cal15 < bound['soft_remainder_bound']
        assert abs(cal9-cal15) < 1e-8
        certificates.append(dict(radius=radius, bounds=bound,
            independent_final_gaussian_remainder_9=cal9,
            independent_final_gaussian_remainder_15=cal15))
        mot = motion(radius, T)
        displacement = -float(mot['relative_mean'][0])
        std = math.sqrt(float(mot['relative_variance'][0])+WIDTH*WIDTH)
        position = cdf((displacement-THRESHOLD)/std)
        clock = math.sin(phase['age']/2)**6
        joint = position*clock  # reference ONLY; transferred with whole-state bound
        decimal_check = independent_decimal_prediction(radius)
        decimal_check['double_precision_difference'] = small(joint-decimal_check['reference_joint'], 'independent decimal reference', 1e-11)
        error = bound['full_prepared_isometry_error']+CONTRACT_ERROR+STAT_ERROR+NUM_ERROR
        assert error < .004
        predictions.append(dict(radius=radius, shortening=displacement,
            effective_cdf_sd=std, reference_position_probability=position,
            reference_clock_probability=clock, reference_joint_probability=joint,
            independent_decimal_check=decimal_check,
            total_per_preparation_error=error,
            conditional_joint_interval=[max(0., joint-error), min(1., joint+error)]))
    delta = predictions[0]['reference_joint_probability']-predictions[1]['reference_joint_probability']
    residual = sum(c['bounds']['full_prepared_isometry_error'] for c in certificates)
    lower = delta-sum(row['total_per_preparation_error'] for row in predictions)
    assert lower > .02
    max_contract = (delta-.02-residual-2*(STAT_ERROR+NUM_ERROR))/2
    # Pre-registered zero-coupling control. Same H and same readout for both
    # preparations; only physical initial translations differ, then subtract R_j.
    free = motion(1000, T, coupling=0.)
    free_sd = math.sqrt(free['relative_variance'][0]+WIDTH*WIDTH)
    free_clock = math.sin((CAPTURE+T)/2)**6
    free_joint = free_clock*cdf(-THRESHOLD/free_sd)
    # Ideal relative Lueders energy cost is not the actual local instrument cost.
    ideal_energy_bound = (1/M+1/MU)/(8*WIDTH**2)
    local_instrument_energy_bound = (1/M+1/MU)/(4*WIDTH**2)
    record_grid = 1e-7
    local_quantization_effect_bound = 2*record_grid/(math.sqrt(2*math.pi)*WIDTH)
    assert local_quantization_effect_bound < CONTRACT_ERROR
    confidence_failure = .01
    repeats = math.ceil(math.log(4/confidence_failure)/(2*STAT_ERROR**2))
    # Direct finite Fisher inequality calibration for Gaussian CDF effect.
    fisher_ratios = []
    for z in np.linspace(-7., 7., 57):
        f = cdf(float(z)); density = math.exp(-z*z/2)/math.sqrt(2*math.pi)
        fisher_ratios.append(density*density/(f*(1-f)))
    assert max(fisher_ratios) < 1.
    for name, expected in SOURCES.items():
        assert digest(BASE/name) == expected, 'source changed: '+name
    return dict(schema='research_round_1169_v1', round=1169, scientific_groups=4,
        sources_sha256=SOURCES, code_sha256=digest(Path(__file__)),
        parameters=dict(m=M, mu=MU, G=G_SOURCE/MU, G_mu=G_SOURCE, softening=SOFT,
            sigma_per_body_axis=SIGMA, wait=T, capture_age=CAPTURE,
            relative_threshold=THRESHOLD, cdf_width=WIDTH),
        group1_complete_original_clock=old,
        group2_moving_two_body_reference=movement,
        group3_all_unknown_input_error_certificates=certificates,
        group4_conditional_actual_joint_prediction=dict(rows=predictions,
            reference_joint_contrast=delta, certified_contrast_lower=lower,
            required_contrast=.02, contract_error_condition=CONTRACT_ERROR,
            max_equal_contract_error_for_strict_threshold=max_contract,
            statistical_error_each=STAT_ERROR, numerical_allowance_each=NUM_ERROR,
            joint_failure_probability=confidence_failure, independent_repeats_per_preparation=repeats,
            zero_coupling_control_joint_probability_each=free_joint,
            zero_coupling_control_contrast=0., relative_ideal_lueders_energy_bound=ideal_energy_bound,
            local_gaussian_records_position_energy_bound=local_instrument_energy_bound,
            finite_grid_spacing=record_grid, finite_grid_effect_error_bound=local_quantization_effect_bound,
            gaussian_binary_fisher_max_sample=max(fisher_ratios),
            readout_contract_is_input_not_implemented_hardware=True,
            actual_local_record_poststate_not_identified_with_relative_lueders=True),
        scope=dict(new_adopted_axioms=0, empirical_calibrations=0, gr_hard_deliveries=0,
            actual_ordered_echo_retained=True, common_moving_gaussian_is_approximation=True,
            new_private_remote_signal_claim=False, post_newton_error_certified=False))


def compare(a, b, path='root'):
    if isinstance(a, dict):
        assert isinstance(b, dict) and a.keys() == b.keys(), path
        for key in a:
            compare(a[key], b[key], path+'.'+key)
    elif isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), path
        for i, (left, right) in enumerate(zip(a, b)):
            compare(left, right, path+f'[{i}]')
    elif isinstance(a, float):
        assert math.isfinite(a) and math.isfinite(b) and math.isclose(a, b, rel_tol=2e-9, abs_tol=3e-10), (path, a, b)
    else:
        assert a == b, (path, a, b)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--write-once', action='store_true')
    mode.add_argument('--preview', action='store_true')
    args = parser.parse_args()
    result = calculate()
    target = HERE/'results.json'
    if args.write_once:
        with target.open('x', encoding='utf-8', newline='\n') as stream:
            json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
            stream.write('\n')
        print('1169 WRITE_ONCE PASS: 4 groups; finite R_CO contract remains input')
    elif args.preview:
        print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    else:
        before = digest(target)
        compare(result, json.loads(target.read_text(encoding='utf-8')))
        assert before == digest(target), 'default must not write results'
        print('1169 READ_ONLY PASS: 4 groups; saved results and all 7 sources unchanged')
        print(json.dumps(result['group4_conditional_actual_joint_prediction'], ensure_ascii=False, indent=2))
