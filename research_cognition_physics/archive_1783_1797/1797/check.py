"""Finite neutral-scalar threshold in a common quark and lepton mass model."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np
from numpy.polynomial import Polynomial


HERE = Path(__file__).resolve().parent
CFG = json.loads((HERE/'inputs.json').read_bytes())
V = CFG['v_GeV']


def quadrature(order, extent):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    edges = np.arange(-extent, extent + 1, 4.0)
    centers = (edges[:-1] + edges[1:])/2
    return (centers[:, None] + 2*nodes).ravel(), np.tile(2*weights, len(centers))


def full_kernel(mass, base, splitting, order=10, extent=76):
    """Positive difference I(A)-I(H), including both same-line gluon insertions."""
    points, weights = quadrature(order, extent)
    u = np.exp(points)[None, :]
    t = 1/(1+np.exp(-points[:, None]))
    other = 1/(1+np.exp(points[:, None]))
    p = t*other
    base2 = base**2
    measure = (weights[:, None]*weights[None, :])*p*(u/(u+p))**2/2
    result = np.zeros((3, 3))
    for i, mi in enumerate(mass):
        for j in range(i+1):
            mj = mass[j]
            h = t*mi**2 + other*mj**2
            val = np.sum(measure*splitting/((h+u*base2)*(h+u*(base2+splitting))))
            result[i,j] = result[j,i] = val*mi*mj/(16*math.pi**2)
    return result


def kernel_tail(masses, base, splitting, extent):
    eta = math.exp(-extent)
    a = base**2
    b = a+splitting
    mi = masses[:,None]
    mj = masses[None,:]
    h0 = np.minimum(mi, mj)**2
    tail = splitting*eta/h0**2*(3.5+np.log1p(2*h0/(a*eta))) + splitting*eta**2/(4*a*b)
    return tail*mi*mj/(16*math.pi**2)


def mass_run(base):
    alpha = CFG['alpha_s_MZ']
    gain = 1.
    for start, end, nf in [(91.1876, 163., 5), (163., base, 6)]:
        beta = 11 - 2*nf/3
        final = alpha/(1 + alpha*beta/(2*math.pi)*math.log(end/start))
        gain *= (final/alpha)**(4/beta)
        alpha = final
    return gain


def common_data():
    previous = runpy.run_path(str(HERE/CFG['mass_source']))
    tensors, ckm = previous['sources']()
    generator = runpy.run_path(str(previous['OLD']['HERE']/previous['OLD']['CFG']['charged_background']))
    flow, _ = generator['charged_background'](.1)
    return tensors, ckm, *(np.array(flow['central'](species))/np.sqrt(2)
                      for species in ('down', 'charged_leptons'))


def all_vacua(poly):
    t = Polynomial([0., 1., -1.])
    square = (poly*Polynomial(poly.coef.conjugate())).coef.real
    normsq = Polynomial(square)
    slope = normsq.deriv()*t-2*normsq*t.deriv()
    extrema = [r.real for r in slope.roots() if abs(r.imag)<1e-8 and 0<r.real<1]
    witness = min(extrema, key=lambda z: abs(poly(z))/t(z))
    # Cover the entire interval, including possible extreme values near its ends.
    half = .5/32768
    mids = (np.arange(32768)+.5)/32768
    low = np.abs(poly(mids))-half*np.abs(poly.deriv()(mids))-half**2*abs(poly.coef[2])
    denom_high = t(mids)+half*np.abs(t.deriv()(mids))+half**2
    floor = float(min(low/denom_high))
    assert floor > 0 and floor <= abs(poly(witness))/t(witness)*(1+1e-12)
    return {'uniform_floor': floor, 'minimizing_s_squared': witness,
            'witness_amplitude': abs(poly(witness))/t(witness)}


def evaluate(d, el, yd, ye, kernel, leading):
    # kernel is positive I(A)-I(H); the overall minus is restored in the signed shift.
    qmass = np.diag(yd)
    poly = Polynomial([np.sum(kernel*d*d.T), -2*np.sum(kernel*d*qmass),
                       np.sum(kernel*qmass*qmass)])
    t = CFG['s']**2
    n = (d-t*qmass)/math.sqrt(t*(1-t))
    value = np.sum(kernel*n*n.T)
    np.testing.assert_allclose(poly(t)/(t*(1-t)), value, rtol=1e-12)
    lep = (el-t*np.diag(ye))/math.sqrt(t*(1-t))
    angle = np.angle(value) - 2*np.angle(lep[2,2])
    root_deg = (-math.degrees(angle)/2) % 90
    amplitude = abs(value)
    tolerance = math.degrees(math.asin(min(1, CFG['theta_budget']/amplitude)))/2
    return {'fixed_s': {'max_abs_strong_threshold': amplitude,
                        'kernel_C_real_imag': [value.real, value.imag],
                        'signed_threshold_at_phi_zero': -value.imag,
                        'heavy_tau_CP_zeros_deg_mod_180': [root_deg, root_deg+90],
                        'single_source_halfwidth_deg': tolerance,
                        'leading_power_relative_difference': abs(np.sum(leading*n*n.T)-value)/amplitude},
            'all_s': all_vacua(poly)}


def kaon_window(d, ckm, mass_a, splitting):
    prev = runpy.run_path(str(HERE.parent/'1796'/'check.py'))
    gamma = np.array([[-10., 1/6], [-40., 34/3]]).T + 16*np.eye(2)
    evals, basis = np.linalg.eig(gamma)
    weights = np.linalg.solve(basis, [1., 0.])
    ref = prev['M']['KA']['oscillation'](prev['M']['KC'])
    gain = mass_run(163.)
    central, _, err = prev['clock'](ckm)
    factor = 1e15 * ref['central_eps_response_per_im']
    obs = prev['EPS']['EPS_EXP']
    err_obs = prev['EPS']['EPS_ERROR']
    kaon_phase = ckm[0,1].conjugate()*ckm[0,0]
    common = kaon_phase.conjugate()/abs(kaon_phase)
    a, b = d[1,0]*common, d[0,1].conjugate()*common
    s = CFG['s']
    norm = gain**2/(s*s*(1-s*s))
    def sll_over_mass2(mass):
        tau = -math.log(mass_run(mass)/gain)/8.
        evolved = basis @ (np.exp(-evals*tau)*weights)
        return np.dot([-.064, -.107], evolved)/mass**2
    mass_h = math.hypot(mass_a, math.sqrt(splitting))
    ss = sll_over_mass2(mass_h)-sll_over_mass2(mass_a)
    const = -norm*.18*(a*b).conjugate()*(1/mass_h**2+1/mass_a**2)
    cos_part = -norm*.5*ss*(a*a+b*b).conjugate()
    sin_part = -norm*.5*ss*(1j*(a*a-b*b)).conjugate()
    amp_im = math.hypot(cos_part.imag, sin_part.imag)
    low, high = (factor*(const.imag+sign*amp_im) for sign in (-1,1))
    # A single SM epsilon within the old marginal box must fit every potential angle.
    common_SM = [max(central-2*err, obs-2*err_obs-low),
                 min(central+2*err, obs+2*err_obs-high)]
    kaon_input = prev['M']['KA']['read'](prev['M']['KC']['kaon'])
    dm, dm_err = kaon_input['ktev']['delta_m_1e6_per_second']
    clock_scale = 1e6*kaon_input['hbar_GeV_s']/2
    re_spread = math.hypot(cos_part.real,sin_part.real)
    real_interval = [(dm-2*dm_err)*clock_scale-const.real+re_spread,
                     (dm+2*dm_err)*clock_scale-const.real-re_spread]
    return {'epsilon_K_shift_range': [low, high], 'common_SM_epsilon_interval': common_SM
            if common_SM[0] <= common_SM[1] else None,
            'common_SM_M12_real_GeV_for_Delta_m_small_CP': real_interval,
            'some_phi_fits_SM_and_data_marginal_boxes':
                bool(central+2*err+high >= obs-2*err_obs and central-2*err+low <= obs+2*err_obs),
            'max_abs_new_real_GeV':abs(const.real)+math.hypot(cos_part.real,sin_part.real)}


def calculate():
    tensors, ckm, down, leptons = common_data()
    out = {}
    for base, quartic in CFG['mass_and_Z4_equal_abs_Z5']:
        splitting = quartic*V**2
        gain = mass_run(base)
        masses = down * V * gain
        k = full_kernel(masses, base, splitting)
        alt = full_kernel(masses, base, splitting, order=14, extent=84)
        np.testing.assert_allclose(k, alt, rtol=2e-8)
        ratios = masses[:,None]/masses[None,:]
        leading = (ratios + 1/ratios)*math.log1p(splitting/base**2)/(32*math.pi**2)
        assert np.max(abs(k/leading-1)) < .003
        result = {name: evaluate(tensors[name][1]*gain, tensors[name][2], down*gain,
                                 leptons, k, leading) for name in CFG['branches']}
        for name in CFG['branches']:
            result[name]['kaon_at_the_same_split'] = kaon_window(tensors[name][1], ckm, base, splitting)
            t = CFG['s']**2
            n = gain*(tensors[name][1]-t*np.diag(down))/math.sqrt(t*(1-t))
            result[name]['fixed_s']['integral_tail_bound'] = np.sum(
                kernel_tail(masses, base, splitting, 76)*abs(n*n.T))
        if base == 24000:
            assert all(case['kaon_at_the_same_split']['common_SM_epsilon_interval']
                       for case in result.values())
            assert result[CFG['branches'][0]]['fixed_s']['max_abs_strong_threshold'] > CFG['theta_budget']
        out[str(base)] = {'QCD_gain':gain,'m_H_minus_A_GeV':splitting/(math.hypot(base,math.sqrt(splitting))+base),
                          'kernel_refinement_max_fraction':np.max(abs(k/alt-1)),
                          'branches':result}
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args()
    data = json.dumps(calculate(), indent=2)+'\n'
    path = HERE/'results.json'
    if args.check:
        assert json.loads(data) == json.loads(path.read_bytes())
    else:
        path.write_bytes(data.encode('utf-8'))
    print('1797',hashlib.sha256(data.encode()).hexdigest()[:12])
