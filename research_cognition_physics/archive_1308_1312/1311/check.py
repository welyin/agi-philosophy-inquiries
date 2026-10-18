"""1311 joint true-kernel domain and negative-b branch degeneracy.

Imports the frozen 1310 kernel; no new apparatus construction is claimed.
Default/--check is read-only. --write-new creates results.json exclusively.
General-domain claims are proved in proof.md, not inferred from these grids.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
RESEARCH = STAGE.parent
ROOT = RESEARCH.parent
OUT = HERE / 'results.json'
KERNEL_PATH = STAGE / '1310' / 'check.py'
spec = importlib.util.spec_from_file_location('frozen_round1310_kernel', KERNEL_PATH)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def integral(x, w, kernel):
    den = kernel['zs'] * x + w
    assert np.all(den > 0)
    return float(kernel['weights'] @ (kernel['ks'] / den**2)) / 4


def fixed_w_root(w, target, kernel):
    # Equality at x=0 is excluded for a finite positive coupling.
    if target >= integral(0.0, w, kernel):
        return None
    c2 = float(kernel['weights'] @ (kernel['ks'] / kernel['zs']**2))
    lo, hi = 0.0, .5 * math.sqrt(c2 / target)
    for _ in range(90):
        mid = (lo + hi) / 2
        if integral(mid, w, kernel) > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def boundary_w(x, target, kernel):
    k0 = kernel['ktot']
    c2 = float(kernel['weights'] @ (kernel['ks'] / kernel['zs']**2))
    xc = .5 * math.sqrt(c2 / target)
    wc = .5 * math.sqrt(k0 / target)
    if x == 0:
        return wc
    if x >= xc:
        return 0.0 if math.isclose(x, xc, rel_tol=2e-15) else None
    lo, hi = 0.0, wc
    for _ in range(90):
        mid = (lo + hi) / 2
        if integral(x, mid, kernel) > target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def boundary_derivatives(x, w, kernel):
    den = kernel['zs'] * x + w
    a = kernel['weights'] * kernel['ks'] / den**3
    wp = -float(a @ kernel['zs']) / float(np.sum(a))
    wpp = 3 * float(np.sum(kernel['weights'] * kernel['ks']
                    * (kernel['zs'] + wp)**2 / den**4)) / float(np.sum(a))
    return wp, wpp


def beta(nc, ng):
    return Fraction(-43 + 2 * ng * (nc + 1), 6)


def lower_x(b, xf, length):
    return xf + max(float(b), 0.) * length / (8 * math.pi**2)


def calibrated_box():
    mu = 2.**-12
    vmin, vmax = .2, .3
    gy, gmax, anchor = .1, .3, .1
    rows = []
    for n in (3, 4, 5):
        neutral = mu * np.arange(1, n + 1)
        charged = mu * np.array([1, 2, 4, 6, 8][:n])
        quark_up = np.array([.0007, .0011, .0016, .0022, .0029][:n])
        quark_down = np.array([.0008, .0012, .0017, .0023, .0028][:n])
        masses = np.r_[neutral, charged, quark_up, quark_down]
        ymax = math.sqrt(2) * float(np.max(masses)) / vmin
        assert np.min(masses) > 0 and np.max(masses) < .003
        assert ymax < .022 and np.max(masses) < anchor
        rows.append({'ng': n, 'max_registered_fermion_mass': float(np.max(masses)),
                     'max_Yukawa_singular_value_over_v_box': ymax,
                     'neutrino_masses': neutral.tolist(),
                     'charged_masses': charged.tolist()})
    mw = gmax * vmax / 2
    mz = math.sqrt(gmax**2 + gy**2) * vmax / 2
    mh = math.sqrt(2 * .04) * vmax
    assert max(mw, mz, mh) < anchor
    phase = old.phase_compatibility()
    assert all(p['fixed_first_two_mass_columns'] for p in phase['branches'])
    lambda_poles = [{'lambda': lam, 'v': .25, 'm_h_squared': 2 * lam * .25**2}
                    for lam in (.02, .03, .04)]
    assert len({p['m_h_squared'] for p in lambda_poles}) == 3
    return {'mu': mu, 'v_closed_box': [vmin, vmax],
            'lambda_closed_box': [.02, .04], 'gY': gy,
            'mu0': anchor, 'mW_upper': mw, 'mZ_upper': mz, 'mH_upper': mh,
            'branches': rows, 'phase_structural_reuse': phase,
            'lambda_changes_physical_mass_at_fixed_task': lambda_poles,
            'all_registered_masses_below_anchor': True,
            'full_total_sources_equal_assumed': False,
            'actual_error_or_device_certified': False}


def compute():
    mu = 2.**-12
    kernel = old.kernel_setup(mu, 2 * mu, .1, order=128)
    dense = old.kernel_setup(mu, 2 * mu, .1, order=256)
    g0, v0, gmax, run_length = .2, .25, .3, 1.
    x0, w0, xf = g0**-2, v0*v0/4, gmax**-2
    sig0 = integral(x0, w0, kernel)
    target = sig0 / 2
    exposure = 2 * math.log(2) / sig0
    assert sig0 > 0 and math.isfinite(exposure)
    quad_rel = abs(sig0 - integral(x0, w0, dense)) / sig0
    assert quad_rel < 1e-10
    assert math.isclose(integral(x0, w0, kernel), old.sigma(g0, v0, kernel), rel_tol=2e-14)
    c2 = float(kernel['weights'] @ (kernel['ks'] / kernel['zs']**2))
    xc, wc = .5*math.sqrt(c2/target), .5*math.sqrt(kernel['ktot']/target)
    assert math.isclose(integral(xc, 0., kernel), target, rel_tol=2e-14)
    assert math.isclose(integral(0., wc, kernel), target, rel_tol=2e-14)
    assert boundary_w(xc*1.001, target, kernel) is None
    wminus, wplus = .2**2/4, .3**2/4
    xminus = fixed_w_root(wminus, target, kernel)
    assert xminus > xf
    Bminus = 8*math.pi**2 * (xminus-xf)/run_length
    boundary_rows = []
    max_fd1, max_fd2, max_level = 0., 0., 0.
    for fraction in (.03, .1, .25, .5, .75, .9, .97):
        x = xc * fraction
        w = boundary_w(x, target, kernel)
        assert w > 0
        wp, wpp = boundary_derivatives(x, w, kernel)
        assert -kernel['z1'] < wp < -kernel['z0'] and wpp > 0
        level_rel = abs(integral(x, w, kernel)/target-1)
        h = xc * 1e-4
        wm, ww = boundary_w(x-h, target, kernel), boundary_w(x+h, target, kernel)
        fd1 = (ww-wm)/(2*h)
        fd2 = (ww-2*w+wm)/(h*h)
        err1 = abs(fd1-wp)/abs(wp)
        err2 = abs(fd2-wpp)/wpp
        assert level_rel < 2e-14 and err1 < 1e-6 and err2 < 2e-5
        max_fd1, max_fd2, max_level = max(max_fd1,err1),max(max_fd2,err2),max(max_level,level_rel)
        boundary_rows.append({'x': x, 'w': w, 'Wprime': wp, 'Wsecond': wpp,
                              'strict_positive_curvature': True})
    # Jensen direction is for the convex BOUNDARY, not its hypograph.
    xa, xb = .2*xc, .8*xc
    wa, wb = boundary_w(xa,target,kernel), boundary_w(xb,target,kernel)
    wm = boundary_w((xa+xb)/2,target,kernel)
    assert wm < (wa+wb)/2
    # Three distinct negative coefficients; all direct grid predicates coincide.
    points = [(g,v) for g in (.08,.13,.19,.2,.21,.3,.34)
                        for v in (.19,.2,.24,.25,.26,.3,.31)]
    branch_rows = []
    signatures = []
    for n in (3,4,5):
        b = beta(3,n)
        assert b < 0 and lower_x(b,xf,run_length) == xf
        signature = []
        for g,v in points:
            x,w = g**-2,v*v/4
            trajectory = [x-float(b)*t/(8*math.pi**2) for t in np.linspace(0,run_length,31)]
            run_ok = min(trajectory) >= xf
            exact_run_ok = x >= lower_x(b,xf,run_length)
            assert run_ok == exact_run_ok
            domain = .2 <= v <= .3 and integral(x,w,kernel) >= target and run_ok
            signature.append(bool(domain))
        signatures.append(signature)
        branch_rows.append({'ng':n,'b_exact':str(b),'x_b':xf,
                            'grid_feasible':sum(signature),'grid_cases':len(points)})
    assert signatures[0] == signatures[1] == signatures[2]
    positive_rows = []
    for b in (Fraction(0),Fraction(5,6),Fraction(10),Fraction(100),Fraction(5000)):
        xlow = lower_x(b,xf,run_length)
        wlim = boundary_w(xlow,target,kernel)
        feasible = integral(xlow,wminus,kernel) >= target
        assert feasible == (xlow <= xminus)
        positive_rows.append({'b':str(b),'x_b':xlow,'positive_v_max_squared':None if wlim is None else 4*wlim,
                              'closed_box_nonempty':bool(feasible)})
    finite_limits = [r['positive_v_max_squared'] for r in positive_rows if r['positive_v_max_squared'] is not None]
    assert all(a>b for a,b in zip(finite_limits,finite_limits[1:]))
    # Exact endpoint semantics for closed v boxes, including x=0 exclusion.
    endpoint_target = integral(xf,wminus,kernel)
    assert integral(xf,wminus,kernel) == endpoint_target
    assert integral(xf*1.001,wminus,kernel) < endpoint_target
    assert integral(xf,wminus*1.001,kernel) < endpoint_target
    assert integral(xf,wminus,kernel) < 1.001*endpoint_target
    assert integral(xf,wminus,kernel) > .999*endpoint_target
    infinite_g_target = integral(0.,wminus,kernel)
    assert fixed_w_root(wminus,infinite_g_target,kernel) is None
    assert fixed_w_root(wminus,infinite_g_target*1.001,kernel) is None
    # Analytic rational rectangle margin, independent of quadrature accuracy.
    ratio_g = Fraction(20,19)**2
    ratio_v = Fraction(26,25)**2
    lower_factor = Fraction(19,20)**4
    assert ratio_v < ratio_g and lower_factor > Fraction(1,2)
    rectangle = []
    for g in (.19,.2,.21):
        for v in (.24,.25,.26):
            value = integral(g**-2,v*v/4,kernel)
            assert value/sig0 >= float(lower_factor)*(1-1e-14)
            probability = -math.expm1(-exposure*value)
            assert probability > .5 and g < gmax
            rectangle.append({'g':g,'v':v,'P':probability,'sigma_over_center':value/sig0})
    prob_lower = 1-2**(-2*float(lower_factor))
    assert prob_lower > .5
    assert math.isclose(-math.expm1(-exposure*sig0),.75,rel_tol=2e-15)
    # Exact integer shape is inherited; floats locate this declared sample's
    # boundary only. No sample value is chosen to single out three generations.
    kmax = math.floor((43+6*Bminus)/24)
    counts, integer_total = [], 0
    for kk in range(1,kmax+1):
        ngmax = math.floor((43+6*Bminus)/(8*kk))
        if ngmax < 3:
            continue
        for ng in range(3,ngmax+1):
            assert lower_x(beta(4*kk-1,ng),xf,run_length) <= xminus
        assert lower_x(beta(4*kk-1,ngmax+1),xf,run_length) > xminus
        integer_total += ngmax-2
        if kk <= 5:
            counts.append({'Nc':4*kk-1,'ng_min':3,'ng_max':ngmax})
    assert integer_total > 3
    files = [HERE/'proof.md',HERE/'check.py',HERE/'NEXT.md',STAGE/'research_note_1311.md',
             STAGE/'_admission/round1311_admission.md',STAGE/'1308/proof.md',
             STAGE/'1309/proof.md',STAGE/'1310/proof.md',KERNEL_PATH,
             STAGE/'1310/results.json',STAGE/'1310/acceptance.json',
             RESEARCH/'archive_1290_1307/1303/proof.md',RESEARCH/'archive_1290_1307/1307/proof.md']
    return {'schema':'round1311_joint_true_kernel_domain_v1',
            'status':'PASS_CONDITIONAL_DOMAIN_AND_BRANCH_DEGENERACY','new_axioms':0,
            'actual_device_certified':False,'all_parameters_classified':False,
            'kernel':{'from':'frozen round1310 module','m':mu,'n':2*mu,'E':.1,
                      'z0':float(kernel['z0']),'z1':float(kernel['z1']),
                      'K_integral':kernel['ktot'],'K_analytic_lower':float(kernel['klower']),
                      'center_quadrature_relative_difference':quad_rel},
            'task':{'g_center':g0,'v_center':v0,'sigma_center':sig0,'S_star':target,
                    'p_star':.5,'exposure_max':exposure,'center_P':.75,
                    'constructed_finite_witness_not_observed_budget':True,
                    'gmax':gmax,'Lrun':run_length},
            'domain':{'x_floor':xf,'x_c':xc,'w_c':wc,'x_at_v_lower':xminus,
                      'B_at_v_lower':Bminus,'boundary_samples':boundary_rows,
                      'max_level_relative_residual':max_level,
                      'max_derivative_fd_relative_difference':max_fd1,
                      'max_second_derivative_fd_relative_difference':max_fd2,
                      'boundary_strictly_convex_not_claiming_convex_domain':True,
                      'positive_b_contraction':positive_rows,
                      'closed_box_empty_singleton_interior_checked':True,
                      'infinite_g_and_vzero_excluded':True},
            'negative_b_degeneracy':{'branches':branch_rows,
                                     'grid_signatures_identical':True,
                                     'universal_claim_supported_by_analytic_eq17_not_grid':True},
            'strict_rectangle':{'exact_F_over_center_lower':str(lower_factor),
                                'analytic_P_lower':prob_lower,'samples':rectangle},
            'local_integer_sample':{'total_count':integer_total,'first_color_blocks':counts,
                                    'not_a_full_global_quantum_classification':True},
            'common_parameter_box':calibrated_box(),
            'source_sha256':{p.relative_to(ROOT).as_posix():sha256(p.read_bytes()).hexdigest() for p in files}}


def compare(actual, saved, path='root'):
    if isinstance(actual,dict):
        assert isinstance(saved,dict) and actual.keys()==saved.keys(),path
        for key in actual:compare(actual[key],saved[key],path+'.'+key)
    elif isinstance(actual,list):
        assert isinstance(saved,list) and len(actual)==len(saved),path
        for i,(a,b) in enumerate(zip(actual,saved)):compare(a,b,f'{path}[{i}]')
    elif isinstance(actual,float):
        assert isinstance(saved,(float,int)) and math.isclose(actual,saved,rel_tol=2e-10,abs_tol=1e-300),(path,actual,saved)
    else:
        assert actual==saved,(path,actual,saved)


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true')
    group.add_argument('--write-new',action='store_true')
    args=parser.parse_args()
    result=compute()
    if args.write_new:
        with OUT.open('x',encoding='utf-8',newline='\n') as handle:
            json.dump(result,handle,ensure_ascii=False,indent=2,allow_nan=False)
            handle.write('\n')
    else:
        compare(result,json.loads(OUT.read_text(encoding='utf-8')))
    print(json.dumps({'status':result['status'],
                      'same_three_four_five_domain':True,
                      'strict_rectangle_P_lower':result['strict_rectangle']['analytic_P_lower'],
                      'local_integer_candidates':result['local_integer_sample']['total_count'],
                      'bindings':len(result['source_sha256'])},ensure_ascii=False))


if __name__=='__main__':
    main()
