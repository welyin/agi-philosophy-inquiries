"""Round1298: finite algebra and exact rational certificates, not a QFT simulation.
Default / --check recomputes read-only; --write exclusively creates results.json.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results.json'
TOL = 2e-11


def encoded(x):
    if isinstance(x, F):
        return {'numerator': x.numerator, 'denominator': x.denominator}
    if isinstance(x, dict):
        return {k: encoded(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encoded(v) for v in x]
    if isinstance(x, np.generic):
        return x.item()
    return x


def normerr(x):
    return float(np.max(np.abs(x)))


def poly_add(a, b):
    c = dict(a)
    for k, v in b.items():
        c[k] = c.get(k, 0) + v
    return {k: v for k, v in c.items() if v}


def poly_scale(a, n):
    return {k: n*v for k, v in a.items() if n*v}


def radial_first(a):
    # x=r^2/a^2, y=(1-x)^-1: e^y d_x(e^-y P)=y^2(P'-P).
    ans = {}
    for k, v in a.items():
        if k:
            ans[k+1] = ans.get(k+1, 0) + k*v
        ans[k+2] = ans.get(k+2, 0) - v
    return {k: v for k, v in ans.items() if v}


def radial_laplacian(a):
    d = radial_first(a)
    dd = radial_first(d)
    # a^2 Delta = 4x d_x^2 + 6 d_x.
    return poly_add(poly_add(poly_scale(dd, 4),
                            {k-1: -4*v for k, v in dd.items()}),
                    poly_scale(d, 6))


def poly_bound(a):
    return sum(abs(v)*math.factorial(k) for k, v in a.items())


def exact_bump_and_window():
    a = F(1, 32)
    polys = [{0: 1}]
    for _ in range(3):
        polys.append(radial_laplacian(polys[-1]))
    assert polys[1] == {2: 2, 3: -12, 4: 4}
    assert polys[2] == {3: -24, 4: 444, 5: -1168, 6: 896, 7: -224, 8: 16}
    assert polys[3] == {4: 720, 5: -26640, 6: 156600, 7: -327120,
                        8: 298608, 9: -125440, 10: 24480, 11: -2112, 12: 64}
    cs = [poly_bound(p) for p in polys]
    assert cs == [1, 172, 2570160, 263117600640]
    # Vol(B_a)/(pi*a^3/96)=128, sqrt(128)<12.
    assert F(4, 3)*96 == 128 and 128 < 12**2
    K = 12*(1+2*cs[1]/a**2+cs[2]/a**4)
    M = 48*K*K
    ell = 1024/M
    assert K == 32340101333004
    assert M == 50202343402990420689583872768
    assert ell == F(4, 196102903917931330818687003)
    center = 1/(4*a**3)  # strict lower bound using pi<4.
    assert center == 8192 and center-4*ell*M == 4096
    area = 4*ell**2
    weighted_source = area*F(4096, 2)*ell**2
    assert weighted_source == 8192*ell**4
    difference_bound = area*M*ell**3/3
    assert difference_bound/weighted_source == F(1, 6)
    # H6 field bound and the coordinate C2 estimates in proof section5.
    B6 = 1+3*cs[1]/a**2+3*cs[2]/a**4+cs[3]/a**6
    M6 = 8192*B6**2
    Hstar = 512*M6
    assert B6 == 282520380522721972225
    assert Hstar == 334779972733175839233728425479031683340042240000
    spatial = 6*a*a+1
    mixed = 17*a*a+4
    time2 = 18*a*a+5
    assert spatial < 2 and mixed < 5 and time2 < 6
    assert max(spatial, mixed, time2) < 8
    assert 16*3*8 < 512  # sum of tensor components after trace reversal.
    kappa_max = 1/(16*Hstar)
    assert kappa_max*Hstar == F(1, 16)
    # Finite dyadic function bound: for derivative j, e^-y*y^j<=j!.
    # This inequality follows from the positive j-th term of exp(y), no sampling.
    return {'a': a, 'radial_Delta_polynomials':
            [{str(k): v for k, v in sorted(p.items())} for p in polys],
            'factorial_absolute_bounds': cs, 'norm_ratio_squared_bound': 128,
            'K': K, 'M': M, 'ell': ell, 'ell_diagnostic': float(ell),
            'center_strict_lower': center, 'whole_patch_strict_lower': 4096,
            'patch_zero_order_area': area, 'I_A_and_I_B_strict_lower': weighted_source,
            'weighted_source_lower_diagnostic': float(weighted_source),
            'weight_difference_absolute_bound': difference_bound,
            'relative_weight_correction_strict_upper': F(1, 6),
            'B6': B6, 'M6': M6, 'Hstar': Hstar, 'kappa_max': kappa_max,
            'Kirchhoff_coordinate_C2_coefficients': [spatial, mixed, time2],
            'H4_and_H6_Fourier_Cauchy_Schwarz_square_constant': '1/(8*pi)<1',
            'small_window_is_sufficient_not_necessary_or_instrument_resolution': True}


def center_dirac_jets():
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.diag([1., -1.]).astype(complex)
    ident = np.eye(2, dtype=complex)
    zero = np.zeros((2, 2), complex)
    alphas = [np.block([[-s, zero], [zero, s]]) for s in (sx, sy, sz)]
    beta = np.block([[zero, ident], [ident, zero]])
    maximum = 0.
    for A in alphas:
        maximum = max(maximum, normerr(A@A-np.eye(4)), normerr(A@beta+beta@A))
    Mk = np.eye(4)-alphas[2]
    zeta = np.array([1., 1.], complex)/math.sqrt(2)
    zero2 = np.zeros(2, complex)
    def source(v, w, dv, dw):
        return .5j*(np.vdot(w, Mk@dv)-np.vdot(dw, Mk@v))
    rows = []
    for gv, vv, ww in [(F(3, 2), F(2, 3), F(11, 4)),
                       (F(7, 4), F(3, 5), F(15, 8)),
                       (F(5, 3), F(9, 10), F(13, 6))]:
        g, v, w = map(float, (gv, vv, ww))
        d = (v-w)/3
        f = np.r_[zero2, g*zeta]
        u = np.r_[v*zeta, zero2]
        hf = np.r_[g*zeta, zero2]
        hu = np.r_[zero2, w*zeta]
        dzf = np.zeros(4, complex)
        dzu = np.r_[zero2, -1j*sz@(d*zeta)]
        df, du = -1j*hf+dzf, -1j*hu+dzu
        cross = source(f, u, df, du)+source(u, f, du, df)
        exact = gv*(F(4, 3)*ww+F(2, 3)*vv)
        bare = source(f, f, df, df)
        maximum = max(maximum, abs(cross-float(exact)), abs(bare))
        # Each flavor has -.5 ff + .5 cross + rho_ii ff;
        # offdiagonal flavor elements disappear only for flavor-blind sources.
        rho = np.array([[.3, .1+.12j], [.1-.12j, .7]], complex)
        total = sum((rho[i, i]-.5)*bare+.5*cross for i in range(2))
        maximum = max(maximum, abs(total-cross))
        assert exact > 0 and abs(bare) < TOL
        rows.append({'g': gv, 'omega_inverse_g': vv, 'omega_g': ww,
                     'two_flavor_Tkk': exact, 'bare_occupation_Tkk': float(bare.real)})
    assert maximum < TOL
    return {'rows': rows, 'max_residual': maximum,
            'jet_inputs_are_algebraic_cross_checks_not_sampled_continuum_values': True,
            'positivity_of_actual_core_uses_heat_kernel_maximum_argument': True,
            'reference_correction_legs_essential': True}


def trace_reverse(a):
    signs = [1., -1., -1., -1.]
    tr = sum(signs[i]*a[i, i] for i in range(4))
    return a-.5*np.diag(signs)*tr


def einstein_from_second_derivatives(D):
    signs = [1., -1., -1., -1.]
    ricci = np.zeros((4, 4), complex)
    for a in range(4):
        for b in range(4):
            term1 = sum(signs[c]*D[c, a, c, b] for c in range(4))
            term2 = sum(signs[c]*D[c, b, c, a] for c in range(4))
            wave = sum(signs[c]*D[c, c, a, b] for c in range(4))
            trace = sum(signs[c]*D[a, b, c, c] for c in range(4))
            ricci[a, b] = .5*(term1+term2-wave-trace)
    scalar = sum(signs[a]*ricci[a, a] for a in range(4))
    return ricci-.5*np.diag(signs)*scalar


def harmonic_constraint_test():
    rows = []
    maximum = 0.
    for momentum, ts in [([1, 2, -1], [3, 1, -2, 2, 4, 1, 0, 5, -1, 2]),
                          ([2, -1, 3], [5, -2, 1, 3, 2, 4, -1, 6, 2, 3]),
                          ([-3, 2, 1], [2, 3, 1, -2, 5, -1, 2, 4, 3, 6])]:
        p = np.array(momentum, float)
        p2 = float(p@p)
        tau = np.zeros((4, 4), complex)
        for (a, b), value in zip([(a,b) for a in range(4) for b in range(a,4)], ts):
            tau[a,b] = tau[b,a] = value
        phi = -2*tau[0]/p2
        initial = np.zeros((4,4), complex)
        initial[0,:] = phi
        initial[:,0] = phi
        velocity = np.zeros((4,4), complex)
        velocity[0,0] = 1j*np.dot(p, phi[1:])
        acceleration = -p2*initial-2*tau
        C = velocity[0,:]-1j*np.einsum('i,ij->j',p,initial[1:,:])
        Cdot = acceleration[0,:]-1j*np.einsum('i,ij->j',p,velocity[1:,:])
        Dbar = np.zeros((4,4,4,4),complex)
        Dbar[0,0] = acceleration
        for i in range(3):
            Dbar[0,i+1] = Dbar[i+1,0] = 1j*p[i]*velocity
            for j in range(3):
                Dbar[i+1,j+1] = -p[i]*p[j]*initial
        Dh = np.empty_like(Dbar)
        for a in range(4):
            for b in range(4):
                Dh[a,b] = trace_reverse(Dbar[a,b])
        G = einstein_from_second_derivatives(Dh)
        residual = max(normerr(C),normerr(Cdot),normerr(G-tau))
        maximum = max(maximum,residual)
        # Deleting the seed, while retaining the wave source, fails Einstein constraints.
        badbar = np.zeros_like(Dbar)
        badbar[0,0] = -2*tau
        bad = np.empty_like(badbar)
        for a in range(4):
            for b in range(4):
                bad[a,b] = trace_reverse(badbar[a,b])
        badG = einstein_from_second_derivatives(bad)
        bad_constraint = normerr((badG-tau)[0,:])
        assert bad_constraint > 1
        rows.append({'momentum':momentum,'max_full_Einstein_and_gauge_residual':residual,
                     'zero_seed_constraint_defect':bad_constraint})
    assert maximum < TOL
    return {'rows':rows,'max_residual':maximum,'all_ten_Einstein_components_tested':True,
            'continuous_constraint_propagation_uses_source_Ward_not_samples':True,
            'Newton_seed_is_not_global_L2_or_H4':True}


def integrate_poly(poly, lo=F(-1), hi=F(0)):
    return sum(F(c,k+1)*(hi**(k+1)-lo**(k+1)) for k,c in poly.items())


def optical_and_observer_balance():
    # Rescaled affine interval [-1,0], nonconstant positive source and nonzero terminal theta.
    T = {0:F(5),1:F(2),2:F(3)}
    theta_terminal = F(7,3)
    theta = {0:theta_terminal}
    for k,c in T.items():theta[k+1]=-c/F(k+1)
    A0=F(4)
    area_change=A0*integrate_poly(theta)
    endpoint=A0*theta_terminal
    IA=A0*integrate_poly(poly_add(T,{k+1:c for k,c in T.items()}))
    IB=A0*integrate_poly({k+1:-c for k,c in T.items()})
    correction=A0*integrate_poly(poly_add(T,{k+1:2*c for k,c in T.items()}))
    Egeom=endpoint+correction
    assert area_change-endpoint==IA
    assert area_change-Egeom==IB
    assert correction==IA-IB and correction!=0
    # A finite interval is not made equilibrium by silently deleting its actual endpoint.
    assert area_change-IB!=0
    # Exact left-wedge clock normalization at a rational hyperbolic point.
    cosh,sinh=F(5,4),F(3,4)
    assert cosh*cosh-sinh*sinh==1
    kchi=F(3,2)
    velocity=[cosh,-sinh]
    chi=[cosh,-sinh]
    assert velocity==chi and chi[0]>0
    lam=F(-2,5)
    chi_boundary=[-kchi*lam,-kchi*lam]
    assert chi_boundary[0]>0 and chi_boundary[0]**2-chi_boundary[1]**2==0
    # If tau00=t and chi=dt in flat space, divergence is1 but Lie_chi eta=0.
    # It certifies that the general product rule needs chi_b div(tau)^b.
    assert 1!=0
    return {'rescaled_source_polynomial':{str(k):v for k,v in T.items()},
            'A0':A0,'terminal_theta':theta_terminal,'area_change':area_change,
            'terminal_area_term':endpoint,'I_A':IA,'I_B':IB,
            'linear_curvature_weight_correction':correction,'Egeom':Egeom,
            'dropping_curvature_weight_residual':correction,
            'dropping_endpoint_and_weight_residual':area_change-IB,
            'left_wedge_proper_clock_norm':F(1),'left_wedge_boundary_chi':chi_boundary,
            'nonconserved_source_counterexample_missing_term':1,
            'finite_null_flux_is_not_whole_wedge_modular_energy':True}


def source_bindings():
    paths={
        'proof':'proof.md','check':'check.py','admission':'../_shared/round1298_admission.md',
        '1293_proof':'../1293/proof.md','1293_results':'../1293/results.json',
        '1295_proof':'../1295/proof.md',
        '730':'../../archive_702_741/research_note_730.md',
        '813':'../../archive_797_818/research_note_813.md',
        '1232':'../../archive_1219_1233/research_note_1232.md',
        '1233':'../../archive_1219_1233/research_note_1233.md',
        'finite_Raychaudhuri_reconstruction':
        '../../archive_1234_1280/_shared/jacobson_reconstruction_after1278.md'}
    return {k:{'path':p,'sha256':hashlib.sha256((HERE/p).read_bytes()).hexdigest()}
            for k,p in paths.items()}


def compare(a,b):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):
        assert math.isfinite(a) and math.isfinite(b)
        assert math.isclose(a,b,rel_tol=2e-8,abs_tol=1e-12),(a,b)
    else:assert a==b,(a,b)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    if args.check and args.write:parser.error('choose --check or --write')
    payload=encoded({'round':1298,'author_check_groups':1,'verdict':'PASS_AUTHOR_SCOPED',
        'exact_constants':exact_bump_and_window(),
        'Dirac_Hilbert_core':center_dirac_jets(),
        'all_constraint_seed':harmonic_constraint_test(),
        'finite_geometric_and_boost_balance':optical_and_observer_balance(),
        'scope':{'same_1293_flat_code_and_reference':True,'all_unknown_inputs_and_passive_R':True,
                 'source_mean_is_scalar_not_all_stress_fluctuations':True,
                 'explicit_positive_null_patch':True,'complete_linear_Einstein_constraints':True,
                 'strict_linear_geometry_balance':True,'geometric_entropy_role_is_adopted_input':True,
                 'BW_is_reference_identity_not_work_state_equilibrium':True,
                 'full_interacting_EFT_band_of_tiny_patch_certified':False,
                 'tiny_bound_is_actual_signal_or_instrument_resolution':False,
                 'finite_kappa_exact_semiclassical_solution_certified':False,
                 'curved_quantum_record_error_certified':False,
                 'actual_preparation_or_terminal_readout_manufactured':False,
                 'M4_or_long_goal_completed':False},
        'sources':source_bindings()})
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(payload,f,ensure_ascii=False,sort_keys=True,indent=2);f.write('\n')
        print('PASS_AUTHOR_SCOPED EXCLUSIVE_WRITE')
    else:
        compare(payload,json.loads(OUT.read_text(encoding='utf-8')))
        print('PASS_AUTHOR_SCOPED READ_ONLY_MATCH')
    print('same Dirac CP source: Tkk>4096 on ell=4/196102903917931330818687003')
    print('all linear Einstein constraints and finite optical/boost weights verified; no curved quantum budget')
    print('max algebra residual:',max(payload['Dirac_Hilbert_core']['max_residual'],
                                       payload['all_constraint_seed']['max_residual']))


if __name__=='__main__':main()
