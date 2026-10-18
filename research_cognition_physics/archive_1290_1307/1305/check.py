"""1305 finite diagnostics; continuous-field statements remain in proof.md.
Default/--check: read only. --write: exclusive first creation of results.json.
No numerical s/kappa or new phase window is presented as an analytic certificate.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results.json'
BASE = HERE.parent


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def seed_rational():
    r0 = F(1, 32)
    K = 12 * (1 + 344 / r0**2 + 2570160 / r0**4)
    assert K == 32340101333004
    bound = 48 * K**2
    ell = 1024 / bound
    assert bound == 50202343402990420689583872768
    derivative_coeff = F(1, 2)*8 + F(3, 2)*16
    assert derivative_coeff == 28 < 48
    center = F(3, 2)*8192
    margin = center - 4*ell*bound
    assert center == 12288 and margin == 8192
    J_lower = margin * (4*ell**2) * ell**2 / 2
    b0 = 8192 * ell**4
    assert J_lower == 2*b0 > b0 > 0
    assert 9*(512/r0**2+2) < 2173**2
    return {'K': str(K), 'derivative_bound_M': str(bound), 'ell': str(ell),
            'single_rank_derivative_coefficient': 8,
            'two_cross_rank_derivative_coefficient': 16,
            'three_mode_derivative_coefficient': str(derivative_coeff),
            'strict_center_lower': str(center),
            'strict_cap_lower': str(margin),
            'strict_J_lower': str(J_lower), 'conservative_b0': str(b0),
            'relative_field_energy_strict_upper': 2173,
            'geometry_is_rebuilt_from_three_mode_source': True}


def annihilation(n, j):
    d = 1 << n
    a = np.zeros((d, d), dtype=complex)
    bit = 1 << (n-1-j)
    for state in range(d):
        if state & bit:
            before = sum((state >> (n-1-k)) & 1 for k in range(j))
            a[state ^ bit, state] = (-1)**before
    return a


def car_partner_diagnostic():
    rng = np.random.default_rng(1305)
    aa = [annihilation(6, j) for j in range(6)]
    ident = np.eye(64)
    car_error = 0.
    for i in range(6):
        for j in range(6):
            target = ident if i == j else np.zeros_like(ident)
            car_error = max(car_error, float(np.max(np.abs(
                aa[i]@aa[j].conj().T + aa[j].conj().T@aa[i] - target))))
    assert car_error < 1e-14
    # Field mode order: encoded F_1,F_2,F_3 then complement Q_1,Q_2,Q_3.
    P = .5*np.block([[np.eye(3), -np.eye(3)], [-np.eye(3), np.eye(3)]])
    V = np.vstack([np.eye(3), np.zeros((3,3))])
    FF = V@V.conj().T
    fock_positions = [1 << (2-a) for a in range(3)]
    max_trace_error = max_source_error = max_partner_error = 0.
    sample_count = 0
    for index in range(3):
        vec = rng.normal(size=9) + 1j*rng.normal(size=9)
        vec /= np.linalg.norm(vec)
        rho = .73*np.outer(vec, vec.conj()) + .27*np.eye(9)/9
        if index == 0:
            vec = np.eye(3).reshape(9)/np.sqrt(3)
            rho = np.outer(vec, vec.conj())
        rho_Q = np.trace(rho.reshape(3,3,3,3), axis1=1, axis2=3)
        rho_R = np.trace(rho.reshape(3,3,3,3), axis1=0, axis2=2)
        delta = -FF@P-P@FF+FF@P@FF+V@rho_Q@V.conj().T
        traced = np.trace(delta.reshape(2,3,2,3), axis1=1, axis2=3)
        expected = np.array([[-.5,1.5],[1.5,0.]])
        max_trace_error = max(max_trace_error, float(np.max(np.abs(traced-expected))))
        # Full CP-prepared Fock state, not a Gaussian reconstruction of rho.
        omega = np.zeros((192,192), dtype=complex)
        for q in range(8):
            embed = np.zeros((192,9), dtype=complex)
            for a in range(3):
                for r in range(3):
                    embed[(fock_positions[a]*8+q)*3+r,a*3+r] = 1.
            omega += embed@rho@embed.conj().T/8
        red_R = np.trace(omega.reshape(64,3,64,3), axis1=0, axis2=2)
        max_partner_error = max(max_partner_error, float(np.max(np.abs(red_R-rho_R))))
        assert abs(np.trace(omega)-1) < 2e-14
        for _ in range(3):
            z = rng.normal(size=(2,2)) + 1j*rng.normal(size=(2,2))
            source2 = (z+z.conj().T)/2
            source = np.kron(source2, np.eye(3))
            zR = rng.normal(size=(3,3)) + 1j*rng.normal(size=(3,3))
            partner = (zR+zR.conj().T)/2
            partner /= max(1., np.linalg.norm(partner, 2))
            dGamma = np.zeros((64,64), dtype=complex)
            for i in range(6):
                for j in range(6):
                    if abs(source[i,j]) > 0:
                        dGamma += source[i,j]*aa[i].conj().T@aa[j]
            observed = np.trace(omega@np.kron(dGamma,partner))
            reference = np.trace(P@source)*np.trace(rho_R@partner)
            predicted = np.trace(expected@source2)*np.trace(rho_R@partner)
            max_source_error = max(max_source_error, float(abs(observed-reference-predicted)))
            sample_count += 1
    assert max_trace_error < 2e-12 and max_source_error < 2e-12 and max_partner_error < 2e-12
    return {'fermion_modes': 6, 'field_fock_dimension': 64, 'partner_dimension': 3,
            'unknown_partner_test_states': 3, 'source_partner_tests': sample_count,
            'CAR_max_error': round(car_error, 14),
            'three_mode_trace_max_error': round(max_trace_error,14),
            'source_partner_max_error': round(max_source_error,14),
            'partner_reduction_max_error': round(max_partner_error,14),
            'source_trace_coefficients': ['-1/2', '3/2', '3/2'],
            'finite_fixture_not_continuum_or_all_R_proof': True}


def hermitian_basis(n):
    out=[]
    for i in range(n):
        b=np.zeros((n,n),complex); b[i,i]=1.; out.append(b)
    for i in range(n):
        for j in range(i+1,n):
            b=np.zeros((n,n),complex); b[i,j]=b[j,i]=1/np.sqrt(2); out.append(b)
            b=np.zeros((n,n),complex); b[i,j]=-1j/np.sqrt(2); b[j,i]=1j/np.sqrt(2); out.append(b)
    return out


def vandermonde(xs):
    v=F(1)
    for i in range(len(xs)):
        for j in range(i+1,len(xs)):
            v *= xs[j]-xs[i]
    return v


def path_diagnostics():
    mu=F(1,4096)
    U=np.exp(2j*np.pi*np.outer(np.arange(3),np.arange(3))/3)/np.sqrt(3)
    records=[]
    for s in [F(0),F(1,16),F(1,4),F(1,2),F(1)]:
        masses=[1+s*(j*mu-1) for j in [1,2,3]]
        charged=[2+s*j for j in [0,1,2]]
        nu=[x*x for x in masses]; ee=[x*x for x in charged]
        vnu,vch=vandermonde(nu),vandermonde(ee)
        off2=(sum(x*x for x in nu)-sum(nu[i]*nu[j] for i in range(3) for j in range(i+1,3)))/9
        invariant_squared=F(1,3)*vnu*vnu*vch*vch
        assert all(mu <= x <= 1 for x in masses)
        if s:
            assert off2>0 and invariant_squared>0 and len(set(masses))==3 and len(set(charged))==3
        else:
            assert off2==0 and invariant_squared==0
        Hnu=U@np.diag([float(x) for x in nu])@U.conj().T
        He=np.diag([float(x) for x in ee])
        cols=[]
        for b in hermitian_basis(3):
            z=np.r_[ (b@Hnu-Hnu@b).ravel(), (b@He-He@b).ravel() ]
            cols.append(np.r_[z.real,z.imag])
        singular=np.linalg.svd(np.column_stack(cols),compute_uv=False)
        rank=int(np.sum(singular > 1e-12)) if s else 0
        if s: assert rank==8, (s,singular)
        a,b=1-s,s*mu
        real_off=[]
        for i,j in [(1,2),(1,3),(2,3)]:
            value=i*j*(2*a*b*(36-7*(i+j))/49+b*b*(14-i*i-j*j)/7)
            if s: assert value>0
            real_off.append(str(value))
        records.append({'s_diagnostic_not_joint_window':str(s), 'neutrino_masses':list(map(str,masses)),
                        'charged_masses':list(map(str,charged)), 'DFT_offdiagonal_modulus_squared':str(off2),
                        'weak_basis_trace_imaginary_invariant_squared':str(invariant_squared),
                        'hermitian_commutant_dimension':9-rank,
                        'real_comparison_offdiagonals':real_off})
    assert [1+1*(j*mu-1) for j in [1,2,3]] == [mu,2*mu,3*mu]
    return {'J_DFT_squared':'1/108','samples':records,
            'simple_irreducible_for_all_positive_s_by_proof':True,
            'sample_s_is_not_a_certified_thermal_geometry_window':True}


def inherited_endpoint():
    p=BASE/'1303'/'check.py'
    spec=importlib.util.spec_from_file_location('round1303_readonly',p)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    old=module.rational_certificate()
    assert F(old['cycle_strict_lower_bound'])>0
    return {'source':'../1303/check.py:rational_certificate',
            'cycle_at_s_1_strict_lower':old['cycle_strict_lower_bound'],
            'all_t_in_original_window':True,
            'fixed_positive_normalized_time_weight_preserves_bound':True,
            'analytic_nonzero_near_s_0_is_proof_not_samples':True,
            'near_zero_sign_inferred_from_endpoint':False}


def split_budgets():
    base=F(1,64)+F(5,2048)+3*F(1,1024)
    perturb=F(1,1024)
    assert base==F(43,2048) and base+perturb==F(45,2048)<F(1,32)
    assert 2*(base+perturb)==F(90,2048)<F(1,8)
    assert 2*(base+perturb)+F(1,32)==F(154,2048)<F(1,8)
    direct_noise=F(5,2048)+3*F(1,1024)
    assert direct_noise+perturb==F(13,2048)<F(45,2048)
    transport=F(1,8); calibration=6*F(1,96); statistics=6*F(1,96)
    assert transport+calibration+statistics==F(1,4)
    assert 2*96**2//4==4608
    return {'units':'thermal entries in m; phase entries in d=abs(Cbar(s)); not equated',
            'coarse_total_bound_in_m':str(base+perturb),
            'direct_fixed_kernel_bound_in_m':str(direct_noise+perturb),
            'coarse_pair_bound_in_m':'90/2048', 'coarse_pair_with_base_thermal_bias_in_m':'154/2048',
            'phase_transport_in_d':str(transport), 'phase_calibration_sum_in_d':str(calibration),
            'phase_statistical_sum_in_d':str(statistics), 'phase_total_in_d':'1/4',
            'Hoeffding_all_18_probabilities':'N >= 4608/d^2 * log(36/eta)',
            'all_budgets_positive_but_permissions_are_inputs':True,
            'actual_selected_s':None,'actual_selected_kappa':None,'actual_calibrated_C':None,
            'source_geometry_residual':'O(kappa^2)+kappa*O(s)',
            'entropy_source_geometry_balance':'O(kappa)+O(s)',
            'actual_temperature_extra_error':'omega_beta(s) -> 0, not asserted O(s)',
            'fixed_Lambda_baseline':True,'unknown_input_dependent_kernel':False}


def source_hashes():
    paths=[HERE/'proof.md',HERE/'check.py',HERE/'NEXT.md',BASE/'research_note_1305.md',
           BASE/'_shared'/'round1305_admission.md',BASE/'1298'/'proof.md',
           BASE/'1300'/'proof.md',BASE/'1300'/'acceptance.json',
           BASE/'1303'/'proof.md',BASE/'1303'/'check.py',BASE/'1303'/'results.json',
           BASE/'_shared'/'nondegenerate_joint_corollary_after1300.md',
           BASE.parent/'archive_1009_1043'/'research_note_1011.md']
    return {str(p.relative_to(BASE.parent)).replace('\\','/'):digest(p) for p in paths}


def compute():
    return {'round':1305,'status':'AUTHOR_FINITE_DIAGNOSTICS_PASS',
            'scope':'same Dirac quadratic background and declared finite effective permissions',
            'seed':seed_rational(),'CAR_partner':car_partner_diagnostic(),
            'path':path_diagnostics(),'endpoint':inherited_endpoint(),
            'budgets':split_budgets(),'frozen_source_hashes':source_hashes(),
            'claims_not_made':['full CP/T experiment','three generations unique','complete SM',
                               'dynamic Gauss realization','actual nonzero four-fermion exchange',
                               'exact nonlinear semiclassical solution','M4 or goal completed']}


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    data=compute()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2,sort_keys=True)
            stream.write('\n')
        mode='exclusive first write'
    elif OUT.exists():
        stored=json.loads(OUT.read_text(encoding='utf-8'))
        assert data==stored,'Saved results or source bindings differ.'
        mode='read-only recomputation and comparison'
    elif args.check:
        raise FileNotFoundError('results.json not yet created')
    else:
        mode='read-only diagnostics; no saved results yet'
    print(json.dumps({'status':'PASS','mode':mode,'round':1305,
                      'source_bindings':len(data['frozen_source_hashes']),
                      'CAR_partner_tests':data['CAR_partner']['source_partner_tests'],
                      'positive_seed_and_analytic_endpoint':True,
                      'invented_s_or_kappa':False},ensure_ascii=False))

if __name__=='__main__':
    main()
