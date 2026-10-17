"""Independent round 1260 checks; imports no author or historical checker.

Two groups: exact rational global estimates; an independently assembled
36-dimensional invariant material corner with an explicit embedding into
all 864 physical coordinates. Continuous propagation claims belong to proof.
Default only compares saved results. --write-results uses exclusive creation.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / 'independent_results.json'
ROOT = HERE.parents[1]


def exact(x):
    return {'fraction': str(x), 'display': float(x)}


def group_rational():
    mass_lo, mass_hi = Q(26), Q(32)
    source, radius, core, sigma = map(Q, (10**24, 10**30, 10**28, 10**9))
    phi, velocity, time_cap = Q(1, 10**6), Q(1, 100), Q(6*10**11)
    umax, slope = Q(3, 2)/core, 1/core**2
    p_axis, p2_all = 1/(2*sigma), 1/sigma**2
    n = 80_000_000_000
    # Independent exact moment bounds for the six real Gaussian coordinates.
    assert 3 < 4 and 15 < 16 and 6 < 9 and 135135 < 384**2
    bounds = {
        'free_SR': time_cap*p2_all*(1/(2*mass_lo)+1/(2*source)),
        'newton_lipschitz': time_cap*mass_hi*source*slope*3*sigma,
        'field_G2_lipschitz': time_cap*mass_hi*source*(mass_hi+source)*umax*slope*3*sigma,
        'field_G2_finite_M_squared': time_cap*phi**2*mass_hi**2/(2*source),
    }
    # The summed divergence self terms use u p^2 - i grad(u).p.
    # Their three weak-derivative terms use a deliberately loose p<=1/sigma.
    self_residual = umax*p2_all + sum(slope/sigma for _ in range(3))
    bounds['velocity_self'] = time_cap*Q(3, 2)*(source/mass_lo+mass_hi/source)*self_residual
    # Cross terms are built one coordinate/index at a time, not from authors' totals.
    scalar_one = umax*p_axis*p_axis + slope*p_axis
    radial_one = (1/core)*p_axis*p_axis + 5*slope*p_axis
    bounds['velocity_scalar'] = time_cap*sum(Q(7, 2)*scalar_one for _ in range(3))
    bounds['velocity_radial'] = time_cap*sum(Q(1, 2)*radial_one for _ in range(9))
    total = sum(bounds.values())
    assert total < Q(1211975, 10**14) < Q(1, 10**7)
    # 3 + 7/2 + 1/2 comes from vector norms, not nine independent worst cases.
    velocity_prefactor = Q(3) + Q(7, 2) + Q(1, 2)
    positive = 1-source*umax*(1+velocity_prefactor*velocity**2)
    assert positive == Q(199969979, 200000000) and positive > Q(999, 1000)
    pi_lo, pi_hi = Q(25, 8), Q(22, 7)
    S_hi = (2*n+Q(1, 2))*pi_hi
    assert S_hi/(1-phi) < time_cap
    delta_lo = (2*n+Q(1, 2))*pi_lo*phi**2/(2*(1-phi))
    assert delta_lo > Q(1, 4)
    assert (Q(n)+Q(1, 4))*phi**2 < (1-phi)/12
    derivative_lo = 3*Q(9, 10)**5*Q(1, 4)
    assert derivative_lo == Q(177147, 400000) and derivative_lo > Q(3, 8)
    coarse_each = Q(1, 10**7)
    gap = Q(3, 32)-2*coarse_each
    assert gap > Q(9, 100)
    same_age_bound = 2*coarse_each
    conditional_gap = gap-2*Q(1, 1000)
    assert conditional_gap > Q(9, 100)
    sr_remainder = time_cap*Q(6, 16)/sigma**6*(mass_lo**-5+source**-5)
    assert sr_remainder < Q(19, 10**51)
    known_static_2pn = time_cap*mass_hi*source*(mass_hi**2+5*mass_hi*source+source**2)/(4*radius**3)
    assert known_static_2pn < Q(4801, 10**9)
    return {
        'bounds': {k: exact(v) for k,v in bounds.items()},
        'total_isometry_bound': exact(total),
        'positive_mass_coefficient': exact(positive),
        'phase_lower': exact(delta_lo),
        'derivative_lower': exact(derivative_lo),
        'fixed_coordinate_probability_gap_lower': exact(gap),
        'same_internal_age_probability_gap_upper': exact(same_age_bound),
        'conditional_extra_001_each_gap_lower': exact(conditional_gap),
        'SR_beyond_1PN_only': exact(sr_remainder),
        'one_known_static_2PN_piece_only': exact(known_static_2pn),
        'continuous_claims_not_proved_by_finite_arithmetic': True,
    }


def group_material():
    # Reconstruct the old clock/gates directly, with no imported checker.
    x = np.array([[0,1],[1,0]], complex)
    y = np.array([[0,-1j],[1j,0]], complex)
    i2 = np.eye(2)
    ur = np.kron((i2-y)/2, i2)-1j*np.kron((i2+y)/2, y)
    um = np.kron(np.diag([1,0]), i2)-1j*np.kron(np.diag([0,1]), x)
    r_aux = np.kron(ur, i2)
    m_aux = np.kron(i2, um)
    aux0 = np.kron(np.array([1,1j])/math.sqrt(2), np.kron([1,0],[1,0]))
    aux = [aux0, aux0, r_aux@aux0, m_aux@r_aux@aux0]
    clock = 1.5*np.eye(4)
    for j in range(3):
        clock[j,j+1] = clock[j+1,j] = math.sqrt((j+1)*(3-j))/2
    swap = np.zeros((9,9))
    for q in range(3):
        for r in range(3): swap[3*r+q,3*q+r] = 1
    # Columns: head, logical qutrit, true external R qutrit.
    embedding = np.zeros((864,36), complex)
    view = embedding.reshape(4,8,3,3,3,36)
    for h in range(4):
        for q in range(3):
            for r in range(3):
                col = 9*h+3*q+r
                qc, qr = (q,0) if h == 0 else (0,q)
                view[h,:,qc,qr,r,col] = aux[h]
    assert np.linalg.norm(embedding.conj().T@embedding-np.eye(36)) < 1e-13
    physical_flag = np.broadcast_to((np.arange(8)%2)[None,:,None,None,None],(4,8,3,3,3)).reshape(864)
    compressed_flag = embedding.conj().T@(physical_flag[:,None]*embedding)
    expected_flag = np.kron(np.diag([0,0,0,1]),np.eye(9))
    flag_error = float(np.linalg.norm(compressed_flag-expected_flag))
    assert flag_error < 1e-13
    # Direct physical nearest-head blocks confirm this is the original FK corner.
    gate1 = np.kron(np.eye(8),swap)
    gate2 = np.kron(ur,np.eye(18))
    gate3 = np.kron(np.kron(i2,um),np.eye(9))
    mstar = 17+2*math.sqrt(19)
    blocks = embedding.reshape(4,72,-1)
    applied = (mstar+1.5)*blocks.copy()
    for j,g in enumerate([gate1,gate2,gate3]):
        c = math.sqrt((j+1)*(3-j))/2
        applied[j+1] += c*(g@blocks[j])
        applied[j] += c*(g.conj().T@blocks[j+1])
    intertwiner_error = float(np.linalg.norm(applied.reshape(864,36)-embedding@np.kron(mstar*np.eye(4)+clock,np.eye(9))))
    assert intertwiner_error < 1e-12
    eig, vec = np.linalg.eigh(clock)
    assert np.max(abs(eig-np.arange(4))) < 1e-14
    delta = math.pi/4
    age = vec@(np.exp(-1j*delta*eig)*(vec.conj().T@np.eye(4)[:,0]))
    load = np.zeros((9,2),complex)
    load[3,0] = load[6,1] = 1
    prepared = np.kron(age[:,None],load)
    projectors = [np.zeros((36,36),complex) for _ in range(6)]
    for j in range(4):
        head_projector = np.outer(vec[:,j],vec[:,j].conj())
        for sign in (1,-1):
            projectors[j+1-sign] += np.kron(head_projector,(np.eye(9)+sign*swap)/2)
    assert np.linalg.norm(sum(projectors)-np.eye(36)) < 1e-13
    weights = []
    for k, mult in enumerate((1,3,4,4,3,1)):
        pullback = prepared.conj().T@projectors[k]@prepared
        assert np.linalg.norm(pullback-mult*np.eye(2)/16) < 1e-13
        weights.append(float(np.trace(pullback).real/2))
    interaction = np.kron(np.eye(4),np.eye(9)-swap)
    assert np.linalg.norm(prepared.conj().T@interaction@prepared-np.eye(2)) < 1e-13
    assert 26 < mstar+1 and mstar+6 < 32
    phi = 1e-6
    source = 1e24
    sage = 2*math.pi*80_000_000_000+math.pi/2
    tminus = sage/(1-phi)
    tplus = sage/(1-phi+phi**2/2)
    phase_shift = tminus*phi**2/2
    base_mass = mstar+1
    def evolve(theta, quadratic):
        # Remove only ONE global scalar phase, never a branch-dependent phase.
        terms = [np.exp(-1j*(theta*k+quadratic*((base_mass+k)**2-base_mass**2)))*projectors[k] for k in range(6)]
        return sum(terms)@prepared
    records = []
    maps = []
    for label,theta,quad in (
        ('minus_fixed_T',math.pi/2,0),
        ('plus_fixed_T',math.pi/2+phase_shift,tminus*phi**2/(2*source)),
        ('plus_same_internal_age',math.pi/2,tplus*phi**2/(2*source))):
        out = evolve(theta,quad)
        effect = out.conj().T@compressed_flag@out
        p = float(np.trace(effect).real/2)
        err = float(np.linalg.norm(effect-p*np.eye(2),2))
        assert err < 1e-13
        assert np.linalg.norm(out.conj().T@out-np.eye(2)) < 1e-13
        maps.append(out)
        records.append({'model':label,'Mr_probability':p,'all_input_effect_scalar_error':err})
    gap = records[1]['Mr_probability']-records[0]['Mr_probability']
    assert gap > 0.18
    assert abs(records[2]['Mr_probability']-records[0]['Mr_probability']) < 1e-13
    # A concrete entangled input is diagnostic only; full-reference quantifier
    # follows analytically from the preparation-isometry estimate.
    entangled = np.array([1,0,0,0,0,1j],complex)/math.sqrt(2)
    ref_identity = np.eye(3)
    entanglement_errors = []
    for out,record in zip(maps,records):
        full = np.kron(out,ref_identity)@entangled
        p = float(np.vdot(full,np.kron(compressed_flag,ref_identity)@full).real)
        entanglement_errors.append(abs(p-record['Mr_probability']))
    assert max(entanglement_errors) < 1e-13
    return {
        'independent_method':'explicit 864-by-36 physical embedding; 4-level clock and SWAP projectors; no checker import',
        'physical_Mr_compression_error':flag_error,
        'original_FK_intertwiner_error':intertwiner_error,
        'complete_six_mass_weights':weights,
        'interaction_source_expectation':1,
        'records':records,
        'fixed_T_reference_gap':gap,
        'entangled_reference_diagnostic_errors':entanglement_errors,
        'arbitrary_reference_proof_is_isometry_bound_not_sample':True,
    }


def calculate():
    return {
        'status':'PASS_INDEPENDENT_WITH_DECLARED_SCOPE',
        'round':1260,'independent_groups':2,
        'independent_code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'old_definition_sha256':hashlib.sha256((ROOT/'archive_1106_1119/1112/check.py').read_bytes()).hexdigest(),
        'group1_rational_full_motion_budget':group_rational(),
        'group2_original_material_actual_Mr':group_material(),
        'scope':{
            'normal_mechanical_preparation_fixed_not_all_orbital_inputs':True,
            'all_original_internal_inputs_and_passive_references':True,
            'no_Gaussian_tail_removed':True,
            'positive_completion_is_public_effective_input':True,
            'same_parent_full_source_geometry_Ward_not_certified':True,
            'actual_common_time_trigger_not_constructed':True,
            'full_GR_remainder_not_certified':True,
            'no_new_axiom_or_CT_adoption':True,
        },
    }


def compare(saved,now,path='root'):
    if isinstance(now,dict):
        assert isinstance(saved,dict) and saved.keys()==now.keys(), path
        for k in now: compare(saved[k],now[k],path+'.'+k)
    elif isinstance(now,list):
        assert isinstance(saved,list) and len(saved)==len(now),path
        for i,(a,b) in enumerate(zip(saved,now)): compare(a,b,f'{path}[{i}]')
    elif isinstance(now,float):
        assert math.isclose(float(saved),now,rel_tol=1e-11,abs_tol=1e-12),path
    else: assert saved==now,path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--write-results',action='store_true')
    args=p.parse_args()
    result=calculate()
    if args.write_results:
        with OUT.open('x',encoding='utf-8') as f: json.dump(result,f,ensure_ascii=False,indent=2); f.write('\n')
        print('WRITE_ONCE_PASS: independent_groups=2')
    else:
        compare(json.loads(OUT.read_text(encoding='utf-8')),result)
        print('READ_ONLY_MATCH: independent_groups=2')
    print('Scope: full normal-state comparison proof plus finite rational/material checks; no actual timing or full GR certificate.')

if __name__=='__main__': main()

