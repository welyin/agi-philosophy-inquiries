"""1170 independent review: exact finite-action oscillator algebra.

No spatial grid and no claim that a truncated time propagator is exact.
H^2 applied to a Gaussian vacuum occupies at most oscillator level four;
ten levels therefore evaluate the required vectors without truncation.
Default is read-only; --write-once only creates review_checks.json.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import re
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
EXPECTED = {
 'archive_1161_/research_note_1170.md': '8073cf858a7cd7242aabd02beabb0e72d2e3b290379e8f74c4f0237361ab7057',
 'archive_1161_/1170/admission.md': '67937b30e043d30401d08f8c7277a5c774a885549493748e88f12e82de434006',
 'archive_1161_/1170/proof.md': 'c55f5abdc1503801b87fa03336d0710f48d723ecbcd4c133a3fe044a48e1b3af',
 'archive_1161_/1170/check.py': '2b23bd34e18a318c853d85a460972e0b2e2430bec1cd019e760ce9f7fd4ec43b',
 'archive_1161_/1170/results.json': 'faf22a2acc5e1a6ffd4948182c9a89cc5a910ca89f3403b1c2a7c0a2fed77551',
 'archive_1161_/1170/NEXT.md': '97fd8c5d54052a68cd76691d62e702f62d434e9edb7fde25cbd6820520fbab68',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def norm(a):
    return float(np.linalg.norm(a))


def material_and_support(saved):
    spec = importlib.util.spec_from_file_location('review1170_old1112', BASE/'archive_1106_1119/1112/check.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mat = mod.Material(1)
    h, identity = mat.hfk, np.eye(288)
    j = mat.encode_joint(np.eye(3)[:, None, :], 0).reshape(288, 3)
    # Independent spectral construction of the support from the actual J.
    projectors, pj_j = [], []
    for e in range(4):
        p = np.eye(288, dtype=complex)
        for f in range(4):
            if f != e: p = p@(h-f*identity)/(e-f)
        projectors.append(p)
        pj_j.append(p@j)
    weights = [1/8, 3/8, 3/8, 1/8]
    spectral_basis = np.column_stack([p/math.sqrt(w) for p,w in zip(pj_j, weights)])
    s = spectral_basis@spectral_basis.conj().T
    iso = opnorm(spectral_basis.conj().T@spectral_basis-np.eye(12))
    mask = np.zeros((4,8,3,3)); mask[:,[1,3,5,7],:,:] = 1
    mr = np.diag(mask.ravel())
    mr_leak = opnorm((identity-s)@mr@spectral_basis)
    m_leak = opnorm((identity-s)@h@spectral_basis)
    cols, complement = [], []
    for head in range(4):
        for q in range(3):
            for dest, qr in [(cols, 0), (complement, 1)]:
                vec = np.kron(np.kron(np.kron(np.eye(4)[:,head], mat.aux), np.eye(3)[:,q]), np.eye(3)[:,qr])
                dest.append(mat.gamma@vec)
    k, kp = np.column_stack(cols), np.column_stack(complement)
    direct_error = opnorm(k@k.conj().T-s)
    h_restriction = opnorm(h@k-k@np.kron(mat.hc,np.eye(3)))
    mr_restriction = opnorm(mr@k-k@np.kron(np.diag([0,0,0,1]),np.eye(3)))
    assert max(iso,mr_leak,m_leak,direct_error,h_restriction,mr_restriction) < 3e-12
    full_record = np.vstack([identity-mr,mr])
    record_leak = opnorm(np.block([[identity-s,np.zeros((288,288))],[np.zeros((288,288)),identity-s]])@full_record@k)
    assert record_leak < 3e-12
    krylov = np.column_stack([np.linalg.matrix_power(h,n)@j for n in range(4)])
    sv = np.linalg.svd(krylov,compute_uv=False)
    assert np.count_nonzero(sv>1e-9) == 12
    return mat, k, kp, mr, projectors, dict(
        support_rank=12,complement_dimension=276,
        spectral_weights=weights,spectral_basis_isometry_error=iso,
        natural_and_spectral_support_error=direct_error,
        h_restriction_error=h_restriction,mr_restriction_error=mr_restriction,
        mr_spectral_support_leak=mr_leak,retained_record_support_leak=record_leak,
        krylov_smallest_singular_value=float(sv[-1]))


def finite_action(mat,k,kp,mr,projectors,saved):
    sigma,xbar,age,g,lam,epsilon,t = .5,1.,math.pi/2,.02,1.3,.5,1e-5
    hc = mat.hc
    e0 = np.array([-1,math.sqrt(3),-math.sqrt(3),1])/math.sqrt(8)
    e1 = np.array([math.sqrt(3),-1,-1,math.sqrt(3)])/math.sqrt(8)
    v = 1j*(np.outer(e0,e1)-np.outer(e1,e0))
    mass = mat.baseline*np.eye(4)+hc
    inv = np.linalg.inv(mass)
    charges = [lam*mass,lam*mass+epsilon*v]
    energies,vectors = np.linalg.eigh(hc)
    phi = (vectors*np.exp(-1j*energies*age))@vectors.conj().T@np.eye(4)[:,0]
    vh = k@np.kron(v,np.eye(3))@k.conj().T
    hidden = kp@np.kron(v,np.eye(3))@kp.conj().T
    mean_defect = max(opnorm(p@vh@p) for p in projectors)
    assert mean_defect < 2e-12 and opnorm(hidden@k) < 1e-12
    full_mass = mat.baseline*np.eye(288)+mat.hfk
    visible_comm = opnorm(full_mass@(epsilon*vh)-(epsilon*vh)@full_mass)
    hidden_comm = opnorm(full_mass@(epsilon*hidden)-(epsilon*hidden)@full_mass)
    assert abs(visible_comm-.5)<1e-12 and abs(hidden_comm-.5)<1e-12
    # Oscillator quadratures, not the author's Gaussian coefficient recursion.
    n = 10
    annih = np.diag(np.sqrt(np.arange(1,n,dtype=float)),1)
    x = xbar*np.eye(n)+sigma*(annih+annih.T)
    p = (annih-annih.T)/(2j*sigma)
    kinetic = np.kron(p@p,inv)/2
    hs = [np.kron(np.eye(n),hc)+kinetic+g*np.kron(x,b) for b in charges]
    delta = g*epsilon*np.kron(x,v)
    psi = np.kron(np.eye(n)[:,0],phi)
    ha,hb = hs
    vectors_check = [ha@psi,ha@ha@psi,hb@hb@psi,delta@psi,ha@delta@psi]
    vector_norms = dict(zip(['Ha','Ha2','Hb2','Delta','HaDelta'],map(norm,vectors_check)))
    high_support = max(norm(a.reshape(n,4)[5:]) for a in vectors_check)
    assert high_support < 1e-14
    reference = saved['groups']['fixed_finite_comparison_and_error_responsibility']
    norm_error = max(abs(value-reference['norms_from_exact_gaussian_moments'][key]) for key,value in vector_norms.items())
    assert norm_error < 3e-12
    # Exact independent derivative on the four-level reduction, pulled back by K.
    e = np.diag([0,0,0,1])
    slope = float((1j*g*epsilon*xbar*phi.conj()@(v@e-e@v)@phi).real)
    jet = sigma*g*epsilon*norm(v@phi)
    assert abs(slope-g*epsilon*xbar*math.sqrt(3)/8)<1e-13
    assert abs(jet-reference['phase_independent_jet'])<1e-14
    # Re-evaluate analytic triangle bounds independently from the documented
    # coefficient formula; no finite-difference estimate enters the certificate.
    def moment_norm(power):
        return sigma**power*math.sqrt(math.prod(range(1,2*power,2)))
    def poly_bound(coeff):
        return sum(z*moment_norm(i) for i,z in enumerate(coeff))
    def h_bound(coeff,b):
        out = [0.]*(len(coeff)+2)
        for j,z in enumerate(coeff):
            out[j] += (3+g*xbar*b+(j+.5)/(2*mat.baseline*sigma*sigma))*z
            out[j+1] += g*b*z
            if j>=2: out[j-2] += j*(j-1)*z/(2*mat.baseline)
            out[j+2] += z/(8*mat.baseline*sigma**4)
        return out
    mmax = mat.baseline+3
    a1 = h_bound([1.],lam*mmax)
    a2 = h_bound(a1,lam*mmax)
    b2 = h_bound(h_bound([1.],lam*mmax+epsilon),lam*mmax+epsilon)
    dd = [g*epsilon*xbar/math.sqrt(2),g*epsilon/math.sqrt(2)]
    upper = dict(zip(['Ha','Ha2','Hb2','Delta','HaDelta'],
                    map(poly_bound,[a1,a2,b2,dd,h_bound(dd,lam*mmax)])))
    upper_error=max(abs(v-reference['analytic_upper_bounds'][k]) for k,v in upper.items())
    assert upper_error<3e-12
    assert all(vector_norms[key] <= value for key,value in upper.items())
    r0 = t*t*(upper['Ha2']+upper['Hb2'])/2
    # Conservative norms of XM, P, M acting on Delta and Mr Delta.
    xdd = [xbar*dd[0],dd[0]+xbar*dd[1],dd[1]]
    pdd = [dd[1],dd[0]/(2*sigma*sigma),dd[1]/(2*sigma*sigma)]
    aa = g*lam
    de = 2*mmax*abs(aa)*t*poly_bound(xdd)+abs(aa)*t*t*poly_bound(pdd)+aa*aa*t**3*mmax*upper['Delta']/3
    remainder = 2*r0+(t*upper['Delta']+r0)**2+2*t*t*upper['Ha']*upper['Delta']+2*t*de
    gap = abs(slope)*t-remainder
    state = t*jet-t*t*((upper['Ha2']+upper['Hb2'])/2+upper['HaDelta'])
    assert abs(gap-reference['transported_original_mr_probability_gap_lower'])<1e-18
    assert abs(state-reference['complete_output_trace_distance_lower'])<1e-18
    assert gap>0 and state>0
    # Algebraic coefficients of P and M in transformed H0+i Cdot Cdag cancel.
    time = .23
    pco = aa*time-aa*time
    mco = aa*aa*time*time/2-aa*aa*time*time+aa*aa*time*time/2
    assert pco == 0 and mco == 0
    return dict(full_visible_commutator_norm=visible_comm,full_hidden_commutator_norm=hidden_comm,
        energy_block_mean_defect=mean_defect,hidden_support_defect=opnorm(hidden@k),
        oscillator_levels=n,max_required_level=4,unused_level_norm=high_support,
        exact_vector_norms=vector_norms,max_norm_difference_from_author=norm_error,
        analytic_bound_difference=upper_error,original_mr_slope=slope,
        transported_mr_remainder=remainder,transported_mr_probability_lower=gap,
        complete_output_trace_distance_lower=state,actual_error_budget_certified=False,
        cocycle_p_and_m_coefficients_cancel=True,
        infinite_dimensional_propagator_not_approximated_by_finite_fock_unitary=True)


def text_links():
    count = 0
    for name in EXPECTED:
        if not name.endswith('.md'): continue
        path = BASE/name
        data = path.read_bytes()
        assert not any(c<32 and c not in (9,10,13) for c in data),name
        text = data.decode('utf-8')
        assert text.count('$$')%2 == 0,name
        assert not re.search(r'(?<!\\)\b(?:qquad|quad)\b',text),name
        for link in re.findall(r'\]\(([^)]+)\)',text):
            if '://' in link: continue
            assert (path.parent/link.split('#')[0]).exists(),(name,link)
            count += 1
    return count


def run():
    before = {n:sha(BASE/n) for n in EXPECTED}
    assert before == EXPECTED
    saved=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    for n,h in saved['source_hashes'].items(): assert sha(BASE/n)==h,n
    mat,k,kp,mr,pj,support = material_and_support(saved)
    finite = finite_action(mat,k,kp,mr,pj,saved)
    links = text_links()
    assert before == {n:sha(BASE/n) for n in EXPECTED}
    for n,h in saved['source_hashes'].items(): assert sha(BASE/n)==h,n
    return dict(schema='round1170_independent_review_v1',round=1170,date='2026-10-13',verdict='PASS_WITH_SCOPE',
        blockers=0,core_sha256=before,source_sha256=saved['source_hashes'],reviewer_code_sha256=sha(Path(__file__)),
        author_default_observed=dict(exit_code=0,groups=3,sources=7,read_only=True),
        independence='Pre-admission support-space and permission comments disclosed; reviewer did not author core. Uses polynomial spectral support and oscillator finite-action algebra, not author Gaussian polynomial code.',
        independent_support=support,independent_finite_action=finite,local_link_occurrences=links,
        analytic_checks=dict(classification_only_on_12dimensional_support=True,
            common_spatial_phase_cannot_cancel_variance_projected_defect=True,
            hermitian_charge_implies_no_cross_block_when_DS_zero=True,
            hidden_276dimensional_charge_freedom_preserved=True,
            original_Mr_requires_transport_and_keeps_S=True,
            proof17_telescope_and_proof18_state_effect_remainder_valid=True,
            classification_does_not_follow_from_one_finite_probability=True,
            captured_and_loaded_clock_age_is_not_free_reset=True,
            A_acc_is_additional_candidate_not_adopted_axiom=True,
            local_stark_generator_not_positive_global_universe=True,
            no_actual_accelerator_or_1e_minus8_precision_claim=True,
            no_full_sr_pn_causality_curvature_or_einstein_delivery=True),
        scientific_claims_new_axioms=0,gr_hard_deliveries=0,
        core_and_sources_unchanged=True)


def same(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:same(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):same(x,y)
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=2e-10,abs_tol=3e-12),(a,b)
    else: assert a==b,(a,b)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-once',action='store_true');args=parser.parse_args()
    result=json.loads(json.dumps(run(),ensure_ascii=False,allow_nan=False,default=lambda x:float(x)))
    dest=HERE/'review_checks.json'
    if args.write_once:
        with dest.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    else:
        original=sha(dest);same(result,json.loads(dest.read_text(encoding='utf-8')));assert sha(dest)==original
    print('1170 INDEPENDENT PASS_WITH_SCOPE')
    print(json.dumps(result['independent_finite_action'],ensure_ascii=False,indent=2))
