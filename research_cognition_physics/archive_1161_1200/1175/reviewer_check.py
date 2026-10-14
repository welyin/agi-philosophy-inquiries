"""1175 independent audit. Default: read-only recomputation; --write: new JSON.

Dynamics come from the reviewer's frozen Taylor implementation, written from
the candidate equations before reading the author's preview implementation.
The author is imported only after independent computation, for comparison.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
import argparse
from decimal import Decimal, localcontext
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ARCHIVE=HERE.parent
SCREEN=ARCHIVE/'_admission'/'after1174_direction'
EXPECTED={
 '../research_note_1175.md':'c3821f148d9b2c9b3db3e6cffc8988d88591bd40f33533b1feb46674f51417d1',
 'admission.md':'35a4124c019f094945210e019844b61e5d47f42dc695d1c409b8ddf4fd410ef3',
 'proof.md':'581427de3c8e03176a1145773d6bd0754930d849524110a22971c0e3f8e861f4',
 'check.py':'1438ef6f89890cf7d4e08b6aa8e3bd9efbe2feeab5a7d43c04ff4306d5de0a54',
 'results.json':'537c8d179be1a46ae7f9385227801f74ac92a09ea486cb07e6b114866d631254',
 'NEXT.md':'5cce5bec4f25d117ac45d0b6fc58283db2055e5565e48ee34e84a38022c56432',
}
OWN_PREVIEW_SHA='a69307d6ca8c254a92168451031bfba3650c635ebbfe74f92ec2d8554bd3a4a2'


def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def module(name,path):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def fourier(axis):
    f=np.eye(4,dtype=complex)/2
    f[0,axis+1]=.25j;f[axis+1,0]=-.25j
    return f


def assemble(blocks):
    u=np.zeros((48,48),complex)
    for p,m in itertools.product(range(4),range(3)):
        i=12*p+4*m;u[i:i+4,i:i+4]=blocks[p,m]
    return u


def taylor32(h):
    total=np.broadcast_to(np.eye(4,dtype=complex),h.shape).copy()
    term=total.copy()
    for k in range(1,33):
        term=(-1j*h)@term/k
        total+=term
    return total


def decimal_resource_and_symbol(k_stored):
    """65-digit independent constants; no huge integer K is instantiated."""
    with localcontext() as ctx:
        ctx.prec=65
        w=Decimal(80245425)/16
        ll=Decimal(709459359375)/100000
        logarithm=(4*ll*2).ln()+6*w*2-Decimal(52275242)*Decimal(10).ln()/2
        sufficient=2*((4*ll*2/Decimal('1e-4')).ln()+6*w*2)/Decimal(10).ln()
        pi_lo=Decimal.from_float(3.141592653589793)
        pi_hi=Decimal.from_float(3.1415926535897936)
        ks=[]
        for pi in (pi_lo,pi_hi):
            x=2*pi/65;term=x;value=x
            for j in range(1,22):
                term=-term*x*x/Decimal((2*j)*(2*j+1));value+=term
            # Next term is < 1e-98 on this small interval; decimal guard dominates.
            ks.append(value*65/16)
        guard=Decimal('1e-60')
        lo,hi=ks[0]-guard,ks[1]+guard
        stored=Decimal.from_float(float(k_stored))
        physical_delta=max(abs(stored-lo),abs(stored-hi))
        return {'w_exact':str(w),'L_exact':str(ll),
                'sufficient_log10_K':str(sufficient),'chosen_etaK_bound':str(logarithm.exp()),
                'physical_low_k_interval_decimal':[str(lo),str(hi)],
                'stored_low_k_absolute_error_bound':float(physical_delta)}


def full_instrument_test(u,helper):
    """48-code + R=3, all sites and three internal projective outcomes."""
    rng=np.random.default_rng(11750048)
    psi=rng.normal(size=(48,3))+1j*rng.normal(size=(48,3))
    psi/=np.linalg.norm(psi)
    out=u@psi
    rhoR=psi.conj().T@psi
    propagatedR=out.conj().T@out
    q0=np.diag([1.,0.,0.]).astype(complex)
    qm=np.zeros((3,3),complex);qm[1:,1:]=helper.PCLK
    qp=np.eye(3)-q0-qm
    qs=[q0,qm,qp]
    assert np.max(np.abs(sum(qs)-np.eye(3)))<1e-15
    rho_after=np.zeros((3,3),complex)
    total=0.;minprob=1.;smooth=np.zeros((3,3))
    n=65;volume=n**3;mom=out.reshape(4,3,4,3)
    for start in range(0,volume,8192):
        ids=np.arange(start,min(start+8192,volume))
        sites=np.stack((ids//(n*n),(ids//n)%n,ids%n),axis=1)
        phase=np.exp(2j*np.pi*(sites@helper.MODES.T)/n)/np.sqrt(volume)
        amp=np.einsum('jp,pmdr->jmdr',phase,mom)
        f=(1+np.sin(2*np.pi*sites/n))/2
        for qi,q in enumerate(qs):
            projected=np.einsum('mn,jndr->jmdr',q,amp)
            probs=np.sum(np.abs(projected)**2,axis=(1,2,3))
            total+=probs.sum();minprob=min(minprob,float(probs.min()))
            rho_after+=np.einsum('jmdr,jmds->rs',projected.conj(),projected)
            smooth[:,qi]+=f.T@probs
    expected=np.zeros_like(smooth)
    for axis in range(3):
        for qi,q in enumerate(qs):
            effect=u.conj().T@np.kron(fourier(axis),np.kron(q,np.eye(4)))@u
            expected[axis,qi]=np.trace(psi.conj().T@effect@psi).real
    # A coherent cross-mass block survives propagation; no mass superselection.
    before=psi.reshape(4,3,4,3).transpose(1,0,2,3).reshape(3,-1)
    after=out.reshape(4,3,4,3).transpose(1,0,2,3).reshape(3,-1)
    coh0=np.linalg.norm(np.outer(before[0],before[1].conj()))
    coh1=np.linalg.norm(np.outer(after[0],after[1].conj()))
    ans={'reference_dimension':3,'fine_outcomes':volume*3,
         'joint_record_normalization_error':float(abs(total-1)),
         'minimum_joint_probability':minprob,
         'R_error_after_propagation':float(np.linalg.norm(propagatedR-rhoR,2)),
         'R_error_after_complete_nonselective_instrument':float(np.linalg.norm(rho_after-rhoR,2)),
         'site_instrument_vs_compressed_effect_max':float(np.max(np.abs(smooth-expected))),
         'initial_cross_mass_coherence_norm':float(coh0),
         'cross_mass_coherence_norm_change':float(abs(coh1-coh0)),
         'note':'Finite diagnostic R and complete instrument test; universal R guarantee comes from operator norm proof.'}
    assert ans['joint_record_normalization_error']<1e-11
    assert ans['R_error_after_complete_nonselective_instrument']<1e-11
    assert ans['site_instrument_vs_compressed_effect_max']<1e-11
    assert coh0>.1 and abs(coh1-coh0)<1e-12
    return ans


def run():
    for name,expected in EXPECTED.items():
        assert digest(HERE/name)==expected,('author core changed',name)
    saved=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
    for name,expected in saved['sources'].items():
        assert digest(ARCHIVE/name)==expected,('dependency changed',name)
    own=SCREEN/'mean_field_review'/'reviewer_preview.py'
    assert digest(own)==OWN_PREVIEW_SHA
    helper=module('independent_frozen_taylor1175',own)
    independent=helper.run()
    helper.compare(independent,json.loads(own.with_suffix('.json').read_text(encoding='utf-8')))
    refined=helper.taylor_run(order=14,steps=80)
    # Exact massless block is obtained from the same independent xi trajectory.
    u0=helper.flat_blocks(refined['xi'],massless=True)[:,:1]
    three=np.concatenate((u0,refined['U'][0]),axis=1)
    independent48=assemble(three)
    instrument=full_instrument_test(independent48,helper)

    print('Independent trajectory finished; compare frozen author endpoint data',flush=True)
    author=module('author_only_after_independent1175',SCREEN/'mean_field_candidate'/'preview.py')
    actual=author.evolve(author.NU)
    masses=np.array([0.,.3,.5])
    h=actual['tau']*author.KIN[:,None]+actual['tau']*masses[None,:,None,None]*author.BETA
    rows=float(np.max(np.sum(np.abs(h),axis=-1)))
    assert rows<2
    flat_taylor=taylor32(h)
    flat_library=author.exact_const(actual['tau'],actual['tau'],masses)
    zero_h=actual['xi']*author.KIN[:,None]
    zero_taylor=taylor32(zero_h)
    zero_library=author.exact_const(0,actual['xi'],np.zeros(1))
    flat_delta=max(np.linalg.norm(a-b,2) for a,b in zip(
        flat_taylor.reshape(-1,4,4),flat_library.reshape(-1,4,4)))
    zero_delta=max(np.linalg.norm(a-b,2) for a,b in zip(
        zero_taylor.reshape(-1,4,4),zero_library.reshape(-1,4,4)))
    fourier_delta=max(np.linalg.norm(fourier(j)-author.site_effect(j),2) for j in range(3))
    resource=decimal_resource_and_symbol(float(author.LMOM[1,0]))

    # Exact massive witness: sixteen coefficients 1/4, no sqrt evaluation.
    massive=np.zeros((4,3,4),complex);massive[:,1:,0]=.25;massive[:,1:,3]=.25
    massive=massive.ravel()
    zero=np.zeros((4,3,4),complex);zero[:,0,:]=helper.CHI/2;zero=zero.ravel()
    clock=np.zeros((3,3),complex);clock[1:,1:]=helper.PCLK
    author_three=np.concatenate((zero_library,actual['u']),axis=1)
    ua=assemble(author_three);uf=assemble(flat_taylor)
    affine_probabilities=[]
    for axis in range(3):
        for c,wit in ((np.eye(3),massive),(clock,massive),(np.eye(3),zero)):
            effect=np.kron(fourier(axis),np.kron(c,np.eye(4)))
            affine_probabilities.append([float(np.vdot(uu@wit,effect@(uu@wit)).real)
                                         for uu in (ua,uf)])
    paper=saved['group3_fixed_menu_and_error_certificate']
    endpoint_y=affine_probabilities[4]
    endpoint_m0=affine_probabilities[2]
    endpoint_probability_difference=max(
        abs(endpoint_y[0]-paper['fixed_y_clock_site_probabilities']['active']),
        abs(endpoint_y[1]-paper['fixed_y_clock_site_probabilities']['proper_time_flat']),
        abs(endpoint_m0[0]-paper['fixed_massless_x_probabilities']['active']),
        abs(endpoint_m0[1]-paper['fixed_massless_x_probabilities']['proper_time_flat']))
    u=2.**-53
    gamma=lambda k:k*u/(1-k*u)
    # Standard fixed-size complex product/sum count; exp(2)<7.39.
    taylor_round=2*gamma(1024)*7.39
    effect_round=64*gamma(4096)
    input_error=4*resource['stored_low_k_absolute_error_bound']+2e-14
    endpoint_guard=1e-10
    arithmetic_sum=taylor_round+effect_round+input_error+2*(flat_delta+zero_delta+fourier_delta)+endpoint_probability_difference
    tail=2.**33/math.factorial(33)/(1-2/34)
    assert arithmetic_sum+tail<endpoint_guard
    num=paper['numerical_enclosure']
    residual=num['probe_own_interpolation_residual_error']
    # Add a deliberately conservative extra 3 delta, beyond the one already paid.
    normalization_extra=3*residual
    assert endpoint_guard+normalization_extra<1e-9
    uniform_numeric=paper['uniform_48_code_numeric_sum_bound']
    assert uniform_numeric<9e-5
    # Recompute the budget with exact decimal arithmetic on archived binary data.
    with localcontext() as ctx:
        ctx.prec=65
        D=Decimal.from_float
        lat=D(paper['lattice_dstar_upper'])*(D(num['stored_tau'])+D(num['clock_tau_error'])+
                                             D(num['stored_xi'])+D(num['clock_xi_error']))
        total=lat+Decimal('.0002')+Decimal('.00004')+Decimal('.00009')
        # Archived total has outward rounding; allow no downward discrepancy.
        assert D(paper['total_two_branch_probability_budget'])>=total
        ly=abs(D(paper['fixed_y_clock_site_probabilities']['signed_gap']))-total
        l0=abs(D(paper['fixed_massless_x_probabilities']['signed_gap']))-total
        assert ly>0 and l0>0
        budgets={'decimal_pair_lattice':str(lat),'decimal_pair_total':str(total),
                 'decimal_y_lower':str(ly),'decimal_massless_lower':str(l0)}
    # 35 scalar comparisons to the independently computed frozen preview.
    comparisons=[]
    for key in ('tau','xi'):
        comparisons.append(abs(independent['active'][key]-saved['group2_autonomous_common_transport']['background'][key]))
    for axis,letter in enumerate('xyz'):
        for typ in ('position','joint'):
            key=letter+('_position' if typ=='position' else '_joint_clock_minusY')
            rr=independent['witness']['lattice'];aa=saved['group2_autonomous_common_transport']['all_three_mobile_clock_effects'][key]
            comparisons.extend([abs(rr['active'][typ][axis]-aa['active']),
                                abs(rr['flat_same_tau'][typ][axis]-aa['proper_time_flat']),
                                abs(rr[typ+'_gap'][axis]-aa['signed_gap']),
                                abs(independent['unknown32_effect_difference_opnorm']['lattice'][2*axis+(typ=='joint')]-aa['complete_input_effect_distance'])])
        aa=saved['group4_delimited_controls']['massless_all_axes'][letter+'_position']
        comparisons.extend([abs(independent['witness']['lattice']['massless_active']['position'][axis]-aa['active']),
                            abs(independent['witness']['lattice']['massless_flat_same_tau']['position'][axis]-aa['proper_time_flat']),
                            abs(independent['witness']['lattice']['massless_gap'][axis]-aa['signed_gap'])])
    assert max(comparisons)<1e-10
    answer={'round':1175,'status':'ACCEPT_WITH_SCOPE','blockers':0,
      'author_six_core_sha256':EXPECTED,'author_dependency_sha256':saved['sources'],
      'reviewer_taylor_dependency_sha256':OWN_PREVIEW_SHA,
      'independence':'Taylor implementation preceded reading author preview.py. This is AI internal independent implementation, not blind or external review.',
      'group1_parent_resources_and_scope':resource,
      'group2_independent_dynamics':{'refinement':independent['refinement'],
        'full_small_grid_Dirac':independent['small_grid_full_dirac'],
        'author_scalar_comparisons':len(comparisons),'max_author_scalar_difference':max(comparisons),
        'author_RK4_vs_independent_Taylor_matrix_max':float(np.max(np.abs(actual['u']-refined['U'][0]))),
        'independent_active':independent['active'],'lapse_control':independent['lapse_flat_max']},
      'group3_48_code_records_and_endpoint':{'full_48_gram_error':float(np.linalg.norm(independent48.conj().T@independent48-np.eye(48),2)),
        'complete_instrument':instrument,'flat_H_max_row_sum':rows,'Taylor32_tail_bound':tail,
        'Taylor32_basic_roundoff_bound':taylor_round,'effect_basic_roundoff_bound':effect_round,
        'physical_coefficient_allowance':input_error,'flat_library_vs_Taylor32':float(flat_delta),
        'massless_library_vs_Taylor32':float(zero_delta),'author_Fourier_vs_analytic':float(fourier_delta),
        'endpoint_probability_recompute_difference':float(endpoint_probability_difference),
        'endpoint_arithmetic_sum_bound':float(arithmetic_sum+tail),'endpoint_guard':endpoint_guard,
        'extra_nonunitary_probability_allowance':normalization_extra,
        'endpoint_plus_normalization_allowance':endpoint_guard+normalization_extra,
        'published_arithmetic_allowance':1e-9},
      'group4_budget_and_stopping_lines':{**budgets,
        'numerical_probability_budget':9e-5,'underlying_uniform_numeric_bound':uniform_numeric,
        'y_strict_but_narrow_margin':paper['certified_y_joint_gap_lower'],
        'massless_margin':paper['certified_massless_x_gap_lower'],
        'finite_K_full_BG_and_complete_discrete_instrument':True,
        'continuum_claim_only_prepared_low_band_smooth_terminal_menu':True,
        'actual_instrument_precision_is_input':True,'formal_hardware_verification':False,
        'Einstein_or_dark_energy_claim':False,
        'note':'Numerical Hermite/Bernstein envelope is code/proof audited and separately reproduced by author default check; this script independently computes Taylor dynamics, endpoints and budget.'}}
    return helper.encode(answer)


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        return max((compare(a[k],b[k]) for k in a),default=0.)
    if isinstance(a,list):
        assert len(a)==len(b)
        return max((compare(x,y) for x,y in zip(a,b)),default=0.)
    if isinstance(a,(float,int)) and not isinstance(a,bool):
        err=abs(a-b);assert err<=2e-11*max(1,abs(a),abs(b));return err
    assert a==b;return 0.


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=run();output=HERE/'review_checks.json'
    if args.write:
        with output.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
        print('WRITE_ONCE',result['status'],flush=True)
    else:
        print('READ_ONLY_MATCH',compare(result,json.loads(output.read_text(encoding='utf-8'))),flush=True)
    print(json.dumps({'status':result['status'],'blockers':result['blockers'],
          'y_lower':result['group4_budget_and_stopping_lines']['y_strict_but_narrow_margin'],
          'endpoint_allowance':result['group3_48_code_records_and_endpoint']['endpoint_plus_normalization_allowance']},indent=2),flush=True)
