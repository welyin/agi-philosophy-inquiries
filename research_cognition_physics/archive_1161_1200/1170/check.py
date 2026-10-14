"""1170: original 288D cyclic support and conditional accelerated continuation.

No old files are written.  --write-once refuses replacement; default recalculates
and compares.  Uniform fields below are local comparison generators, not an
implemented accelerator or a globally positive gravitational model.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
SOURCES = {
 'archive_1106_1119/1112/check.py': '6e49d71ec5535fd6999432481938663a5543e2704db993ec7a00b8f475b6cd9f',
 'archive_1106_1119/1119/proof.md': 'b4a97154a8302a49c4f4f6788f615909b4849820bfa42e3bf11574c16808008d',
 'archive_935_955/research_note_941.md': 'a12e2dd67212a77538f850908efc6057e450d4def00dccd59aaa7cfaea346b7f',
 'archive_956_989/research_note_978.md': '0f1481e7e8602a6980855a3ec1a3efe7c18b945d9d03b79b10a6258580a304be',
 'archive_956_989/research_note_979.md': '2dafe3a18cedff7b87fe10744e5611ba90310cdd68086e62c6d85435be4f3b5f',
 'archive_1161_/1169/proof.md': '33f4384f5a851a5b8e35f816e577fc4101eefcf8365d9024bab0a70e739d5871',
 'archive_1161_/1169/acceptance.json': '0c1958cad480aaf001a4fd1eade7489a9871f1206893501b419c1915fa802819',
}
SIGMA, XBAR, AGE, G, LAMBDA, EPSILON, TIME = .5, 1., math.pi/2, .02, 1.3, .5, 1e-5


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm(x):
    return float(np.linalg.norm(x, 2))


def small(x, label, tol=3e-10):
    value = norm(x)
    assert value < tol, (label, value)
    return value


def moment(n):
    if n % 2:
        return 0.
    return SIGMA**n * math.prod(range(1, n, 2))


def polynomial_norm(poly):
    gram = sum(moment(j+k)*a.conj().T@b for j,a in enumerate(poly) for k,b in enumerate(poly))
    return float(math.sqrt(max(0., np.linalg.eigvalsh(gram)[-1])))


def mul_x(poly):
    out = np.zeros((len(poly)+1,)+poly.shape[1:], complex)
    out[:-1] += XBAR*poly
    out[1:] += poly
    return out


def derivative_wave(poly):
    """P times polynomial(z) times real Gaussian; common Gaussian factored out."""
    out = np.zeros((len(poly)+1,)+poly.shape[1:], complex)
    for j, a in enumerate(poly):
        if j:
            out[j-1] -= 1j*j*a
        out[j+1] += 1j*a/(2*SIGMA**2)
    return out


def h_poly(poly, charge, h, inv_mass):
    """H minus its common scalar rest baseline, on the Schwartz core."""
    out = np.zeros((len(poly)+2,)+poly.shape[1:], complex)
    for j,a in enumerate(poly):
        out[j] += h@a + G*XBAR*charge@a
        out[j+1] += G*charge@a
        out[j] += (j+.5)/SIGMA**2/2*inv_mass@a
        if j >= 2:
            out[j-2] -= j*(j-1)/2*inv_mass@a
        out[j+2] -= inv_mass@a/(8*SIGMA**4)
    return out


def envelope_h(poly, charge_bound, m_min):
    """Analytic coefficient bounds: ||h||<=3, ||M^-1||<=1/m_min."""
    out = np.zeros(len(poly)+2)
    for j,a in enumerate(poly):
        out[j] += (3+G*abs(XBAR)*charge_bound+(j+.5)/(2*m_min*SIGMA**2))*a
        out[j+1] += G*charge_bound*a
        if j>=2:
            out[j-2] += j*(j-1)*a/(2*m_min)
        out[j+2] += a/(8*m_min*SIGMA**4)
    return out


def envelope_norm(poly):
    return sum(float(a)*math.sqrt(moment(2*j)) for j,a in enumerate(poly))


def envelope_x(poly):
    out=np.zeros(len(poly)+1)
    out[:-1]+=abs(XBAR)*np.asarray(poly)
    out[1:]+=poly
    return out


def envelope_p(poly):
    out=np.zeros(len(poly)+1)
    for j,a in enumerate(poly):
        if j:
            out[j-1]+=j*a
        out[j+1]+=a/(2*SIGMA**2)
    return out


def calculate():
    for rel, expected in SOURCES.items():
        assert digest(BASE/rel)==expected, ('frozen source mismatch', rel)
    spec=importlib.util.spec_from_file_location('source1112', BASE/'archive_1106_1119/1112/check.py')
    old=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    material=old.Material(1)
    h, baseline=material.hfk, material.baseline
    mass=baseline*np.eye(288)+h
    inv_mass=np.linalg.inv(mass)
    columns=[]
    other=[]
    for head in range(4):
        for q in range(3):
            for dest, qr in [(columns,0),(other,1)]:
                bare=old.tensor(np.eye(4)[:,head],material.aux,np.eye(3)[:,q],np.eye(3)[:,qr])
                dest.append(material.gamma@bare)
    k, kp=np.column_stack(columns), np.column_stack(other)
    support=k@k.conj().T
    j=material.encode_joint(np.eye(3)[:,None,:],0).reshape(288,3)
    js=material.encode_joint(np.eye(3)[:,None,:],AGE).reshape(288,3)
    mask=np.zeros((4,8,3,3)); mask[:,[1,3,5,7],:,:]=1
    mr=np.diag(mask.ravel())
    krylov=np.column_stack([np.linalg.matrix_power(h,n)@j for n in range(4)])
    sv=np.linalg.svd(krylov,compute_uv=False)
    rank=int(np.count_nonzero(sv>1e-8))
    assert rank==12
    support_result={
      'full_internal_dimension':288,'prepared_columns':3,'cyclic_rank':rank,'unconstrained_dimension':276,
      'krylov_smallest_singular_value':float(sv[-1]), 'rest_baseline':baseline,
      'basis_isometry_error':small(k.conj().T@k-np.eye(12),'S basis'),
      'direct_fk_restriction_error':small(h@k-k@np.kron(material.hc,np.eye(3)),'full FK restriction'),
      'mr_restriction_error':small(mr@k-k@np.kron(np.diag([0,0,0,1]),np.eye(3)),'original Mr restriction'),
      'original_mr_free_probability':float((js.conj().T@mr@js)[0,0].real),
      'full_mr_mass_commutator_norm':norm(mass@mr-mr@mass),
    }
    assert support_result['full_mr_mass_commutator_norm']>.1
    # Exact energy projectors are fixed polynomials of the original 4-level hc.
    es=np.linalg.eigvalsh(material.hc)
    small(es-np.arange(4), 'head energy labels')
    e0=np.array([-1,math.sqrt(3),-math.sqrt(3),1])/math.sqrt(8)
    e1=np.array([math.sqrt(3),-1,-1,math.sqrt(3)])/math.sqrt(8)
    small(material.hc@e0,'explicit energy zero vector')
    small(material.hc@e1-e1,'explicit energy one vector')
    vhead=1j*(np.outer(e0,e1)-np.outer(e1,e0))
    v=k@np.kron(vhead,np.eye(3))@k.conj().T
    hidden=kp@np.kron(vhead,np.eye(3))@kp.conj().T
    visible_b=LAMBDA*mass+EPSILON*v
    hidden_b=LAMBDA*mass+EPSILON*hidden
    eig_diag=[]
    for e in range(4):
        pj=np.eye(288,dtype=complex)
        for f in range(4):
            if f!=e:
                pj=pj@(h-f*np.eye(288))/(e-f)
        eig_diag.append(small(pj@v@pj,'all eigenstate mean loophole'))
    # Original complete two-outcome retained record, not an averaged poststate.
    record=np.vstack([np.eye(288)-mr,mr])
    small(record.conj().T@record-np.eye(288),'retained record isometry')
    future_errors=[]
    for age in [.23,.71,1.31]:
        ja=material.encode_joint(np.eye(3)[:,None,:],age).reshape(288,3)
        future_errors.append(small((hidden_b-LAMBDA*mass)@ja,'hidden free prefix'))
        future_errors.append(small(np.block([[hidden,np.zeros((288,288))],[np.zeros((288,288)),hidden]])@record@ja,'hidden recorded future'))
    charge_result={
      'lambda':LAMBDA,'g':G,'acceleration':G*LAMBDA,'epsilon':EPSILON,
      'visible_noncommuting_norm':norm(mass@visible_b-visible_b@mass),
      'hidden_noncommuting_norm':norm(mass@hidden_b-hidden_b@mass),
      'visible_support_defect':norm((visible_b-LAMBDA*mass)@k),
      'hidden_support_defect':small((hidden_b-LAMBDA*mass)@k,'invisible complement'),
      'all_mass_eigenstate_mean_defects':eig_diag,
      'hidden_complete_future_defects':future_errors,
      'zero_branch_support_defect':small(hidden@k,'lambda zero remains allowed'),
      'visible_charge_positive_lower_bound':LAMBDA*baseline-EPSILON,
      'visible_js_norm':norm(v@js),
    }
    # Fixed finite witness. No spatial lattice and no replacement by <M>.
    aa=G*LAMBDA
    poly=js[None,:,:]
    delta=(G*EPSILON)*mul_x(np.array([v@js]))
    ha=h_poly(poly,LAMBDA*mass,h,inv_mass)
    haa=h_poly(ha,LAMBDA*mass,h,inv_mass)
    hb=h_poly(poly,visible_b,h,inv_mass)
    hbb=h_poly(hb,visible_b,h,inv_mass)
    had=h_poly(delta,LAMBDA*mass,h,inv_mass)
    mmax=baseline+3
    bound_a, bound_b=LAMBDA*mmax,LAMBDA*mmax+EPSILON
    a1=envelope_h([1.],bound_a,baseline)
    a2=envelope_h(a1,bound_a,baseline)
    b2=envelope_h(envelope_h([1.],bound_b,baseline),bound_b,baseline)
    # Exact ||V Js||=1/sqrt(2), independent of unknown q and free phase.
    denv=[G*EPSILON*abs(XBAR)/math.sqrt(2),G*EPSILON/math.sqrt(2)]
    adenv=envelope_h(denv,bound_a,baseline)
    A1,A2,B2,D0,AD=map(envelope_norm,[a1,a2,b2,denv,adenv])
    actual=[polynomial_norm(z) for z in [ha,haa,hbb,delta,had]]
    analytic=[A1,A2,B2,D0,AD]
    assert all(x<=y+1e-10 for x,y in zip(actual,analytic))
    # Relative output: project off the common initial spatial Gaussian.
    jet=SIGMA*G*EPSILON/math.sqrt(2)
    cstate=(A2+B2)/2+AD
    output_lower=TIME*jet-TIME**2*cstate
    # Same transported Mr: E_a=C_a Mr C_a^dagger in both descriptions.
    slope_matrix=1j*G*EPSILON*XBAR*js.conj().T@(v@mr-mr@v)@js
    slope=float(slope_matrix[0,0].real)
    small(slope_matrix-slope*np.eye(3),'held-out slope independent of q')
    assert abs(slope-G*EPSILON*XBAR*math.sqrt(3)/8)<1e-12
    # D_E bounds ||(E_a(t)-Mr) Delta Psi||, including both mass phases.
    XM=2*mmax*envelope_norm(envelope_x(denv))
    PP=2*envelope_norm(envelope_p(denv))
    MM=2*mmax*D0
    de=abs(aa)*TIME*XM+abs(aa)*TIME**2/2*PP+aa**2*TIME**3/6*MM
    r0=TIME**2*(A2+B2)/2
    probability_remainder=2*r0+(TIME*D0+r0)**2+2*TIME**2*A1*D0+2*TIME*de
    probability_lower=abs(slope)*TIME-probability_remainder
    assert output_lower>0 and probability_lower>0
    finite={
      'position_mean':XBAR,'position_sigma':SIGMA,'clock_age':AGE,'time':TIME,
      'momentum_rms':1/(2*SIGMA),'rms_velocity_upper':1/(2*SIGMA*baseline),
      'comparison_generator_only':True,'physical_accelerator_constructed':False,
      'phase_independent_jet':jet,'state_bound_coefficient':cstate,
      'complete_output_trace_distance_lower':output_lower,
      'transported_original_mr_slope':slope,
      'transported_original_mr_remainder':probability_remainder,
      'transported_original_mr_probability_gap_lower':probability_lower,
      'norms_from_exact_gaussian_moments':dict(zip(['Ha','Ha2','Hb2','Delta','HaDelta'],actual)),
      'analytic_upper_bounds':dict(zip(['Ha','Ha2','Hb2','Delta','HaDelta'],analytic)),
      'actual_two_run_total_probability_error_must_be_less_than':probability_lower,
      'actual_error_budget_certified':False,
    }
    for rel, expected in SOURCES.items():
        assert digest(BASE/rel)==expected, ('source changed during calculation',rel)
    return {'schema':'research_round1170_v1','round':1170,'status':'PASS_WITH_SCOPE',
      'new_adopted_axioms':0,'new_empirical_calibrations':0,'gr_hard_deliveries':0,
      'groups':{'original_cyclic_support':support_result,'charge_classification_and_weak_tests':charge_result,
                'fixed_finite_comparison_and_error_responsibility':finite},
      'source_hashes':SOURCES,'code_sha256':digest(Path(__file__))}


def compare(a,b,where='root'):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),where
        for key in a:
            compare(a[key],b[key],where+'.'+key)
    elif isinstance(a,list):
        assert len(a)==len(b),where
        for i,(x,y) in enumerate(zip(a,b)):
            compare(x,y,f'{where}[{i}]')
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=2e-10,abs_tol=2e-12),(where,a,b)
    else:
        assert a==b,(where,a,b)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-once',action='store_true')
    parser.add_argument('--preview',action='store_true')
    args=parser.parse_args()
    result=calculate()
    dest=HERE/'results.json'
    if args.write_once:
        with dest.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
    elif not args.preview:
        compare(result,json.loads(dest.read_text(encoding='utf-8')))
    print(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
