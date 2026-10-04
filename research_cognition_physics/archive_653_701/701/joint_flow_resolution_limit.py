"""701: joint scale/reflection limits of the actual free magnetic field from681.

All spatial lattice modes and their two physical polarizations are retained.
This is the declared Abelian comparison branch, not the original chiral/H_b
candidate. Analytic bounds in the note prove the limits; finite sums diagnose
the implementation. Only stdlib and the existing NumPy are used.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
import joint_gauge_flow_time_interface as prior

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_flow_resolution_limit_results.json'


def heat_defect(t,rho,omega,n=128):
    x,w=np.polynomial.legendre.leggauss(n)
    v=(x+1)/2
    return -math.sqrt(2*rho/math.pi)*float((w/2)@np.exp(
        -2*rho*omega**2*v*v-t*t/(8*rho*v*v)))


def L(z):
    return math.sqrt(2/math.pi)*math.exp(-z*z/8)+z*math.erf(z/math.sqrt(8))/2


def kernel_and_boundary_checks():
    identities=[]
    for rho,omega in ((.02,.7),(.08,1.),(.3,2.3)):
        for t in (0.,.1,.5,1.):
            closed=prior.heat_closed(t,rho,omega)
            zero=math.exp(-omega*t)/(2*omega)
            d=heat_defect(t,rho,omega)
            fine=heat_defect(t,rho,omega,192)
            assert abs(d-fine)<2e-12
            assert abs(closed-zero-d)<2e-12
            bound=math.sqrt(2*rho/math.pi)*math.exp(-t*t/(8*rho))
            assert -bound-2e-14<=d<=2e-14
            identities.append(dict(rho=rho,omega=omega,time=t,
                exact_identity_residual=abs(closed-zero-d),defect=d,bound=bound))
    z=.5;omega=1.
    A=L(2*z)-2*L(3*z)+L(4*z)
    lower=z*z*math.exp(-2*z*z)/math.sqrt(8*math.pi)
    assert A>lower>0
    witnesses=[]
    for rho in (1e-2,1e-3,1e-4,1e-6,1e-8):
        t=z*math.sqrt(rho)
        cov=np.array([prior.heat_closed(j*t,rho,omega) for j in (2,3,4)])
        raw=float(cov[0]-2*cov[1]+cov[2])
        scaled=raw/math.sqrt(rho)
        assert raw<0
        witnesses.append(dict(rho=rho,times=[t,2*t],source_amplitude=rho**(-.25),
            raw_reflection_norm=raw,rescaled_reflection_norm=scaled,
            strict_negative_limit=-A,limit_error=abs(scaled+A)))
    assert witnesses[-1]['limit_error']<3e-5
    return dict(heat_identity=identities,contact_witnesses=witnesses,
        strict_A_lower_bound=lower,A=A,
        contact_negative_limit_is_analytic_not_a_numerical_extrapolation=True)


def wick(cov,indices):
    if not indices:return 1.
    if len(indices)%2:return 0.
    first=indices[0]
    return sum(cov[first,indices[j]]*wick(cov,indices[1:j]+indices[j+1:])
               for j in range(1,len(indices)))


def polynomial_gram(same,cross):
    cov=np.block([[same,cross],[cross.T,same]])
    basis=((),(0,),(1,),(0,0),(0,1),(1,1))
    return np.array([[wick(cov,left+tuple(j+2 for j in right))
                      for right in basis] for left in basis])


def mode_sum(n,rho,times,lattice=True):
    # Midpoint link phases are removed unitarily; q_i=2 sin(a k_i/2)/a.
    # f_k=exp(-.15 |k|^2) v is a smooth fixed physical source, not a cutoff.
    a=2*math.pi/n;cap=(n-1)//2
    v=np.array([1.,.3,-.2]);same=np.zeros((len(times),len(times)))
    cross=same.copy();zero=cross.copy();S={j:0. for j in (1,2,3,4)}
    max_symbol_residual=0.
    for tup in itertools.product(range(-cap,cap+1),repeat=3):
        k=np.array(tup,dtype=float);r=float(np.linalg.norm(k))
        if r==0:continue  # magnetic k=0 is identically zero, not removed charge.
        q=2*np.sin(a*k/2)/a if lattice else k
        omega=float(np.linalg.norm(q));weight=math.exp(-.3*r*r)
        magnetic=omega**2*np.eye(3)-np.outer(q,q)
        assert np.linalg.norm(magnetic@q)<1e-10
        eigen=np.linalg.eigvalsh(magnetic)
        max_symbol_residual=max(max_symbol_residual,
            float(np.max(abs(eigen-np.array([0.,omega**2,omega**2])))))
        factor=float(v@magnetic@v)*weight
        for j in S:S[j]+=r**j*float(v@v)*weight
        for i,t in enumerate(times):
            for j,u in enumerate(times):
                def c(s):
                    return prior.heat_closed(s,rho,omega) if rho else math.exp(-omega*s)/(2*omega)
                same[i,j]+=factor*c(abs(t-u))
                cross[i,j]+=factor*c(t+u)
                zero[i,j]+=factor*math.exp(-omega*(t+u))/(2*omega)
    return same,cross,zero,S,max_symbol_residual


def joint_physical_sources():
    times=np.array([.3,.7]);coeff=np.array([1.,-.4]);tau=min(times);T=max(times)
    # Continuum calibration includes all modes through |k_i|<=16. The theorem
    # handles the infinite tail; this finite comparison is not a rigorous tail certificate.
    same,cross,_,sums,_=mode_sum(33,0.,times,lattice=False)
    target=float(coeff@cross@coeff)
    target_gram=polynomial_gram(same,cross)
    assert np.linalg.eigvalsh(target_gram)[0]>-1e-10
    assert np.linalg.norm(mode_sum(29,0.,times,lattice=False)[1]-cross)<1e-11
    bfactor=float(np.sum(abs(coeff))**2)
    Sj={j:sums[j]*bfactor for j in sums}
    rows=[]
    for n in (7,11,17,25):
        a=2*math.pi/n;rho=.15*a*a
        sam,ref,unflow,_,symbol_error=mode_sum(n,rho,times)
        value=float(coeff@ref@coeff);unflow_value=float(coeff@unflow@coeff)
        flow_bound=math.sqrt(2*rho/math.pi)*math.exp(-tau*tau/(2*rho))*Sj[2]
        lattice_bound=a*a*(Sj[3]/12+T*Sj[4]/24)+a**3*Sj[4]/(2*math.pi**3)
        assert abs(value-unflow_value)<=flow_bound+1e-10
        assert abs(value-target)<=flow_bound+lattice_bound+1e-9
        gram=polynomial_gram(sam,ref)
        # The full ordinary covariance remains positive; RP of each flowed
        # regulator is not asserted. All finite polynomial entries converge.
        assert np.linalg.eigvalsh(np.block([[sam,ref],[ref.T,sam]]))[0]>-1e-9
        rows.append(dict(spatial_sites=n**3,spacing=a,flow=rho,
            retained_nonzero_modes=n**3-1,physical_polarizations_per_mode=2,
            reflected_value=value,unflowed_value=unflow_value,
            continuum_calibration=target,error=abs(value-target),
            flow_defect_bound=flow_bound,joint_bound=flow_bound+lattice_bound,
            polynomial_gram_error=float(np.max(abs(gram-target_gram))),
            flowed_polynomial_gram_minimum=float(np.linalg.eigvalsh(gram)[0]),
            transverse_symbol_residual=symbol_error))
    assert rows[-1]['error']<rows[0]['error']
    assert rows[-1]['polynomial_gram_error']<rows[0]['polynomial_gram_error']
    return dict(rows=rows,source_moments={str(j):Sj[j] for j in Sj},
        continuum_calibration_spatial_modes=33**3-1,
        polynomial_basis=['1','B(f1)','B(f2)','B(f1)^2','B(f1)B(f2)','B(f2)^2'],
        full_magnetic_tensor_and_both_polarizations_retained=True,
        all_mode_limit_proof_uses_S4_and_explicit_tail_bound=True,
        quadrature_and_finite_sums_are_diagnostics_not_interval_certificates=True)


def run():
    deps=('research_note_382.md','research_note_383.md','research_note_384.md','research_note_386.md',
          'research_note_425.md','research_note_522.md','research_note_523.md','research_note_524.md',
          'research_note_646.md','research_note_677.md','research_note_681.md','research_note_699.md',
          'research_note_700.md','joint_gauge_flow_time_interface.py',
          'round701_drafts/joint_scale_entry.md','round701_drafts/entry_checks.json')
    return dict(date='2026-10-02',round=701,tests_run=2,failures=0,errors=0,
        heat_and_contact=kernel_and_boundary_checks(),
        full_magnetic_sources=joint_physical_sources(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(inherited_681_free_Abelian_physical_branch=True,
            full_spatial_modes_and_gauge_invariant_sources=True,
            continuous_physical_time_input=True,spatial_lattice_and_flow_joint_limit=True,
            fixed_smooth_polynomial_observables_have_RP_limit=True,
            shrinking_rescaled_sources_can_keep_negative_norm=True,
            original_Hb_Gauss_S9_chiral_candidate_continuum_unresolved=True,
            no_new_spatial_dimension_or_autonomous_instrument_claim=True,
            existing_523_524_instruments_not_reproved=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    elif TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=701,tests=2,all_checks_passed=True,
        contact_negative_limit=-result['heat_and_contact']['A'],
        final_source_error=result['full_magnetic_sources']['rows'][-1]['error'])))
