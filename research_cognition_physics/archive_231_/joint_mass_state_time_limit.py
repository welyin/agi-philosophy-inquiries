"""654: original full mass, candidate temporal Pfaffian, state and common limit.

This is an explicitly chosen physical mass insertion, not a claimed complete
four-dimensional overlap Yukawa prescription. No spatial continuum is tested.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as old
import joint_spinor_subgroup_mass as dictionary
import joint_gauss_fermion_influence as car
from joint_subgroup_measure_source import pfaffian

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_mass_state_time_limit_results.json'
PHIS=np.array([[.45,.31,-.17,.26,.62],[-.22,.53,.35,.12,-.41],
               [.37,-.19,.28,-.44,.55]])
W=dictionary.ph_matrix()
BLOCKS=[list(range(4*c,4*c+4))+[12+2*c,13+2*c,18+2*c,19+2*c]
        for c in range(3)]+[list(range(24,32))]


def norm(x):return float(np.linalg.norm(x,2))


def delta(phi):
    b=W@dictionary.bdg(phi)@W.conj().T
    assert norm(b[:32,:32])<1e-13
    return b[:32,32:]


def delta_source(phi,j):
    unit=np.eye(5)[j];f=old.original.F(phi)
    return (delta(unit)*np.sqrt(old.original.F(unit)/f)
            +delta(phi)*phi[j]/(6*f))


def creation(d):
    n=len(d);out=np.zeros((1<<n,1<<n),complex)
    for i in range(n):
        for j in range(i+1,n):
            if abs(d[i,j])==0:continue
            for s in range(1<<n):
                t,sign=car.word(s,[(j,True),(i,True)])
                if sign:out[t,s]+=sign*d[i,j]
    return out


def nil_exp(c,a):
    out=np.eye(len(c),dtype=complex);term=out.copy()
    for k in range(1,(len(c).bit_length()-1)//2+1):
        term=(-a/k)*(c@term);out+=term
    return out


def transfer(d,a):
    c=creation(d);l=nil_exp(c,a)
    return l@l.conj().T


def bdg_transfer(d,a):
    n=len(d);eye=np.eye(n)
    return np.block([[eye+a*a*d@d.conj().T,-a*d],[-a*d.conj().T,eye]])


def effective(d,a):
    vals,vec=np.linalg.eigh(bdg_transfer(d,a))
    assert min(vals)>0
    return (vec*(-np.log(vals)/a))@vec.conj().T


def covariance(b,beta):
    e,v=np.linalg.eigh(b)
    return (v/(1+np.exp(-beta*e)))@v.conj().T


def log_partition(d,a,beta):
    energies=2/a*np.arcsinh(a*np.linalg.svd(d,compute_uv=False)/2)
    return float(np.sum(np.logaddexp(beta*energies/2,-beta*energies/2)))


def raw_pf(ds,a):
    n=len(ds);r=len(ds[0]);s=np.zeros((n,n))
    for j in range(n):s[j,(j+1)%n]=-1 if j==n-1 else 1
    kinetic=np.kron((np.eye(n)-s.T)/2,np.eye(r))
    pair=np.zeros((n*r,n*r),complex)
    for j,d in enumerate(ds):pair[r*j:r*(j+1),r*j:r*(j+1)]=a*d/2
    anti=np.block([[pair,-kinetic.T],[kinetic,-pair.conj()]])
    sign=(-1)**((r*n)*(r*n+1)//2)
    return sign*pfaffian(anti)*2.**(r*n)


def conditional(ds,a):
    value=1.+0j
    for ids in (BLOCKS[0],BLOCKS[-1]):
        u=np.eye(256,dtype=complex)
        for d in ds:u=transfer(d[np.ix_(ids,ids)],a)@u
        z=np.trace(u)
        value*=z**(3 if ids==BLOCKS[0] else 1)
    return value


def pair(z):return [float(z.real),float(z.imag)]


def full_pfaffian_check():
    ds=[delta(p) for p in PHIS];rows=[];a=.37
    block_error=0.
    for d in ds:
        reconstructed=np.zeros_like(d)
        for ids in BLOCKS:reconstructed[np.ix_(ids,ids)]=d[np.ix_(ids,ids)]
        block_error=max(block_error,norm(reconstructed-d))
        for ids in BLOCKS[1:3]:
            block_error=max(block_error,norm(d[np.ix_(ids,ids)]-d[np.ix_(BLOCKS[0],BLOCKS[0])]))
    assert block_error<1e-13
    for n in (1,2,3):
        actual=raw_pf(ds[:n],a);expected=conditional(ds[:n],a)
        relative=float(abs(actual-expected)/abs(expected))
        assert relative<2e-11
        rows.append(dict(time_sites=n,normalized_Pfaffian=pair(actual/2**32),
                         normalized_CAR_trace=pair(expected/2**32),relative_error=relative))
    # Exact covariance in all 32 modes, not a fixed-background Gauss projection.
    rng=np.random.default_rng(654);errors=[]
    for _ in range(5):
        c=old.gauge.group_exp(rng.normal(size=8)*.2,3)
        w=old.gauge.group_exp(rng.normal(size=3)*.2,2);z=np.exp(.2j*rng.normal())
        r=np.kron(dictionary.left_rep(c,w,z),np.eye(2))
        rn=np.block([[r,np.zeros_like(r)],[np.zeros_like(r),r.conj()]])
        p=PHIS[0];h=z**3*w@(p[:2]+1j*p[2:4]);pnew=np.r_[h.real,h.imag,p[4]]
        errors.append(norm(bdg_transfer(delta(pnew),a)-rn@bdg_transfer(ds[0],a)@rn.conj().T))
    assert max(errors)<2e-12
    return dict(actual_all32_mode_rows=rows,block_reconstruction_error=block_error,
                original_group_covariance_error=max(errors),original_complex_Y_preserved=True,
                time_dependent_conditional_trace_need_not_be_real=True,
                nonzero_scalar_mass_insertion_is_declared_candidate=True)


def full_state_source_check():
    p=PHIS[0];d=delta(p);zero=np.zeros_like(d)
    b=np.block([[zero,d],[d.conj().T,zero]]);beta=2.
    target=covariance(b,beta);rows=[]
    for a in (1.,.5,.25,.125):
        ba=effective(d,a)
        first=np.block([[-a*d@d.conj().T/2,d],[d.conj().T,a*d.conj().T@d/2]])
        s=np.linalg.svd(d,compute_uv=False);e=2/a*np.arcsinh(a*s/2)
        eig_expected=np.sort(np.r_[e,-e])
        spectrum_error=float(np.max(abs(np.linalg.eigvalsh(ba)-eig_expected)))
        assert spectrum_error<2e-13
        rows.append(dict(time_step=a,quasiparticle_energy_error=float(max(abs(e-s))),
            covariance_operator_error=norm(covariance(ba,beta)-target),
            corrected_BdG_remainder=norm(ba-first),spectrum_identity_error=spectrum_error))
    assert all(rows[i+1]['covariance_operator_error']<.55*rows[i]['covariance_operator_error'] for i in range(3))
    assert all(rows[i+1]['corrected_BdG_remainder']<.27*rows[i]['corrected_BdG_remainder'] for i in range(3))
    # Independently insert the derivative into the complete original block traces.
    sources=[];a=.5;n=4;step=2e-5
    ts={};cs={}
    for ids in (BLOCKS[0],BLOCKS[-1]):
        key=tuple(ids);ts[key]=transfer(d[np.ix_(ids,ids)],a)
    for j in range(5):
        dd=delta_source(p,j);source=0.+0j
        for ids in (BLOCKS[0],BLOCKS[-1]):
            t=ts[tuple(ids)];dc=creation(dd[np.ix_(ids,ids)])
            dt=-a*(dc@t+t@dc.conj().T)
            term=n*np.trace(np.linalg.matrix_power(t,n-1)@dt)/np.trace(np.linalg.matrix_power(t,n))
            source+=(3 if ids==BLOCKS[0] else 1)*term
        unit=np.eye(5)[j]*step
        fd=(log_partition(delta(p+unit),a,beta)-log_partition(delta(p-unit),a,beta))/(2*step)
        exact_derivative_error=norm(dd-(delta(p+unit)-delta(p-unit))/(2*step))
        assert abs(source.imag)<1e-11 and abs(source.real-fd)<2e-8 and exact_derivative_error<2e-9
        sources.append(dict(component=j,transfer_insertion=float(source.real),spectral_derivative=fd,
            insertion_error=float(abs(source.real-fd)),original_mass_derivative_error=exact_derivative_error))
    # Independent density matrix on the original neutral Dirac/Majorana block.
    ids=[24,25,30,31];dp=delta(np.array([.45,0.,0.,0.,.62]))[np.ix_(ids,ids)]
    tp=transfer(dp,.5);rho=np.linalg.matrix_power(tp,4);rho/=np.trace(rho)
    ops=[]
    for j in range(4):
        op=np.zeros((16,16),complex)
        for s in range(16):
            target,sign=car.word(s,[(j,False)])
            if sign:op[target,s]=sign
        ops.append(op)
    ops+=[op.conj().T for op in ops]
    actual=np.array([[np.trace(rho@x@y.conj().T) for y in ops] for x in ops])
    density_error=norm(actual-covariance(effective(dp,.5),2.))
    assert density_error<2e-12
    return dict(all32_state_rows=rows,all5_original_scalar_sources=sources,
                independent_four_mode_density_covariance_error=density_error,
                eigenvalue_convergence_does_not_imply_equal_finite_step_state=True,
                common_limit_proved_on_compact_scalar_sets_not_full_configuration_domain=True)


def neutral_state_check():
    p=np.array([0.,0.,0.,0.,1.2]);mu=delta(p)[30,31];d=mu*dictionary.EPS
    c=creation(d);h=c+c.conj().T;beta=2.;rows=[]
    for a in (1.,.5,.25):
        n=int(beta/a);t=transfer(d,a);e,v=np.linalg.eigh(t)
        ha=(v*(-np.log(e)/a))@v.conj().T
        rho=np.linalg.matrix_power(t,n);rho/=np.trace(rho)
        energy=2/a*np.arcsinh(a*abs(mu)/2)
        naive=car.exp_h(h*energy/abs(mu),beta);naive/=np.trace(naive)
        exact=car.exp_h(h,beta);exact/=np.trace(exact)
        # The diagonal O(a) term includes its vacuum constant.
        correction=h+a*(c.conj().T@c-c@c.conj().T)/2
        identity=4**-n*(2+2*np.cosh(beta*energy))
        pf=raw_pf([d]*n,a)/4**n
        assert abs(pf-identity)<1e-13
        rows.append(dict(time_step=a,original_mass=abs(mu),effective_energy=float(energy),
            exact_energy_only_refit_state_trace_distance=float(np.sum(abs(np.linalg.eigvalsh(rho-naive)))/2),
            original_state_trace_distance=float(np.sum(abs(np.linalg.eigvalsh(rho-exact)))/2),
            first_order_H_remainder=norm(ha-correction),
            vacuum_probability_shift=float((rho-naive)[0,0].real)))
    assert min(r['exact_energy_only_refit_state_trace_distance'] for r in rows)>1e-3
    return dict(original_neutral_rows=rows,finite_step_Bogoliubov_rotation_would_change_sources=True,
                no_state_match_inferred_from_only_partition_spectrum=True)


def harmonic_characters():
    weights=[-2,-2,-2,3,3,2,2,2,-3,-3]
    h=[{0:1}]+[{} for _ in range(8)]
    for w in weights:
        for n in range(1,9):
            for q,m in list(h[n-1].items()):h[n][q+w]=h[n].get(q+w,0)+m
    out=[]
    for l in range(9):
        c=h[l].copy()
        if l>=2:
            for q,m in h[l-2].items():c[q]=c.get(q,0)-m
        assert all(m>=0 for m in c.values())
        out.append(c)
    return out


def common_auxiliary_limit_check():
    saved=json.loads((HERE/'joint_chiral_character_state_results.json').read_text('utf8'))
    spec=saved['temporal_gluing']['spherical_harmonic_spectrum']
    lam0=Fraction(spec[0]['eigenvalue']);ratios=[Fraction(x['eigenvalue'])/lam0 for x in spec]
    assert ratios[1]==Fraction(8,17) and all(ratios[i+1]<ratios[i] for i in range(8))
    chars=harmonic_characters();assert [sum(c.values()) for c in chars]==[x['dimension'] for x in spec]
    rows=[];theta=.21;beta=2.
    for n in (2,4,8,16,32,64):
        a=beta/n;tail=sum(x['dimension']*r**n for x,r in zip(spec[1:],ratios[1:]))
        bound=35749*ratios[1]**n;assert tail<=bound
        value=0j;derivative=0j
        for r,c in zip(ratios,chars):
            value+=float(r**n)*sum(m*np.exp(1j*q*theta) for q,m in c.items())
            derivative+=float(r**n)*sum(1j*q*m*np.exp(1j*q*theta) for q,m in c.items())
        assert abs(value-1)<=float(tail)+1e-13 and abs(derivative)<=24*float(tail)+1e-13
        rows.append(dict(time_sites=n,time_step=a,auxiliary_excited_trace=float(tail),
            auxiliary_state_trace_distance=float(tail/(1+tail)),uniform_tail_bound=float(bound),
            first_auxiliary_gap=float(np.log(17/8)/a),holonomy_weight=pair(value),
            holonomy_derivative=pair(derivative)))
    assert rows[-1]['auxiliary_state_trace_distance']<2e-20
    return dict(lambda0=str(lam0),largest_excited_ratio='8/17',rows=rows,
                exact_singlet_low_energy_projection=True,
                finite_gap_nontrivial_auxiliary_dynamics_not_preserved_by_fixed_kernel_refinement=True,
                discarded_vacuum_energy_per_step=float(32*np.log(2)-np.log(float(lam0))),
                vacuum_subtraction_not_a_solution_of_gravitational_energy=True)


def run():
    deps=('research_note_598.md','research_note_612.md','research_note_614.md','research_note_615.md',
          'research_note_643.md','research_note_646.md','research_note_653.md',
          'joint_fermion_gauss_completion.py','joint_spinor_subgroup_mass.py',
          'joint_gauss_fermion_influence.py','joint_subgroup_measure_source.py',
          'joint_chiral_character_state_results.json')
    return dict(date='2026-10-02',round=654,tests_run=4,failures=0,errors=0,
        full_mass_pfaffian=full_pfaffian_check(),state_and_sources=full_state_source_check(),
        neutral_state=neutral_state_check(),common_auxiliary_limit=common_auxiliary_limit_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Declared spatially trivial physical mass insertion: all original 32 modes and complex masses, positive temporal transfer, common state/source limit on compact scalar backgrounds. Original auxiliary kernel projects to its singlet in that limit. Not full spatial/chiral/GR equivalence.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
