"""672: off-diagonal signed matter histories and original-group electric sewing.

Boundary-conditioned finite diagnostics, not a Gauss-integrated no-go theorem.
The analytic sewing criterion and series bounds are in research_note_672.md.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_nonflat_mass_measure as prior
old=prior.old;internal=prior.internal;mass=prior.mass
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_gauge_history_kernel_results.json'


def pack(z):return dict(real=z.real.tolist(),imag=z.imag.tolist())


def configuration(left,right,lam=.37,audit=False,rephase=False):
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    links[0,:2]=left;links[0,2:]=right
    u,v,d,h=prior.kernel(links)
    if rephase:
        u=u*np.exp(1j*np.linspace(-.3,.4,128))
        v=v*np.exp(1j*np.linspace(.2,-.7,128))
    jm=np.kron(np.kron(np.eye(4),internal.VM),np.eye(16))
    jp=np.kron(np.kron(np.eye(4),internal.VP),np.eye(16))
    m=np.kron(np.eye(4),np.kron(internal.B,internal.T[0]))
    mb=np.kron(np.eye(4),np.kron(internal.B,internal.T[0].conj().T))
    pair,_=mass.mass_pairing(4,mass.car.PHI[0])
    q=prior.joint(u,v,d,jm,jp,m,mb,pair,lam=lam,local_check=audit)
    reflect=np.eye(4)[[2,3,0,1]];tr=np.kron(reflect,np.eye(32))
    theta=np.block([[np.zeros_like(tr),tr],[tr,np.zeros_like(tr)]])
    idx=np.r_[np.arange(64,128),np.arange(192,256)]
    q['cross']=(theta@q['covariance'])[np.ix_(idx,idx)]
    q['gap']=float(min(abs(np.linalg.eigvalsh(h))))
    q['frame_phase']=np.linalg.det(q['S'])
    return q


def gram_checks(k):
    herm=old.err(k-k.conj().T)
    assert herm<2e-11
    kh=(k+k.conj().T)/2
    ev,vec=np.linalg.eigh(kh);w=vec[:,0]
    return dict(kernel=pack(k),hermiticity_error=herm,eigenvalues=ev.tolist(),
        witness=pack(w),rayleigh=old.cpair(np.vdot(w,k@w)),
        eigen_residual=old.norm(k@w-ev[0]*w))


def family(links):
    rng=np.random.default_rng(67282)
    packets=rng.normal(size=(128,3))+1j*rng.normal(size=(128,3))
    packets/=np.linalg.norm(packets,axis=0)
    scalar=np.zeros((3,3),complex);linear=np.zeros_like(scalar)
    gaps=[];phases=[];raw=[];diagonal_minima=[]
    for i in range(3):
        for j in range(3):
            q=configuration(links[i],links[j]);gaps.append(q['gap'])
            scalar[i,j]=q['weight']
            linear[i,j]=q['weight']*np.vdot(packets[:,i],q['cross']@packets[:,j])
            raw.append(old.cpair(q['weight']));phases.append(old.cpair(q['frame_phase']))
            if i==j:
                assert old.err(q['cross']-q['cross'].conj().T)<2e-12
                diagonal_minima.append(float(min(np.linalg.eigvalsh((q['cross']+q['cross'].conj().T)/2))))
    assert max(abs(np.diag(scalar).imag))<1e-12*max(abs(np.diag(scalar)))
    assert min(np.diag(scalar).real)>0
    scale=np.sqrt(np.diag(scalar).real)
    k=scalar/scale[:,None]/scale[None,:]
    w=linear/scale[:,None]/scale[None,:]
    assert min(np.linalg.eigvalsh((k+k.conj().T)/2))<-.1
    checked=configuration(links[0],links[2],audit=True)
    phase=configuration(links[0],links[2],rephase=True)
    errors=dict(direct_Pfaffian=checked['local_weight_error'],
        direct_physical_source=checked['local_covariance_error'],
        triangular=checked['triangular_error'],
        frame_rephasing=float(abs(phase['weight']/checked['weight']-1)))
    assert max(errors.values())<2e-10
    return k,dict(scalar=gram_checks(k),linear_sample=gram_checks(w),
        minimum_Wilson_gap=min(gaps),raw_signed_weights=raw,frame_determinants=phases,
        diagonal_physical_Gram_minima=diagonal_minima,identity_errors=errors,
        no_entrywise_absolute_value_or_pairwise_renormalization=True)


def tail(a,power,n):
    """Upper bound sum_{j>n}(j+1)^power exp(-a*j*j), power 0,2,4."""
    assert n>=max(1,math.sqrt(power/(2*a)))
    integ=math.sqrt(math.pi)/(2*math.sqrt(a))*math.erfc(math.sqrt(a)*n)
    for k in range(2,power+1,2):
        integ=n**(k-1)*math.exp(-a*n*n)/(2*a)+(k-1)*integ/(2*a)
    return 2**power*integ


def series_bound(a,power,n):
    return sum((j+1)**power*math.exp(-a*j*j) for j in range(n+1))+tail(a,power,n)


def electric_kernel(tau,delta):
    """Bi-invariant electric heat on original (SU3 x SU2 x U1)/Z6.

    Original integer charges, Casimirs C3 and l(l+2)/4, and circle n^2.
    Evaluated on [I3,I2,exp(i delta)]. Equal positive coefficients are a
    declared diagonal metric branch. No restriction to a standalone U1 theory.
    """
    # Keep Gaussian exponents near 150 to avoid underflow in a claimed tail bound.
    cutoff3=math.ceil(math.sqrt(450/tau))
    cutoff2=math.ceil(math.sqrt(600/tau))
    cutoff1=math.ceil(math.sqrt(150/tau))
    a3=np.zeros(3);a2=np.zeros(2)
    for p in range(cutoff3+1):
        for q in range(cutoff3+1):
            dim=(p+1)*(q+1)*(p+q+2)/2
            cas=(p*p+q*q+p*q+3*p+3*q)/3
            a3[(p+2*q)%3]+=dim*dim*math.exp(-tau*cas)
    for ell in range(cutoff2+1):
        a2[ell%2]+=(ell+1)**2*math.exp(-tau*ell*(ell+2)/4)
    residue=np.zeros(6)
    for triality in range(3):
        for parity in range(2):
            residue[(-2*triality-3*parity)%6]+=a3[triality]*a2[parity]
    charges=np.arange(-cutoff1,cutoff1+1)
    weights=np.exp(-tau*charges**2)*residue[charges%6]
    value=float(np.dot(weights,np.cos(charges*delta)))
    diagonal=float(sum(weights))
    s4=series_bound(tau/3,4,cutoff3);s2=series_bound(tau/3,2,cutoff3)
    b3=s4*s2;b2=series_bound(tau/4,2,cutoff2)
    b1=1+2*sum(math.exp(-tau*j*j) for j in range(1,cutoff1+1))+2*tail(tau,0,cutoff1)
    t3=tail(tau/3,4,cutoff3)*s2+tail(tau/3,2,cutoff3)*s4
    t2=tail(tau/4,2,cutoff2);t1=2*tail(tau,0,cutoff1)
    omitted=t3*b2*b1+b3*t2*b1+b3*b2*t1
    # Infinite-series tail only; IEEE rounding is separately tested, not certified.
    ratio_error=2*omitted/(diagonal-omitted)
    assert 0<omitted<diagonal and 0<ratio_error<1e-35
    return value/diagonal,dict(tau=tau,delta=delta,ratio=value/diagonal,
        absolute_series_tail_bound=omitted,ratio_series_tail_bound=ratio_error,
        cutoffs=[cutoff3,cutoff2,cutoff1],rounding_error_not_interval_certified=True)


def run():
    full=[prior.rep(*prior.group(67281,s)) for s in (0.,.1,.3)]
    _,nonabelian=family(full)
    angles=np.array([0.,.04,.12])
    links=[prior.rep(np.eye(3),np.eye(2),np.exp(1j*t)) for t in angles]
    k,abelian=family(links)
    r=float(abs(k[0,2]));necessary_log=math.log(r)
    rows=[]
    for tau in (.04,.25):
        b=np.ones((3,3));samples=[]
        for i in range(3):
            for j in range(i+1,3):
                h,meta=electric_kernel(tau,float(angles[j]-angles[i]))
                # Two spatial links, and two AP-time interfaces; bosons periodic.
                b[i,j]=b[j,i]=h**4;samples.append(meta)
        combined=b*k
        row=dict(tau=tau,electric_kernel=b.tolist(),electric_eigenvalues=np.linalg.eigvalsh(b).tolist(),
            combined=gram_checks(combined),pair02_minor=float((1-abs(combined[0,2])**2)),
            pair02_electric_log_cost=float(-math.log(b[0,2])),heat_samples=samples)
        assert min(row['electric_eigenvalues'])>-1e-12
        rows.append(row)
    assert min(rows[0]['combined']['eigenvalues'])>0
    assert min(rows[1]['combined']['eigenvalues'])<-.01 and rows[1]['pair02_minor']<0
    deps=('joint_nonflat_mass_measure.py','joint_static_gauge_reflection.py','research_note_643.md',
          'research_note_670.md','research_note_671.md','round672_drafts/dynamic_commutator_probe_results.json')
    return dict(date='2026-10-02',round=672,tests_run=2,failures=0,errors=0,
        original_SM_half_history=nonabelian,original_group_electric_sewing=dict(
            angles=angles.tolist(),mass_scale=.37,nt=2,nx=2,auxiliary_E='e0 at each site',
            actual_16_channels=abelian,required_pair02_log_cost=necessary_log,
            pair02_maximum_electric_coherence=1/r,rows=rows),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope='Original signed finite matter weight on independent gauge half histories, fixed shared temporal links and original scalar/auxiliary data. A bare kernel has a robust numerical negative witness; analytic Cauchy-Schwarz criterion constrains any multiplicative bosonic sewing. Original quotient-group electric heat factor in a declared diagonal-metric two-step candidate can pass or fail the sampled scalar test. Not the full Gauss/boundary-integrated process, not a certified floating-point proof, not all-polynomial RP or a no-go for the original Hamiltonian/continuum/GR.',
        all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=672,tests_run=2,all_checks_passed=True,
        heat_minima=[x['combined']['eigenvalues'][0] for x in result['original_group_electric_sewing']['rows']])))
