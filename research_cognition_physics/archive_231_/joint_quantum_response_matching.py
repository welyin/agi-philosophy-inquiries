"""602: match state, static susceptibility, causal response and source noise.

Exact frozen-background quadratic CAR sector of the original 598 masses.
No claim of a full Gauss-projected interacting state or continuum limit.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as matter
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_quantum_response_matching_results.json'
_,VAC,_=matter.original.lattice.scalar.parameters()
HIGGS=float(np.sqrt(VAC[0])); S0=float(np.sqrt(VAC[1]))
TOL=1e-10

def matrices(s,indices=None):
    phi=np.array([0.,HIGGS,0.,0.,s])
    h,d=matter.mass_matrices(phi);f=float(matter.original.F(phi))
    hnum=h*np.sqrt(f);dnum=d*np.sqrt(f)
    ds=np.zeros_like(dnum)
    ds[matter.SLICES['nu'],matter.SLICES['nu']]=matter.Y['s']*np.array([[0,1],[-1,0]])
    if indices is not None:
        ix=np.ix_(indices,indices)
        hnum,dnum,ds=hnum[ix],dnum[ix],ds[ix]
    u=f**-.5;up=s/(6*f**1.5);upp=1/(6*f**1.5)+s*s/(12*f**2.5)
    pairs=[(hnum*u,dnum*u),(hnum*up,ds*u+dnum*up),
           (hnum*upp,2*ds*up+dnum*upp)]
    def bdg(h,d):return np.block([[h,d],[d.conj().T,-h.T]])
    return [bdg(*pair) for pair in pairs],pairs

def thermal(H,beta):
    e,v=np.linalg.eigh(H)
    p=np.exp(-beta*(e-e.min()));p/=p.sum()
    return e,v,p,(v*p)@v.conj().T

def response(s=S0,n=1.,beta=2.,indices=None):
    (B,J,C),_=matrices(s,indices)
    H=n*B;gs=[B,n*J];cs=[[np.zeros_like(B),J],[J,n*C]]
    e,v=np.linalg.eigh(H);f=1/(1+np.exp(beta*e))
    gg=np.array([v.conj().T@g@v for g in gs])
    de=e[None,:]-e[:,None];df=f[:,None]-f[None,:]
    same=abs(de)<TOL;weight=np.zeros_like(de)
    np.divide(df,de,out=weight,where=~same)
    pop=(f*(1-f))[:,None]*same
    mean=np.array([.5*np.sum(f*np.diag(g)).real for g in gg])
    contact=np.array([[.5*np.sum(f*np.diag(v.conj().T@c@v)).real for c in row] for row in cs])
    ch=np.zeros((2,2));D=np.zeros_like(ch)
    for a in range(2):
        for b in range(2):
            products=(gg[a]*gg[b].T).real
            ch[a,b]=.5*np.sum(weight*products)
            D[a,b]=.5*np.sum(pop*products)
    free=-.5/beta*np.sum(np.logaddexp(0,-beta*e))
    return dict(mean=mean,contact=contact,chi_R=ch,D=D,chi_E=ch+beta*D,
                K_E=contact-ch-beta*D,K_R=contact-ch,free=float(free),
                e=e,f=f,g=gg,de=de,weight=weight,same=same)

def source_matching_check():
    r=response();fd=np.empty((2,2));eps=2e-5
    for b in range(2):
        args1=dict(s=S0,n=1.);args2=dict(args1)
        key=('n','s')[b];args1[key]+=eps;args2[key]-=eps
        fd[:,b]=(response(**args1)['mean']-response(**args2)['mean'])/(2*eps)
    error=float(np.max(abs(fd-r['K_E'])))
    assert error<2e-8
    (B,J,C),_=matrices(S0)
    bp=matrices(S0+eps)[0][0];bm=matrices(S0-eps)[0][0]
    derivative_error=float(np.max(abs((bp-bm)/(2*eps)-J)))
    second_error=float(np.max(abs((bp-2*B+bm)/eps**2-C)))
    assert derivative_error<1e-9 and second_error<2e-6
    assert np.min(np.linalg.eigvalsh(r['D']))>1e-5
    assert abs(r['chi_R'][0,0])<1e-24 and r['chi_R'][1,1]>1e-5
    return dict(frozen_phi=[0.,HIGGS,0.,0.,S0],F=float(matter.original.F(np.array([0,HIGGS,0,0,S0]))),
        CAR_modes=32,BdG_dimension=64,beta=2.,mean_sources=r['mean'].tolist(),
        static_free_energy_Hessian=r['K_E'].tolist(),isolated_zero_frequency_kernel=r['K_R'].tolist(),
        conserved_noise_matrix=r['D'].tolist(),noise_eigenvalues=np.linalg.eigvalsh(r['D']).tolist(),
        susceptibility_E=r['chi_E'].tolist(),susceptibility_R_zero=r['chi_R'].tolist(),
        Hessian_finite_difference_error=error,mass_derivative_error=derivative_error,
        mass_second_derivative_error=second_error)

def annihilators(m):
    out=[]
    for i in range(m):
        a=np.zeros((2**m,2**m),complex)
        for n in range(2**m):
            if n&(1<<i):a[n^(1<<i),n]=(-1)**((n&((1<<i)-1)).bit_count())
        out.append(a)
    return out

def fock(pair):
    h,d=pair;aa=annihilators(len(h));B=np.zeros_like(aa[0])
    for i in range(len(h)):
        for j in range(len(h)):
            B+=h[i,j]*aa[i].conj().T@aa[j]
            term=.5*d[i,j]*aa[i].conj().T@aa[j].conj().T
            B+=term+term.conj().T
    assert np.max(abs(B-B.conj().T))<1e-13
    return B

def fock_check():
    ids=[24,25,30,31]
    _,pairs=matrices(S0,ids)
    B,J,C=[fock(p) for p in pairs]
    e,v,p,rho=thermal(B,2.)
    g=v.conj().T@J@v;mean=float(np.trace(rho@J).real)
    dc=g-mean*np.eye(len(e));de=e[None,:]-e[:,None];same=abs(de)<TOL
    weight=np.zeros_like(de);np.divide(p[:,None]-p[None,:],de,out=weight,where=~same)
    chi=float(np.sum(weight*abs(g)**2))
    D=float(np.sum(p[:,None]*same*abs(dc)**2))
    q=response(indices=ids)
    errors=[abs(mean-q['mean'][1]),abs(chi-q['chi_R'][1,1]),abs(D-q['D'][1,1])]
    times=[0.,.5,2.,7.];noise=[];step=[]
    for t in times:
        direct_noise=float(np.sum(p[:,None]*abs(dc)**2*np.cos(de*t)))
        bdg_noise=float(.5*np.sum(q['f'][:,None]*(1-q['f'][None,:])*abs(q['g'][1])**2*np.cos(q['de']*t)))
        errors.append(abs(direct_noise-bdg_noise))
        noise.append(dict(t=t,exact_Fock_noise=direct_noise,BdG_noise=bdg_noise))
        if not t:continue
        contact=float(np.trace(rho@C).real)
        expected=contact-float(np.sum(weight*abs(g)**2*(1-np.cos(de*t))))
        outputs=[];epsilon=1e-5
        for change in (epsilon,-epsilon):
            _,perturbed=matrices(S0+change,ids)
            H1,J1,_=[fock(p1) for p1 in perturbed]
            ee,vv=np.linalg.eigh(H1);U=(vv*np.exp(-1j*ee*t))@vv.conj().T
            outputs.append(float(np.trace(U@rho@U.conj().T@J1).real))
        got=(outputs[0]-outputs[1])/(2*epsilon)
        assert abs(got-expected)<2e-8
        step.append(dict(t=t,actual_unitary_step_derivative=got,causal_prediction=expected,error=abs(got-expected)))
    assert max(errors)<1e-13
    return dict(exact_Fock_dimension=16,retained_neutral_CAR_modes=ids,
                maximum_BdG_Fock_error=max(errors),noise_rows=noise,unitary_step_rows=step,
                initial_thermal_state_held_fixed_during_quench=True)

def state_and_spectrum_check():
    rows=[]
    for beta in (.7,2.,5.):
        q=response(beta=beta)
        # Nonzero-frequency FDT weights. Degenerate blocks are kept separately.
        fn=q['f'][:,None];fm=q['f'][None,:]
        active=(q['de']>TOL)&(abs(q['g'][1])>1e-8)
        lhs=fn*(1-fm)+fm*(1-fn)
        rhs=np.zeros_like(lhs)
        rhs[active]=(fn-fm)[active]/np.tanh(beta*q['de'][active]/2)
        fdt_error=float(np.max(abs(lhs-rhs)[active]))
        assert fdt_error<1e-13
        assert np.min(np.linalg.eigvalsh(q['D']))>=-1e-12
        # Keep proper beta fixed when rescaling a uniform lapse.
        n=1.23;aligned=response(n=n,beta=beta/n)
        source_error=abs(aligned['mean'][0]-q['mean'][0])
        wrong=response(n=n,beta=beta)
        wrong_change=wrong['mean'][0]-q['mean'][0]
        assert source_error<1e-13 and abs(wrong_change)>1e-3
        rows.append(dict(beta=beta,static_causal_gap=(beta*q['D']).tolist(),
            nonzero_frequency_FDT_weight_error=fdt_error,
            lapse_rescaled_to=n,proper_temperature_aligned_source_error=float(source_error),
            fixed_coordinate_beta_energy_change=float(wrong_change)))
    return dict(rows=rows,zero_frequency_noise_is_explicit_not_discarded=True,
                equilibrium_preparation_not_derived=True,thermalization_not_assumed=True)

def run():
    a=source_matching_check();b=fock_check();c=state_and_spectrum_check()
    deps=('research_note_529.md','research_note_558.md','research_note_586.md','research_note_590.md',
          'research_note_591.md','research_note_598.md','research_note_600.md','research_note_601.md',
          'joint_fermion_gauss_completion.py')
    return dict(round=602,tests_run=3,failures=0,errors=0,source_matching=a,Fock_validation=b,
        state_and_spectrum=c,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_all_32_mode_mass_sector=True,static_bosonic_background_declared=True,
            full_interacting_Gauss_state_not_computed=True,source_state_and_noise_jointly_matched=True,
            established_Lehmann_identity_applied_not_new_fundamental_theorem=True,
            no_quantum_continuum_or_GR_completion=True))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=602,tests=3,all_passed=True,noise_eigenvalues=result['source_matching']['noise_eigenvalues'])))
