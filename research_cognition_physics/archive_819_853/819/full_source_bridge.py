"""819: full-space covariance support and curved-coefficient source diagnostics.
Finite Clifford identities hold for a correlated state including partner and
remote modes. The variable-coefficient Dirac calculation is a diagnostic, not
a solution of the original future Einstein--matter equations.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
import logical_input_source as early
bridge=early.bridge;original=early.original
TARGET=HERE/'full_source_bridge_results.json'

def full_covariance_support():
    units=bridge.prep.data()[6]
    a=bridge.prep.previous.old.car(6);F=bridge.record_frame(a);xi=a+[c.conj().T for c in a]
    bilinear=np.array([[c@d.conj().T for d in xi] for c in xi])
    cov=lambda state:np.einsum('ab,ijba->ij',state,bilinear)
    rng=np.random.default_rng(819)
    v=rng.normal(size=64)+1j*rng.normal(size=64)
    v[[i for i in range(64) if i.bit_count()%2]]=0;v/=np.linalg.norm(v)
    average=np.zeros((64,64),complex)
    for U in units:
        w=np.kron(np.eye(2),U)@v;average+=np.outer(w,w.conj())/256
    block=F.conj().T@average@F
    corner=np.einsum('aras->rs',block.reshape(2,32,2,32))
    Pbar=cov(average);S=[0,1,6,7];mask=np.ones((12,12),bool);mask[np.ix_(S,S)]=False
    results={};full={}
    for name,p in zip(('x','y','z'),bridge.PAULI[1:]):
        family=[]
        for sign in (1.,-1.):
            rho=(np.eye(2)+sign*p)/2
            state=F@np.kron(rho,corner)@F.conj().T
            delta=cov(state)-Pbar
            family.append(cov(state));full[name+str(sign)]=cov(state)
            assert np.max(abs(delta[mask]))<1e-12
        D=family[0]-family[1]
        results[name]=dict(outside_S_max=float(np.max(abs(D[mask]))),
            diagonal_max=float(np.max(abs(np.diag(D)))),
            difference_norm=float(np.linalg.norm(D)),
            block_real=D[np.ix_(S,S)].real.tolist(),
            block_imag=D[np.ix_(S,S)].imag.tolist())
    Dz=full['z1.0']-full['z-1.0']
    expected=np.zeros((12,12));expected[0,0]=1;expected[6,6]=-1
    assert np.max(abs(Dz-expected))<1e-12
    assert results['x']['diagonal_max']<1e-12 and results['y']['diagonal_max']<1e-12
    return dict(arbitrary_correlated_even_reference_used=True,
        twirl_and_reset_include_partner_and_remote_mode=True,
        logical_source_blocks=results,
        Z_rank_two_difference_error=float(np.max(abs(Dz-expected))))

def common_polarization():
    rng=np.random.default_rng(819)
    gamma=[g[np.ix_([30,31],[30,31])] for g in original.GAMMA]
    rows=[]
    for i in range(12):
        n=rng.normal(size=3);n/=np.linalg.norm(n)
        Gn=sum((n[k]*gamma[k] for k in range(3)),np.zeros((2,2),complex))
        k=rng.normal(size=(4,4));k=(k+k.T)/2
        M=np.dot(k[0,1:],n)*np.eye(2,dtype=complex)
        vec=(k[1:,1:]+k[0,0]*np.eye(3))@n
        M+=sum((vec[j]*gamma[j]/2 for j in range(3)),np.zeros((2,2),complex))
        e,v=np.linalg.eigh(Gn);j=max(range(2),key=lambda j:np.linalg.norm(M@v[:,j]))
        active=float(np.linalg.norm(M@v[:,j]));energy=float(np.vdot(v[:,j],Gn@v[:,j]).real)
        assert active>1e-5 and abs(abs(energy)-1)<1e-12
        rows.append(dict(active_shear_norm=active,normal_energy_principal=energy))
    return dict(cases=rows,scope='Both helicities cannot lie in the kernel of a nonzero shear matrix.')

def variable_coefficient_energy():
    x,w=np.polynomial.legendre.leggauss(256)
    chi,dchi=early.compact(x,0.,.8)
    norm=np.sqrt(np.dot(w,chi*chi));chi/=norm;dchi/=norm
    v=1+.2*np.sin(x);dv=.2*np.cos(x)
    lapse=1+.1*np.cos(x);dlapse=-.1*np.sin(x)
    speed=lapse*v;dspeed=dlapse*v+lapse*dv
    e,vec=np.linalg.eigh(original.GAMMA[2][np.ix_([30,31],[30,31])])
    s=np.zeros(64,complex);s[[30,31]]=vec[:,-1]
    slope=float(np.dot(w,speed*chi*chi))
    intercept=float(.17*np.dot(w,speed*np.cos(x)*chi*chi))
    rows=[]
    for k in (4.,8.,16.,32.):
        phase=np.exp(1j*(k*x+.17*np.sin(x)))
        f=chi[:,None]*phase[:,None]*s
        df=(dchi+1j*(k+.17*np.cos(x))*chi)[:,None]*phase[:,None]*s
        kinetic=-1j*(speed[:,None]*(df@original.GAMMA[2].T)+.5*dspeed[:,None]*(f@original.GAMMA[2].T))
        Hf=kinetic+lapse[:,None]*(f@original.MASS.T)
        energy=np.dot(w,np.einsum('ti,ti->t',f.conj(),Hf))
        predicted=k*slope+intercept
        assert abs(energy.real-predicted)<1e-11 and abs(energy.imag)<1e-11
        rows.append(dict(kappa=k,measured_energy=float(energy.real),principal_plus_lower_order=predicted))
    return dict(principal_coefficient=slope,fixed_lower_order_intercept=intercept,
        source_values=rows,original_mass_retained=True,
        original_coupled_curved_background_solved=False)

def run():
    return dict(round=819,all_checks_passed=True,
        full_covariance=full_covariance_support(),
        common_polarization=common_polarization(),
        energy=variable_coefficient_energy(),
        original_full_constraint_compensation_computed=False,
        geometric_change_alone_forced=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
