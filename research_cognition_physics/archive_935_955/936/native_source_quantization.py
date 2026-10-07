"""936: native joint Hamiltonian, sources and canonical force.

The explicit symbol is a calibration, not the full Einstein--SM Hamiltonian.
The old material-clock coefficients fix only its energy/source units.
"""
from pathlib import Path
import argparse,hashlib,json,math
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'native_source_quantization_results.json'
EPS=.6;RADIUS=2.;NMAX=33
I=np.eye(4,dtype=complex)
def annihilation(k):
    c=np.zeros((4,4),complex)
    for s in range(4):
        if (s>>k)&1:c[s^(1<<k),s]=(-1)**((s&((1<<k)-1)).bit_count())
    return c
C=[annihilation(k) for k in range(2)]
N1=C[0].conj().T@C[0];N2=C[1].conj().T@C[1]
HOP=C[0].conj().T@C[1]+C[1].conj().T@C[0]
def norm(a):return float(np.linalg.norm(a,2))
def coherent(q,p,n):
    a=(q+1j*p)/np.sqrt(2*EPS)
    c=np.empty((len(a),n),complex);c[:,0]=np.exp(-abs(a)**2/2)
    for k in range(1,n):c[:,k]=c[:,k-1]*a/np.sqrt(k)
    return c
def quadrature(nr,nt):
    x,w=np.polynomial.legendre.leggauss(nr)
    r=(x+1)*RADIUS/2;wr=w*RADIUS/2
    theta=2*np.pi*np.arange(nt)/nt
    q=(r[:,None]*np.cos(theta)).ravel();p=(r[:,None]*np.sin(theta)).ravel()
    weight=np.broadcast_to((wr*r)[:,None]/(EPS*nt),(nr,nt)).ravel().copy()
    return q,p,weight,coherent(q,p,NMAX)
def sylvester(S_eig,V,Y):
    y=V.conj().transpose(0,2,1)@Y@V
    sol=y/(S_eig[:,:,None]+S_eig[:,None,:])
    return V@sol@V.conj().transpose(0,2,1)
def symbol(q,p,j1=0.,j2=0.,derivatives=False):
    r2=q*q+p*p;u=r2/RADIUS**2
    chi=np.exp(-u/(1-u));chi_q=chi*(-2*q/RADIUS**2)/(1-u)**2
    geom=.06*r2/(1+r2)
    B=geom[:,None,None]*I+(.12+j1)*N1+(.16+j2)*N2+.05*np.tanh(q)[:,None,None]*HOP
    ev,V=np.linalg.eigh(I-B);assert ev.min()>.60
    se=np.sqrt(ev);S=(V*se[:,None,:])@V.conj().transpose(0,2,1)
    out={'h':chi[:,None,None]*(I-S)}
    if derivatives:
        bq=(.12*q/(1+r2)**2)[:,None,None]*I+(.05/np.cosh(q)**2)[:,None,None]*HOP
        s1=sylvester(se,V,-np.broadcast_to(N1,B.shape))
        s2=sylvester(se,V,-np.broadcast_to(N2,B.shape))
        sq=sylvester(se,V,-bq)
        s12=sylvester(se,V,-(s1@s2+s2@s1))
        out.update(h1=-chi[:,None,None]*s1,h12=-chi[:,None,None]*s12,
                   hq=chi_q[:,None,None]*(I-S)-chi[:,None,None]*sq)
    return out
def assemble(data,j1=0.,j2=0.,derivatives=False):
    q,p,w,c=data
    pair=(c[:,:,None]*c.conj()[:,None,:]).reshape(len(q),NMAX*NMAX).T
    out={}
    for key,value in symbol(q,p,j1,j2,derivatives).items():
        tmp=pair@(w[:,None]*value.reshape(len(q),16))
        mat=tmp.reshape(NMAX,NMAX,4,4).transpose(0,2,1,3).reshape(4*NMAX,4*NMAX)
        assert norm(mat-mat.conj().T)<1e-12
        out[key]=(mat+mat.conj().T)/2
    return out
def evolve(h,t):
    e,v=np.linalg.eigh(h);return (v*np.exp(-1j*t*e/EPS))@v.conj().T
def canonical(n):
    a=np.diag(np.sqrt(np.arange(1,n)),1)
    return np.kron(np.sqrt(EPS/2)*(a+a.T),I),np.kron(-1j*np.sqrt(EPS/2)*(a-a.T),I)
def tail_upper(n):
    mu=RADIUS**2/(2*EPS)
    # Poisson tail from n, bounded by its first term times a geometric series.
    assert mu<n+1
    return math.exp(-mu+n*math.log(mu)-math.lgamma(n+1))/(1-mu/(n+1))
def run():
    bg=json.loads((HERE.parent/'861/magnetic_reduced_hamiltonian_results.json').read_text('utf-8'))['original_859_background_clock_reduction'][1]
    alpha,v=bg['a'],bg['actual_clock_speed'];qstar=v*v/(4*alpha);hstar=v/(2*alpha)
    inherited=2*alpha/v**3*qstar*qstar/hstar
    assert abs(inherited-.25)<1e-15
    grids=[(40,64),(64,96)]
    operators=[assemble(quadrature(*g),derivatives=True) for g in grids]
    fine=operators[-1];quad={k:norm(fine[k]-operators[0][k]) for k in fine}
    assert max(quad.values())<2e-7
    data=quadrature(*grids[-1]);fd=[]
    for step in (.004,.002):
        hp=assemble(data,step,0)['h'];hm=assemble(data,-step,0)['h']
        mix=sum(a*b*assemble(data,a*step,b*step)['h'] for a in (-1,1) for b in (-1,1))/(4*step*step)
        fd.append(dict(step=step,first_error=norm((hp-hm)/(2*step)-fine['h1']),mixed_error=norm(mix-fine['h12'])))
    assert 3.8<fd[0]['first_error']/fd[1]['first_error']<4.2
    assert 3.8<fd[0]['mixed_error']/fd[1]['mixed_error']<4.2
    # Exact infinite-dimensional identity, checked after retaining the one
    # omitted oscillator level needed by the canonical momentum matrix.
    Qbig,Pbig=canonical(NMAX);Hbig=fine['h']
    rows=[]
    matter=np.array([0,1,1j,0],complex)/np.sqrt(2)
    finite_states=[]
    for n in (8,12,20,32):
        d=4*n;h=Hbig[:d,:d];hq=fine['hq'][:d,:d];Q,P=canonical(n)
        raw=1j*(h@P-P@h)/EPS+hq
        corrected=(1j*(Hbig@Pbig-Pbig@Hbig)/EPS)[:d,:d]+hq
        # Without the off-subspace boundary, finite CCR cannot be exact.
        assert norm(corrected)<2e-7
        c=coherent(np.array([.4]),np.array([.2]),n)[0];c/=np.linalg.norm(c)
        psi=np.kron(c,matter);unitary=evolve(h,4.);final=unitary@psi
        energy_error=abs(np.vdot(final,h@final)-np.vdot(psi,h@psi))
        assert energy_error<1e-13 and norm(unitary.conj().T@unitary-np.eye(d))<1e-12
        joint=final.reshape(n,4);rho_m=joint.conj().T@joint
        # Transpose convention has no effect on purity/eigenvalues.
        purity=float(np.trace(rho_m@rho_m).real)
        pchange=float((np.vdot(final,P@final)-np.vdot(psi,P@psi)).real)
        assert purity<.999999 and abs(pchange)>1e-5
        # Integral support area/(2*pi*EPS) = RADIUS^2/(2*EPS).
        # ||h|| <= 1 on support, giving a deliberately conservative bound.
        bound=2*(RADIUS**2/(2*EPS))*math.sqrt(tail_upper(n))
        finite_states.append((n,final))
        rows.append(dict(oscillator_states=n,dimension=d,projection_boundary_force_norm=norm(raw),
            boundary_corrected_force_error=norm(corrected),energy_conservation_error=float(energy_error),
            matter_purity=purity,canonical_momentum_change=pchange,
            analytic_generator_tail_bound=bound,unitary_error_bound_T4=4*bound/EPS,
            minimum_H_eigenvalue=float(np.linalg.eigvalsh(h).min())))
    reference=finite_states[-1][1]
    for row,(n,vec) in zip(rows,finite_states):
        padded=np.zeros_like(reference);padded[:len(vec)]=vec
        row['state_difference_from_32']=float(np.linalg.norm(padded-reference))
    assert rows[0]['projection_boundary_force_norm']>1e-6
    assert rows[-1]['analytic_generator_tail_bound']<1e-8
    assert rows[1]['state_difference_from_32']<rows[0]['state_difference_from_32']
    # Retain source Hessian from the same H and connect to a finite record.
    n=20;d=4*n;h=fine['h'][:d,:d];h1=fine['h1'][:d,:d]
    c=coherent(np.array([.4]),np.array([.2]),n)[0];c/=np.linalg.norm(c)
    psi=np.kron(c,matter);effect=np.kron(np.eye(n),N1)
    e,V=np.linalg.eigh(h);f=np.exp(-4j*e/EPS)
    diff=e[:,None]-e[None,:];divided=np.empty_like(diff,dtype=complex)
    np.divide(f[:,None]-f[None,:],diff,out=divided,where=abs(diff)>1e-10)
    near=abs(diff)<=1e-10
    divided[near]=(-4j/EPS*np.exp(-2j*(e[:,None]+e[None,:])/EPS))[near]
    du=V@(divided*(V.conj().T@h1@V))@V.conj().T;uu=(V*f)@V.conj().T
    response=float(2*np.vdot(uu@psi,effect@du@psi).real)
    def record(j):
        hh=assemble(data,j,0)['h'][:d,:d];z=evolve(hh,4.)@psi
        return float(np.vdot(z,effect@z).real)
    delta=.0005;difference=(record(delta)-record(-delta))/(2*delta)
    assert abs(response-difference)<2e-7 and abs(response)>1e-4
    paths=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/native_model_reuse_decision.md',
           HERE.parent/'research_note_861.md',HERE.parent/'861/magnetic_reduced_hamiltonian_results.json',
           HERE.parent/'research_note_873.md',HERE.parent/'research_note_935.md']
    return dict(round=936,date='2026-10-07',all_scientific_checks_passed=True,
        same_native_H_for_evolution_source_and_force=True,epsilon=EPS,phase_support_radius=RADIUS,
        old_clock_alpha=alpha,old_clock_speed=v,source_unit=qstar,energy_unit=hstar,
        original873_normalized_mixed_coefficient=inherited,
        quadrature_comparison=quad,source_derivative_rows=fd,finite_process_rows=rows,
        record_source_derivative=response,record_source_finite_difference=difference,
        quadrature_comparison_is_not_a_rigorous_integration_certificate=True,
        source_independent_chart_and_cutoffs_required=True,
        calibration_symbol_is_not_full_Einstein_SM=True,
        full_physical_symbol_and_matching_still_required=True,
        arbitrary_H_native_quantization_not_claimed_to_generate_GR=True,
        original_Q_or_E_rejected=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
