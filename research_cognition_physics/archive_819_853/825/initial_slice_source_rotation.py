"""825: local positive-mass mode orthogonal to transported original modes.

All original 64-component principal/mass matrices are retained. Propagation
here uses the existing constant auxiliary coefficients on a periodic grid;
the actual curved Cauchy argument is analytic in the note, not this grid.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original
sys.path.insert(0,str(HERE.parent/'819'))
import logical_input_source as packets
TARGET=HERE/'initial_slice_source_rotation_results.json'

def run():
    N=129;x=np.linspace(-np.pi,np.pi,N,endpoint=False);dx=2*np.pi/N
    k=np.fft.fftfreq(N)*N;C=original.CHARGE;Ns=original.vertex.NUMERATORS[4]
    s=np.zeros(64,complex)
    eig,vec=np.linalg.eigh(original.GAMMA[2][np.ix_([30,31],[30,31])]);s[[30,31]]=vec[:,-1]
    nu=np.linalg.norm(Ns@s);positive=(s+Ns@s/nu)/np.sqrt(2)
    assert np.linalg.norm(Ns@positive-nu*positive)<1e-13
    eigens=[np.linalg.eigh(original.hamiltonian((0.,0.,float(w)))) for w in k]
    def evolve(f,t):
        coeff=np.fft.fft(f,axis=0);out=np.empty_like(coeff)
        for i,(values,vectors) in enumerate(eigens):
            out[i]=vectors@(np.exp(-1j*t*values)*(vectors.conj().T@coeff[i]))
        return np.fft.ifft(out,axis=0)
    def charge(f):return np.einsum('ij,nj->ni',C,f.conj())
    def inner(f,g):return dx*np.vdot(f,g)
    star=[]
    for center in (-.9,.9):
        chi,_=packets.compact(x,center,.45);chi/=np.sqrt(dx*np.dot(chi,chi))
        star.append(chi[:,None]*np.exp(6j*x[:,None])*s)
    dt=.17;old=[evolve(f,-dt) for f in star]
    basis=np.stack([old[0],old[1],charge(old[0]),charge(old[1])],axis=-1)
    gram=lambda a,b:dx*np.einsum('nik,nij->kj',a.conj(),b)
    assert np.max(abs(gram(basis,basis)-np.eye(4)))<1e-12
    assert np.linalg.norm(evolve(charge(old[0]),dt)-charge(star[0]))<1e-11
    # Twelve local smooth profiles versus four exact complex constraints.
    chi,_=packets.compact(x,-.8,.6)
    pol=np.polynomial.legendre.legvander((x+.8)/.6,11)
    profiles=chi[:,None]*pol
    candidates=np.einsum('nk,i->nik',profiles,positive)
    constraints=gram(basis,candidates)
    _,_,vh=np.linalg.svd(constraints,full_matrices=True)
    coeff=vh.conj().T[:,-1]
    g=np.einsum('nik,k->ni',candidates,coeff)
    g/=np.sqrt(inner(g,g).real)
    residual=max(abs(inner(basis[:,:,i],g)) for i in range(4))
    assert residual<1e-12 and abs(inner(charge(g),g))<1e-12
    assert np.linalg.norm(np.einsum('ij,nj->ni',Ns,g)-nu*g)<1e-12
    assert np.max(abs(g[chi==0]))==0
    eta=np.zeros(N);a=abs(x+.8);eta[a<=.7]=1
    transition=(a>.7)&(a<1.1);eta[transition]=.5*(1+np.cos(np.pi*(a[transition]-.7)/.4))
    source=lambda f:inner(f,eta[:,None]*np.einsum('ij,nj->ni',Ns,f)).real
    qa=source(old[0]);qc=source(g)
    cross=inner(old[0],eta[:,None]*np.einsum('ij,nj->ni',Ns,g))
    assert abs(qc-nu)<1e-12 and qc>0
    # A self-adjoint C-odd finite coefficient diagnostic for the continuity
    # of a fixed first-order response, not the original curved process K_j.
    q=1+.3*np.cos(x)+.2*np.sin(2*x)
    momentum=lambda f:np.fft.ifft(k[:,None]*np.fft.fft(f,axis=0),axis=0)
    def response(f):
        spatial=.5*(q[:,None]*momentum(f)+momentum(q[:,None]*f))
        return -.05*np.einsum('ij,nj->ni',original.GAMMA[2],spatial)+.04*q[:,None]*np.einsum('ij,nj->ni',original.MASS,f)
    rows=[];source_error=0.;orth_error=0.
    for theta in (0.,.001,.01,.03,.06):
        f=np.cos(theta)*old[0]+np.sin(theta)*g
        new0=[f,old[1],charge(f),charge(old[1])]
        Q0=np.stack(new0,axis=-1)
        orth_error=max(orth_error,float(np.max(abs(gram(Q0,Q0)-np.eye(4)))))
        Qs=np.stack([evolve(v,dt) for v in new0],axis=-1)
        JQ=np.stack([response(Qs[:,:,i]) for i in range(4)],axis=-1)
        out=JQ-np.einsum('nik,kj->nij',Qs,gram(Qs,JQ))
        leakage=np.linalg.eigvalsh(gram(out,out))
        assert leakage.min()>1e-7
        actual=source(f)
        predicted=qa*np.cos(theta)**2+qc*np.sin(theta)**2+cross.real*np.sin(2*theta)
        source_error=max(source_error,abs(actual-predicted))
        assert abs(actual)>1e-9
        assert np.linalg.norm(response(charge(Qs[:,:,0]))+charge(response(Qs[:,:,0])))<1e-10
        rows.append(dict(theta=theta,original_initial_slice_bilinear=actual,
            fixed_response_Gram_min=float(leakage.min()),
            propagated_mode_norm=float(inner(Qs[:,:,0],Qs[:,:,0]).real)))
    assert source_error<1e-12 and orth_error<1e-12
    return dict(round=825,all_checks_passed=True,grid=N,
        original_64_masses_and_principal_matrices_used=True,
        local_positive_mode_bilinear=float(qc),local_positive_mode_Ns_eigenvalue=float(nu),
        orthogonality_to_transported_old_CAR_space_error=float(residual),
        all_rotated_CAR_errors=orth_error,source_trigonometric_formula_error=float(source_error),
        old_source_coefficient=float(qa),cross_source_real=float(cross.real),rows=rows,
        new_mode_support_preserved_in_construction=True,
        original_curved_propagation_numerically_computed=False,
        original_initial_slice_nonzero_source_and_noise_connection_is_analytic=True,
        separated_packet_XY_source_property_preserved=False,
        actual_autonomous_mode_rotation_proven=False,formal_test_groups_added=1)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
