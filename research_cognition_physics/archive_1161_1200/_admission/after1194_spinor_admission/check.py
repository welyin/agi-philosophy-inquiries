"""Spinor-frame fixed-input kernel diagnostic; not a quantum-gravity implementation."""
from pathlib import Path
import argparse, json
import numpy as np

def run():
    I=np.eye(2,dtype=complex)
    X=np.array([[0,1],[1,0]],complex)
    Y=np.array([[0,-1j],[1j,0]],complex)
    Z=np.diag([1,-1]).astype(complex)
    gamma=[np.kron(Z,I)]+[1j*np.kron(Y,s) for s in (X,Y,Z)]
    alpha=[gamma[0]@g for g in gamma[1:]]
    g5=1j*gamma[0]@gamma[1]@gamma[2]@gamma[3]
    eta=np.diag([1.,-1.,-1.,-1.])
    mass=p=1.;k=.5;b=-k/2;T=.5
    L=8*np.pi;n=64;z=np.arange(n)*L/n
    generator=gamma[1]@gamma[2]
    S=np.cos(k*z/2)[:,None,None]*np.eye(4)+np.sin(k*z/2)[:,None,None]*generator
    phi=np.array([1.,1.,0.,0.],complex)/np.sqrt(2)
    psi=np.einsum('zij,j->zi',S,phi)
    hrot=p*alpha[0]+mass*gamma[0]+b*g5
    def evolve(h,v,t):
        e,u=np.linalg.eigh(h)
        return u@(np.exp(-1j*t*e)*(u.conj().T@v))
    twist=np.einsum('zij,j->zi',S,evolve(hrot,phi,T))
    coeff=np.fft.fft(psi,axis=0,norm='ortho')
    momenta=2*np.pi*np.fft.fftfreq(n,d=L/n)
    propagated=np.empty_like(coeff)
    h2psi=np.empty_like(coeff)
    for j,pz in enumerate(momenta):
        h=p*alpha[0]+pz*alpha[2]+mass*gamma[0]
        propagated[j]=evolve(h,coeff[j],T)
        h2psi[j]=h@h@coeff[j]
    flat=np.fft.ifft(propagated,axis=0,norm='ortho')
    Etheta2=mass**2+(p+b)**2
    Eflat2=mass**2+p**2+b**2
    overlap=lambda v:abs(np.vdot(psi,v)/n)**2
    Pt=1-(p+b)**2/Etheta2*np.sin(np.sqrt(Etheta2)*T)**2
    P0=1-(p*p+b*b)/Eflat2*np.sin(np.sqrt(Eflat2)*T)**2
    # Only modes +/- k/2 are present; this is exact Fourier support, not a UV estimate.
    keep=np.isclose(abs(momenta),k/2)
    support_error=float(np.linalg.norm(coeff[~keep])/np.sqrt(n))
    beta_mean=np.vdot(psi,np.einsum('ij,zj->zi',gamma[0],psi)).real/n
    result={
        'status':'PASS_DIAGNOSTIC_ONLY',
        'parameters':{'m':mass,'p_x':p,'k':k,'b':b,'T':T,'L':'8*pi','kL':'4*pi'},
        'same_original_input':'psi(z,x)=S(z) exp(i p_x x) phi; beta phi=phi, Sigma_x phi=phi',
        'spin_structure_preserved':True,
        'exact_probabilities':{'twisted':'1-9/25*sin(5/8)^2','flat':'1-17/33*sin(sqrt(33)/8)^2'},
        'probabilities':{'twisted':float(Pt),'flat':float(P0),'difference':float(Pt-P0)},
        'verification':{
            'clifford_residual':float(max(np.linalg.norm(gamma[i]@gamma[j]+gamma[j]@gamma[i]-2*eta[i,j]*np.eye(4)) for i in range(4) for j in range(4))),
            'axial_identity_residual':float(np.linalg.norm(alpha[2]@np.kron(I,Z)-g5)),
            'hrot_hermitian_residual':float(np.linalg.norm(hrot-hrot.conj().T)),
            'fourier_support_error':support_error,
            'flat_square_identity_residual':float(np.linalg.norm(np.fft.ifft(h2psi,axis=0,norm='ortho')-Eflat2*psi)/np.sqrt(n)),
            'twisted_probability_residual':float(abs(overlap(twist)-Pt)),
            'flat_probability_residual':float(abs(overlap(flat)-P0)),
            'beta_mean':float(beta_mean)},
        'scope':'fixed-background effective Dirac kernel, known upper-spinor input; NOT a positive-energy Fock preparation or actual detector',
        'not_proved':['composite-frame positive real-time parent','same complete stress and backreaction','controlled loop error','quantum gravity or dark sector']}
    assert max(v for key,v in result['verification'].items() if key!='beta_mean')<1e-12
    assert Pt>P0
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--write-results',action='store_true')
    args=ap.parse_args()
    data=run()
    if args.write_results:
        with Path(__file__).with_name('results.json').open('x',encoding='utf-8') as f:
            json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps(data,ensure_ascii=False,indent=2))
