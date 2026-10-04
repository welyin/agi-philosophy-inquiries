"""738: full h/s response at the original632 matched constant vacuum.

Keeps the mixed unequal-Majorana threshold, all charged species and half weights.
No physical gauge-mode count or nonlinear semiclassical existence is inferred.
"""
import argparse
import hashlib
import json
from functools import lru_cache
from pathlib import Path
import numpy as np
import joint_background_contact_matching as old
import joint_tensor_stress_spectrum as tensor
import joint_causal_log_inverse as causal
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_mixed_neutral_response_results.json'


def maxabs(a):return float(np.max(np.abs(a)))


@lru_cache(maxsize=1)
def original_data():
    F=old.matter.original.M-np.sum(old.U0)/6
    x=np.sqrt(old.U0/F)
    yn,ys=old.matter.Y['nu'],old.matter.Y['s']
    derivatives=np.array([[[0,yn],[yn,0]],[[0,0],[0,ys]]],complex)
    M=np.einsum('a,aij->ij',x,derivatives)
    beta=-np.angle(ys)/2;alpha=-np.angle(yn)-beta
    phase=np.diag(np.exp(1j*np.array([alpha,beta])))
    real=phase.T@M@phase
    eigen,O=np.linalg.eigh(real.real)
    U=phase@O@np.diag(np.where(eigen<0,1j,1.))
    m=np.abs(eigen)
    T=np.array([U.T@d@U for d in derivatives])
    assert maxabs(U.T@M@U-np.diag(m))<2e-15
    assert maxabs(np.einsum('a,aij->ij',x,T)-np.diag(m))<2e-15
    compressed,weights=old.canonical_mass_data(old.U0)
    assert maxabs(m-compressed[-2:])<2e-15
    return x,m,T,compressed,weights


@lru_cache(maxsize=1)
def terms():
    x,m,T,compressed,_=original_data()
    items=[]
    for k,(name,n) in enumerate((('u',3),('d',3),('e',1))):
        C=np.diag([n*abs(old.matter.Y[name])**2/(4*np.pi**2),0.])
        items.append((name,2*compressed[k],0.,'S',C))
    for i in range(2):
        for j in range(2):
            for kind,vector in (('S',T[:,i,j].real),('P',T[:,i,j].imag)):
                # Keep every term including floating-point near-zero entries.
                C=np.outer(vector,vector)/(8*np.pi**2)
                items.append((f'nu{i}{j}_{kind}',m[i]+m[j],abs(m[i]-m[j])/(m[i]+m[j]),kind,C))
    return items


def profile(x,r,kind):
    a=np.maximum(0.,1-x*x);b=np.maximum(0.,1-r*r*x*x)
    return np.sqrt(a*b)*(a if kind=='S' else b)


def spectrum(omega,off_diagonal=True):
    total=np.zeros((2,2))
    for label,a,r,kind,C in terms():
        if not off_diagonal and label.startswith(('nu01','nu10')):continue
        if omega>a:total+=np.pi*omega**2*C*profile(a/omega,r,kind)
    return total


def direct_cut(omega,direction):
    x,m,T,charged,_=original_data()
    alpha,beta=tensor.dirac_matrices()
    kinetic=np.einsum('a,aij->ij',direction,alpha)
    gamma5=np.block([[np.zeros((2,2)),np.eye(2)],[np.eye(2),np.zeros((2,2))]])
    answer=np.zeros((2,2))
    for k,(name,n) in enumerate((('u',3),('d',3),('e',1))):
        mass=charged[k]
        if omega<=2*mass:continue
        E=omega/2;p=np.sqrt(E*E-mass*mass)
        plus=(np.eye(4)+(p*kinetic+mass*beta)/E)/2;minus=np.eye(4)-plus
        vertex=abs(old.matter.Y[name])*beta
        answer[0,0]+=n*p*E/(2*np.pi)*np.trace(minus@vertex@plus@vertex).real
    for i in range(2):
        for j in range(2):
            if omega<=m[i]+m[j]:continue
            Ei=(omega*omega+m[i]*m[i]-m[j]*m[j])/(2*omega)
            Ej=omega-Ei;p=np.sqrt(max(0.,Ei*Ei-m[i]*m[i]))
            minus=(np.eye(4)-(p*kinetic+m[i]*beta)/Ei)/2
            plus=(np.eye(4)+(p*kinetic+m[j]*beta)/Ej)/2
            vertices=[beta*t.real+1j*beta@gamma5*t.imag for t in T[:,i,j]]
            for a in range(2):
                for b in range(2):
                    answer[a,b]+=.5*p*Ei*Ej/(np.pi*omega)*np.trace(minus@vertices[a]@plus@vertices[b]).real
    return answer


@lru_cache(maxsize=4)
def quad(n=160):
    z,w=np.polynomial.legendre.leggauss(n)
    theta=(z+1)*np.pi/4
    return np.sin(theta),w*np.pi/4*np.cos(theta)


def kernel(s,n=160):
    x,w=quad(n);result=np.zeros((2,2),complex)
    for _,a,r,kind,C in terms():
        result+=C*np.dot(w,s*s*x*profile(x,r,kind)/(a*a+s*s*x*x))
    return result


def decomposition(s,n=160):
    x,w=quad(n);Ctotal=np.zeros((2,2));B=np.zeros((2,2));R=np.zeros((2,2),complex)
    for _,a,r,kind,C in terms():
        gap=1-profile(x,r,kind)
        Ctotal+=C
        B+=C*(-np.log(a)-np.dot(w,gap/x))
        R+=C*(.5*np.log(1+(a/s)**2)+np.dot(w,gap*a*a/(x*(a*a+s*s*x*x))))
    return Ctotal,B,R


def mixed_spectrum_check():
    x,m,T,compressed,weights=original_data()
    omega=float((2*m[1]+m.sum())/2)
    crossing=spectrum(omega)-spectrum(omega,False)
    assert omega>m.sum() and omega<2*m[1]
    assert np.linalg.norm(crossing)>1e-5
    assert abs(x@crossing@x)<2e-16
    errors=[];projection=[]
    probes=(.25,omega,.9,1.6)
    for w in probes:
        analytic=spectrum(w)
        for n in (np.array([1.,0,0]),np.array([1.,2.,3.])/np.sqrt(14)):
            errors.append(maxabs(analytic-direct_cut(w,n)))
        previous=float(causal.scalar.density(np.array([w]),compressed,weights)[0])
        projection.append(abs(x@analytic@x-previous))
        assert min(np.linalg.eigvalsh(analytic))>-2e-15
    assert max(errors)<3e-14 and max(projection)<3e-15
    C,_,_=decomposition(1.2+.3j)
    Y=old.matter.Y
    expected=np.diag([(3*abs(Y['u'])**2+3*abs(Y['d'])**2+abs(Y['e'])**2+abs(Y['nu'])**2)/(4*np.pi**2),
                      abs(Y['s'])**2/(8*np.pi**2)])
    assert maxabs(C-expected)<2e-16 and np.linalg.det(C)>0
    return dict(original_mass_coordinates=x.tolist(),original_neutral_masses=m.tolist(),
                mixed_pair_threshold=float(m.sum()),probe_between_mixed_and_heavy_threshold=omega,
                missing_mixed_cut_matrix=crossing.tolist(),
                missing_mixed_cut_eigenvalues=np.linalg.eigvalsh(crossing).tolist(),
                common_scaling_cancels_mixed_pair=float(x@crossing@x),
                direct_unequal_mass_projector_error=max(errors),
                old_common_scaling_spectrum_projection_error=max(projection),
                full_log_coefficient=C.tolist(),full_log_eigenvalues=np.linalg.eigvalsh(C).tolist(),
                original_heat_kinetic_Gram_error=maxabs(C-expected),
                one_generation_and_Majorana_half_weights_retained=True)


def matrix_inverse_check():
    maximum=0.;quadrature=0.;positive=[]
    for s in (1.1+.3j,2.7+1.4j,.4+1.1j,3.2-.7j):
        C,B,R=decomposition(s)
        exact=kernel(s)
        maximum=max(maximum,maxabs(exact-C*np.log(s)-B-R))
        quadrature=max(quadrature,maxabs(exact-kernel(s,256)))
        positive.append(float(min(np.linalg.eigvalsh(exact.imag/(s*s).imag))))
    assert maximum<2e-13 and quadrature<2e-13 and min(positive)>0
    eig,V=np.linalg.eigh(C);W=(V*(1/np.sqrt(eig)))@V.T
    T=.01
    Bnorm=float(np.linalg.norm(W@B@W,2))
    remainder=0.
    for _,a,r,kind,Cj in terms():
        assert a*T<1
        remainder+=np.linalg.norm(W@Cj@W,2)*a*a*T*T*(7/4+np.log(1/(a*T)))
    upper=causal.inverse_norm_upper(T)
    rate=upper*(Bnorm+remainder)
    assert rate<1
    return dict(exact_matrix_decomposition_error=maximum,quadrature_refinement_error=quadrature,
                positive_measure_sample_minimum=min(positive),
                exact_constant_matrix=B.tolist(),normalized_constant_norm=Bnorm,
                time_window=T,normalized_memory_L1_upper=float(remainder),
                inverse_log_L1_upper=upper,sufficient_contraction_upper=float(rate),
                original_finite_mass_thresholds_retained=True,
                additional_local_gravity_terms_and_constraints_still_separate=True,
                short_time_inverse_not_a_nonlinear_spacetime=True)


def run():
    names=('mixed_spectrum_check','matrix_inverse_check')
    results={name:globals()[name]() for name in names}
    deps=('research_note_599.md','research_note_630.md','research_note_632.md',
          'research_note_730.md','research_note_736.md','research_note_737.md',
          'joint_fermion_gauss_completion.py','joint_background_contact_matching.py',
          'joint_tensor_stress_spectrum.py','joint_causal_log_inverse.py')
    return dict(round=738,tests_run=2,failures=0,errors=0,checks=list(names),results=results,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Full two-neutral-source free-fermion response on the original matched constant vacuum, including mixed Majorana thresholds, positive matrix logarithmic coefficient and exact short-time memory inverse. Full gauge constraints, interacting loops, actual time-dependent reference and nonlinear semiclassical existence remain open.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
