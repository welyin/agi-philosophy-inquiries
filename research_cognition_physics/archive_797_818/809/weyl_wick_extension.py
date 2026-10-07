"""809: noncommuting Weyl/polynomial Gaussian Gram and analytic-vector calibration.

No oscillator cutoff, new state for the original theory, or PDE computation.
The finite CCR covariance is diagnostic; the CAR marginal uses original matrices.
"""
from pathlib import Path
import argparse,json,math,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'weyl_wick_extension_results.json'
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original

def coefficient(ls,lt,qs,qt,qst,m,n):
    c=np.zeros((m+1,n+1),complex);c[0,0]=1
    for b in range(1,n+1):
        c[0,b]=(lt*c[0,b-1]+(2*qt*c[0,b-2] if b>=2 else 0))/b
    for a in range(1,m+1):
        for b in range(n+1):
            c[a,b]=(ls*c[a-1,b]+(2*qs*c[a-2,b] if a>=2 else 0)
                    +(qst*c[a-1,b-1] if b else 0))/a
    return c[m,n]*math.factorial(m)*math.factorial(n)/1j**(m+n)

def run():
    rng=np.random.default_rng(809);V=np.array([[1.2,.15],[.15,.9]])
    J=np.array([[0.,1.],[-1.,0.]])
    floor=float(np.linalg.eigvalsh(V+.5j*J).min());assert floor>0
    occ=float(original.occupied(np.array([0.,1.,0.]))[30,30].real)
    rho=np.diag([1-occ,occ])
    items=[]
    for i in range(12):
        items.append((rng.normal(size=2)*.6,rng.normal(size=2)*.5,i%4,
            rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))))
    def boson(i,j):
        ki,ai,m,_=items[i];kj,aj,n,_=items[j];d=kj-ki
        pref=np.exp(.5j*ki@J@kj-.5*d@V@d)
        return pref*coefficient(-ai@V@d-.5j*ai@J@d,
            -aj@V@d-.5j*d@J@aj,-.5*ai@V@ai,-.5*aj@V@aj,
            -ai@V@aj-.5j*ai@J@aj,m,n)
    gram=np.array([[boson(i,j)*np.trace(rho@items[i][3].conj().T@items[j][3])
                    for j in range(12)] for i in range(12)])
    herm=float(np.max(abs(gram-gram.conj().T)))
    ev=np.linalg.eigvalsh((gram+gram.conj().T)/2)
    assert herm<1e-12 and ev.min()>-1e-12
    positives=[]
    for _ in range(12):
        z=rng.normal(size=12)+1j*rng.normal(size=12)
        val=float(np.vdot(z,gram@z).real);positives.append(val);assert val>0
    # The imaginary covariance is essential for ordered CCR contractions.
    a=np.array([1.,0.]);b=np.array([0.,1.])
    ordered=a@(V+.5j*J)@b;reverse=b@(V+.5j*J)@a
    assert abs(ordered-reverse-1j)<1e-14
    phases=np.exp(-.5j*a@J@b),np.exp(.5j*a@J@b)
    phase_defect=float(abs(phases[0]-phases[1]));assert phase_defect>.5
    # Analytic-vector example: x^3 Omega, checked in Gaussian L2, no Fock cutoff.
    nodes,weights=np.polynomial.hermite.hermgauss(100)
    x=np.sqrt(2)*nodes;weights=weights/np.sqrt(np.pi);alpha=.7;rows=[]
    polynomial=np.ones_like(x,dtype=complex);term=polynomial.copy()
    for n in range(1,33):
        term*=1j*alpha*x/n;polynomial+=term
        if n in (4,8,16,32):
            error=float(np.sqrt(np.dot(weights,abs(x**3*(polynomial-np.exp(1j*alpha*x)))**2)))
            rows.append(dict(Taylor_degree=n,gaussian_vector_error=error))
    assert rows[-1]['gaussian_vector_error']<1e-11
    assert all(b['gaussian_vector_error']<a['gaussian_vector_error'] for a,b in zip(rows,rows[1:]))
    return dict(round=809,all_checks_passed=True,Gaussian_uncertainty_floor=floor,
        original_CAR_marginal_occupation=occ,extended_Gram_dimension=12,
        max_Weyl_polynomial_degree=3,Gram_hermiticity_residual=herm,
        Gram_min_eigenvalue=float(ev.min()),positive_sample_min=min(positives),
        ordered_CCR_residual=float(abs(ordered-reverse-1j)),
        omitted_Weyl_order_phase_defect=phase_defect,
        analytic_vector_calibration=rows,original_boson_covariance_used_numerically=False,
        original_continuous_signal_numerically_computed=False,
        scope='Finite mixed CCR/CAR calibration; original distributional extension and formal positivity are proved in the note.',
        new_numbered_test_groups=1)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
