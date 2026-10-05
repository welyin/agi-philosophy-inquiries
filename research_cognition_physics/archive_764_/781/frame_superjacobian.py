"""781: explicit finite-point super-Jacobian of the 780 frame dictionary.

No continuum determinant, anomaly coefficient, or measure is inferred from
these finite matrices. Two exact jet checks track the surviving density and
the already-known 769 Wick contact in this specific frame map.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement
import argparse
import importlib.util
import json
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frame780',HERE.parent/'780/frame_quartet_dictionary.py')
frame=importlib.util.module_from_spec(spec)
spec.loader.exec_module(frame)
ETA,I=frame.ETA,frame.I4
PAIRS=list(combinations(range(4),2))
SYMM=list(combinations_with_replacement(range(4),2))


def bases():
    metric,lorentz=[],[]
    for a,b in SYMM:
        value=np.zeros((4,4))
        value[a,b]=value[b,a]=1
        metric.append(value)
    for a,b in PAIRS:
        value=np.zeros((4,4))
        value[a,b],value[b,a]=1,-1
        lorentz.append(ETA@value)
    return metric,lorentz


def point_jacobian(g,q):
    metric,lorentz=bases()
    e=frame.section(g)
    lam,_=frame.exp_derivative(q)
    ilam=np.linalg.inv(lam)
    sylvester=np.kron(I,e)+np.kron(e.T,I)
    columns=[]
    for h in metric:
        db=np.linalg.solve(sylvester,(ETA@h).reshape(-1,order='F')).reshape((4,4),order='F')
        columns.append((lam@db).reshape(-1))
    jr=[]
    for generator in lorentz:
        _,dl=frame.exp_derivative(q,generator)
        columns.append((dl@e).reshape(-1))
        lower=ETA@dl@ilam
        jr.append([lower[a,b] for a,b in PAIRS])
    return np.column_stack(columns),np.column_stack(jr)


def complete_point_factor():
    rng=np.random.default_rng(78101)
    baseline,_=point_jacobian(ETA,np.zeros((4,4)))
    constant=abs(np.linalg.det(baseline))
    assert abs(constant-1/16)<1e-14
    errors=dict(bosonic_factor=0.,super_factor=0.,odd_cotangent_square=0.,spin_determinant=0.)
    omitted=[]
    gamma=frame.gamma_matrices()
    for _ in range(32):
        h=rng.normal(size=(4,4))*.025
        g=ETA+h+h.T
        q=frame.lorentz_generator(rng,.12)
        je,jr=point_jacobian(g,q)
        db,dg=abs(np.linalg.det(je)),np.linalg.det(jr)
        predicted=constant/np.sqrt(abs(np.linalg.det(g)))
        errors['bosonic_factor']=max(errors['bosonic_factor'],abs(db-predicted*dg))
        errors['super_factor']=max(errors['super_factor'],abs(db/dg-predicted))
        # Full odd cotangent lift: even (E,lambda*) and odd (lambda,E*).
        zero=np.zeros((16,6))
        even=np.block([[je,zero],[zero.T,np.linalg.inv(jr).T]])
        odd=np.block([[jr,zero.T],[zero,np.linalg.inv(je).T]])
        full=abs(np.linalg.det(even)/np.linalg.det(odd))
        errors['odd_cotangent_square']=max(errors['odd_cotangent_square'],abs(full-(db/dg)**2))
        lower=ETA@q
        spin_generator=sum(lower[a,b]*gamma[a]@gamma[b]/4 for a in range(4) for b in range(4))
        spin=frame.exp_derivative(spin_generator)[0]
        errors['spin_determinant']=max(errors['spin_determinant'],abs(np.linalg.det(spin)-1))
        omitted.append(abs(dg-1))
    errors={key:float(value) for key,value in errors.items()}
    assert max(errors.values())<2e-12,errors
    assert max(omitted)>.01
    return dict(samples=32,normalized_flat_component_constant=constant,
                maximum_residuals=errors,omit_ghost_maximum_relative_factor_error=max(omitted),
                surviving_relative_berezinian='sqrt(abs(det(g0))/abs(det(g)))',
                lorentz_coordinate_density_cancels=True,
                pointwise_spin_determinant_one_is_not_anomaly_proof=True)


def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),Q(0)) for j in range(len(b[0]))] for i in range(len(a))]


def tr(a): return sum(a[i][i] for i in range(len(a)))
def eye(n): return [[Q(i==j) for j in range(n)] for i in range(n)]
def show(a): return [[str(x) for x in row] for row in a]


def local_density_jet():
    # A=g0^-1 h has nontrivial off-diagonal mixing. Characteristic-polynomial
    # coefficients give an independent determinant calculation for trace-log.
    a=[[Q(-1,7),Q(1,5),0,0],[Q(-1,5),Q(1,3),0,0],
       [0,0,Q(2,9),Q(1,4)],[0,0,Q(1,4),Q(-1,6)]]
    # Faddeev--LeVerrier computes det(I+t A) exactly, independent of logdet.
    powers=[eye(4)]
    for n in range(1,7): powers.append(mm(powers[-1],a))
    coefficients=[Q(1)]
    for n in range(1,5):
        coefficients.append(sum(((-1)**(k-1)*coefficients[n-k]*tr(powers[k]) for k in range(1,n+1)),Q(0))/n)
    order=6
    detseries=coefficients+[Q(0)]*(order+1-len(coefficients))
    # Recursively solve J(t)^2 det(I+tA)=1, J(0)=1.
    jac=[Q(1)]+[Q(0)]*order
    for n in range(1,order+1):
        known=sum((jac[i]*jac[j]*detseries[k] for i in range(n+1) for j in range(n+1-i)
                   for k in [n-i-j] if i<n and j<n),Q(0))
        jac[n]=-known/2
    log=[Q(0)]+[Q((-1)**n,2*n)*tr(powers[n]) for n in range(1,order+1)]
    exp=[Q(1)]
    for n in range(1,order+1):
        exp.append(sum((k*log[k]*exp[n-k] for k in range(1,n+1)),Q(0))/n)
    assert exp==jac
    h1,h2=log[1],2*log[2]
    assert h1==-tr(a)/2 and h2==tr(mm(a,a))/2
    return dict(order=order,trace_log_equals_independent_determinant_recursion=True,
                relative_jacobian_coefficients=[str(x) for x in jac],
                log_density_first_derivative=str(h1),log_density_second_derivative=str(h2),
                functional_trace_coefficient_assigned=False)


def concrete_wick_contact():
    # Two metric modes, same positive finite covariance. This is a smooth
    # finite-mode/state-difference diagnostic, never a bare W(x,x) evaluation.
    eta=[[Q(-1),0,0,0],[0,Q(1),0,0],[0,0,Q(1),0],[0,0,0,Q(1)]]
    h1=[[Q(1,3),Q(1,5),0,0],[Q(1,5),Q(1,7),0,0],[0,0,Q(-1,4),0],[0,0,0,Q(1,6)]]
    h2=[[Q(1,8),0,Q(1,9),0],[0,Q(-1,5),0,0],[Q(1,9),0,Q(1,3),0],[0,0,0,Q(1,11)]]
    modes=[mm(eta,h1),mm(eta,h2)]
    covariance=[[Q(2),Q(1,3)],[Q(1,3),Q(3,2)]]
    quadratic=[[sum((covariance[a][b]*mm(modes[a],modes[b])[i][j] for a in range(2) for b in range(2)),Q(0))
                 for j in range(4)] for i in range(4)]
    contact=[[-x/8 for x in row] for row in quadratic]
    # From E^T eta E=g, the 2nd-order contraction must vanish identically.
    transpose=lambda a:list(map(list,zip(*a)))
    k1,k2=mm(transpose(contact),eta),mm(eta,contact)
    fluct=[[sum((covariance[a][b]*mm(mm(transpose(modes[a]),eta),modes[b])[i][j]/4
                for a in range(2) for b in range(2)),Q(0)) for j in range(4)] for i in range(4)]
    residual=[[k1[i][j]+k2[i][j]+fluct[i][j] for j in range(4)] for i in range(4)]
    assert not any(any(row) for row in residual)
    assert any(any(row) for row in fluct)
    # A scalar component source is enough to exhibit a nonzero response loss.
    source=contact[1][1]
    assert source
    return dict(metric_modes=2,covariance=show(covariance),
                frame_mean_contact=show(contact),
                metric_reconstruction_contact_identity_exact=True,
                omit_frame_mean_contact_metric_defect=show(fluct),
                omit_composite_source_mean_response_defect=str(source),
                absolute_hadamard_coincidence_used=False)


def run():
    return dict(round=781,groups={
        'complete_frame_ghost_berezinian':complete_point_factor(),
        'surviving_local_density_jet':local_density_jet(),
        'same_composite_source_wick_contact':concrete_wick_contact()},
        all_checks_passed=True,continuum_measure_or_anomaly_computed=False,
        original_full_N1_N2_proven=False,nonlinear_quantum_transport_complete=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    result=run()
    destination=HERE/'frame_superjacobian_results.json'
    if args.write: destination.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if args.check: assert result==json.loads(destination.read_text(encoding='utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
