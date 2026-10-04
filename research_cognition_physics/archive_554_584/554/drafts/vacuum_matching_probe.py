"""Unnumbered interface check: finite tadpoles, additive matching and Goldstone IR.

This is not a quantum vacuum/stability proof. All couplings are fixed-scale
MS proxies inherited from 544, in Landau gauge and on flat classical geometry.
"""
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ARCHIVE = HERE.parent
sys.path.insert(0, str(ARCHIVE))
import joint_singlet_common_mass_rg as rg

TARGET = HERE/'vacuum_matching_probe_results.json'


def potential_and_gradient(x, C, L, couplings, mu2=1.):
    X,Y = x; lh,p,ls = L[0,0],L[0,1],L[1,1]
    q,S,z,gy,gw = couplings
    M = np.array([[-C[0]+3*lh*X+p*Y,2*p*math.sqrt(X*Y)],
                  [2*p*math.sqrt(X*Y),-C[1]+p*X+3*ls*Y]])
    vals,vecs = np.linalg.eigh(M)
    dM = [np.array([[3*lh,p*math.sqrt(Y/X)],[p*math.sqrt(Y/X),p]]),
          np.array([[p,p*math.sqrt(X/Y)],[p*math.sqrt(X/Y),3*ls]])]
    radial_grad = np.array([[vecs[:,i]@dm@vecs[:,i] for dm in dM] for i in range(2)])
    masses = np.array([*vals,gw*X/4,(gy+gw)*X/4,q*X/2,z*Y+S*X/2])
    grads = np.vstack((radial_grad,[[gw/4,0],[(gy+gw)/4,0],[q/2,0],[S/2,z]]))
    n = np.array([1.,1.,6.,3.,-12.,-4.])
    k = np.array([1.5,1.5,5/6,5/6,1.5,1.5])
    assert np.all(masses > 0)
    log = np.log(masses/mu2)
    value = np.sum(n*masses**2*(log-k))/(64*math.pi**2)
    grad = (n*masses*(2*log-2*k+1))@grads/(64*math.pi**2)
    # Three Landau-gauge Goldstones. At G=0 the value and first derivative
    # have finite zero limits; no finite second derivative is asserted.
    G = -C[0]+lh*X+p*Y
    if abs(G) > 1e-12:
        value += 3*G*G*(math.log(abs(G)/mu2)-1.5)/(64*math.pi**2)
        grad += 3*G*(2*math.log(abs(G)/mu2)-2)*np.array([lh,p])/(64*math.pi**2)
    return float(value),grad,float(G),masses


def run():
    source = ARCHIVE/'joint_singlet_common_mass_rg_results.json'
    saved = json.loads(source.read_text('utf8')); rows=[]
    for index in (2,4,5):
        row = saved['examples'][index]; a = row['state']
        q,S,z,lh,p,ls,ch,cs = [a[k] for k in ('q','S','z','lambda_H','p','lambda_s','x','y')]
        C=.25*np.array([ch,cs]); L=np.array([[lh,p],[p,ls]])
        x=np.linalg.solve(L,C)
        gy,gw,gc=rg.gauge_squared(-row['u'],rg.XSTAR)
        couplings=(q,S,z,gy,gw)
        v,grad,G,masses=potential_and_gradient(x,C,L,couplings)
        assert abs(G)<1e-14
        # First-order formal expansion; no exact one-loop minimizer is assumed.
        dx=-2*np.linalg.solve(L,grad)
        residual=np.max(abs(L@dx/2+grad))
        assert residual<1e-15
        assert np.linalg.norm(grad)>1e-6  # A constant alone cannot remove tadpoles.
        fd=[]
        for i in range(2):
            e=np.zeros(2); e[i]=1e-5
            plus=potential_and_gradient(x+e,C,L,couplings)[0]
            minus=potential_and_gradient(x-e,C,L,couplings)[0]
            fd.append((plus-minus)/(2*e[i]))
        diff=np.max(abs(np.array(fd)-grad)); assert diff<1e-9
        # Second derivative of 3 G^2(log|G|-3/2)/(64 pi^2), for G linear in x.
        direction=np.array([lh,p]); norm=np.linalg.norm(direction)
        tangent=np.array([-p,lh])/norm; normal=direction/norm
        goldstone=[]
        for exponent in (2,4,8,16,32):
            Hess=3*(-exponent*math.log(10))*np.outer(direction,direction)/(32*math.pi**2)
            assert abs(tangent@Hess@tangent)<1e-14
            goldstone.append(dict(log10_abs_G_over_mu2=-exponent,
                normal_Hessian=float(normal@Hess@normal)))
        assert all(goldstone[j+1]['normal_Hessian']<goldstone[j]['normal_Hessian'] for j in range(4))
        rows.append(dict(source_index=index,q0=row['q0'],tree_x=x.tolist(),
            one_loop_value_at_mu2_one=v,finite_vacuum_constant_shift=-v,
            one_loop_x_gradient=grad.tolist(),formal_first_order_dx=dx.tolist(),
            relative_first_order_dx=(dx/x).tolist(),formal_tree_F_change=float(-np.sum(dx)/6),
            tadpole_linear_equation_residual=float(residual),
            independent_gradient_difference=float(diff),Goldstone_Hessian_probe=goldstone))
    return dict(status='interface_diagnostic_not_completed_research_round',tests_run=3,failures=0,
        checks=['off_shell_mass_gradient_and_independent_finite_difference',
                'formal_tadpole_shift_additive_vacuum_matching_is_not_parameter_reduction',
                'Goldstone_second_derivative_diverges_in_a_rank_one_direction'],
        assumptions=dict(flat_classical_metric=True,one_loop_Landau_MS=True,
            fixed_dimensionless_mu_squared=1,real_part_for_finite_difference_only=True,
            no_pole_mass_or_full_stability_claim=True,no_finite_curvature_matching_claim=True),
        examples=rows,dependency_hashes={str(source.relative_to(ARCHIVE)):
            hashlib.sha256(source.read_bytes()).hexdigest(),
            'joint_singlet_common_mass_rg.py':hashlib.sha256((ARCHIVE/'joint_singlet_common_mass_rg.py').read_bytes()).hexdigest()})


if __name__=='__main__':
    result=run()
    if TARGET.exists():
        assert json.loads(TARGET.read_text('utf8'))==result
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
