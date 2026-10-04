"""Unfrozen 556 probe: radial reduction of the same single-cell quartic model.

Internal Higgs components are not spatial dimensions. This is a global SU(2)
singlet quantum-mechanical sector, not local electroweak gauge-field dynamics.
"""
import json
import hashlib
import math
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_matter_energy_moment_control as poly
TARGET=HERE/'radial_domain_probe_results.json'


def exponential_integral_small(x):
    # E1(x)=-gamma-log x-sum_{k>=1}(-x)^k/(k k!), used only x<.1.
    assert 0<x<.1
    term=-x; total=term
    for k in range(2,80):
        term*=(-x)/k
        total+=term/k
    return -0.5772156649015328606-math.log(x)-total


def radial_moment_cut(k,delta,omega):
    x=omega*delta*delta
    if k==-2:
        return omega**2*exponential_integral_small(x)
    n=k+2
    gamma=math.factorial(n-1)*math.exp(-x)*sum(x**j/math.factorial(j) for j in range(n))
    return omega**(-k)*gamma


def run():
    omega=1.3;a=.5;g=3*a/4
    L,C,u=poly.matter();d=5
    q=[poly.var(i,d) for i in range(d)]
    R2=poly.add(*(poly.mul(q[i],q[i]) for i in range(4)))
    s2=poly.mul(q[4],q[4])
    delta=[poly.add(R2,poly.const(-u[0],d)),poly.add(s2,poly.const(-u[1],d))]
    V=poly.scale(poly.quadratic(delta,L),.25)
    T=poly.add(poly.const(a*d*omega,d),poly.scale(poly.add(R2,s2),-a*omega**2))
    H=poly.add(T,V)
    E1=poly.real(poly.expect(H,1/(2*omega)))
    E2=poly.norm2(H,1/(2*omega))
    assert math.isfinite(E2) and E2>=E1**2
    rows=[];nodes,weights=np.polynomial.legendre.leggauss(180)
    maxres=0.
    for cutoff in (1e-1,1e-2,1e-3,1e-4,1e-5):
        rm={k:radial_moment_cut(k,cutoff,omega) for k in (-2,-1,0,1,2)}
        cent=g*g*rm[-2]
        flat=a*a*(9*rm[-2]/16-6*omega*rm[-1]+17.5*omega**2*rm[0]
                  -8*omega**3*rm[1]+omega**4*rm[2])
        cross=-2*a*g*(.75*rm[-2]-4*omega*rm[-1]+omega**2*rm[0])
        whole=a*a*(16*omega**2*rm[0]-8*omega**3*rm[1]+omega**4*rm[2])
        assert abs(flat+cent+cross-whole)<1e-11
        low=math.log(cutoff);high=math.log(10/math.sqrt(omega))
        t=(high-low)*nodes/2+(high+low)/2;r=np.exp(t)
        measure=2*omega**2*r**4*np.exp(-omega*r*r)
        integral=(high-low)/2*np.sum(weights*measure*(g/r**2)**2)
        maxres=max(maxres,abs(integral-cent))
        assert abs(integral-cent)<1e-9
        rows.append(dict(radial_cutoff=cutoff,centrifugal_piece_norm_squared=cent,
                         flat_radial_piece_norm_squared=flat,cross_term=cross,
                         full_radial_kinetic_norm_squared=whole))
    slope=9*a*a*omega*omega/8
    measured=(rows[-1]['centrifugal_piece_norm_squared']-rows[-2]['centrifugal_piece_norm_squared'])/math.log(10)
    assert abs(measured-slope)<1e-7
    assert abs(rows[-1]['full_radial_kinetic_norm_squared']-6*a*a*omega**2)<1e-7
    # Global SU(2) acts transitively on the complex-doublet norm spheres.
    rng=np.random.default_rng(556);maxinv=0.
    for _ in range(24):
        z=rng.normal(size=4);z/=np.linalg.norm(z)
        alpha=z[0]+1j*z[1];beta=z[2]+1j*z[3]
        U=np.array([[alpha,beta],[-np.conj(beta),np.conj(alpha)]])
        X=rng.normal(size=4);Z=np.array([X[0]+1j*X[1],X[2]+1j*X[3]])
        transformed=U@Z
        maxinv=max(maxinv,abs(np.vdot(transformed,transformed).real-X@X))
    assert maxinv<1e-12
    return dict(status='preliminary_556_not_a_completed_round',omega=omega,a=a,
                inherited_L=L.tolist(),C=C.tolist(),u=u.tolist(),
                full_Cartesian_E1=E1,full_Cartesian_E2=E2,
                radial_operator_extra_coefficient=g,
                predicted_log_divergence_coefficient=slope,measured_last_decade_slope=measured,
                exact_finite_full_radial_kinetic_second_moment=6*a*a*omega**2,
                cutoff_diagnostics=rows,independent_log_quadrature_residual=maxres,
                SU2_norm_invariance_residual=maxinv,
                frozen_dependency_hashes={name:hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest()
                    for name in ('joint_matter_energy_moment_control.py',
                                 'joint_singlet_common_mass_rg_results.json')},
                scope=dict(single_cell_global_SU2_singlet_only=True,
                           cutoff_only_diagnoses_norm_divergence_not_new_boundary_dynamics=True,
                           no_contradiction_to_round555_smooth_full_space_theorem=True,
                           complete_local_gauge_model_and_actual_radial_POVM_not_yet_constructed=True))


if __name__=='__main__':
    result=run()
    if TARGET.exists():
        assert json.loads(TARGET.read_text('utf8'))==result
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('status','full_Cartesian_E2',
                 'predicted_log_divergence_coefficient','independent_log_quadrature_residual')}))
