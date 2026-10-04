"""600: full-stress and Yukawa-force completion of the inherited EFT map.

Tensor checks are local algebra, not a solved semiclassical stress tensor.
Actual mass derivatives use the 598 model and the original 574 potential.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_scalar_effective_geometry as scalar_loop
import joint_fermion_gauss_completion as fermion
import joint_fermion_scalar_loop_matching as mass_loop

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_full_matter_operator_reduction_results.json'


def tensor_factorization():
    rng=np.random.default_rng(600);rows=[]
    for case in range(18):
        X=rng.normal(size=(4,5))*.25;B=X@X.T;MX=X.T@X
        A=rng.normal(size=(5,5))*.15;A=(A+A.T)/2
        U=float(original.node_potential(rng.normal(size=5)*.3))
        tau=rng.normal(size=(4,4))*.2;tau=(tau+tau.T)/2
        if case%3==0:tau-=np.trace(tau)*np.eye(4)/4
        t=float(np.trace(tau));S=float(np.trace(B));Q=float(np.trace(B@B))
        Ric=rng.normal(size=(4,4))*.3;Ric=(Ric+Ric.T)/2
        W=tau-t*np.eye(4)/2;C=B+U*np.eye(4)+W;r0=S+4*U-t
        if case%3==1:Ric=C.copy()
        R=float(np.trace(Ric));ric2=float(np.trace(Ric@Ric))
        riem2=1.31+2*ric2-R*R/3;euler=riem2-4*ric2+R*R
        full=scalar_loop.entry.invariant_coefficient(X,A,R)+5*((riem2-ric2)/180+R*R/72)
        oldmet=(euler/36-7*S*S/216+11*Q/108+U*S/18+U*U
                -2*U*np.trace(A)/3+np.trace(A@A)/2-np.trace(A@MX)/6)
        extra=(np.sum(B*tau)/6+np.sum(tau*tau)/12
               -t*S/18-U*t/2+t*np.trace(A)/6+t*t/24)
        D=(Ric+C)/12+np.eye(4)*((R+r0)/24-np.trace(A)/6-S/9)
        E=Ric-C;e=float(np.trace(E));H=-(E-e*np.eye(4)/2)/2
        T=-2*D+np.trace(D)*np.eye(4)
        metric_error=float(abs(full-oldmet-extra-np.sum(H*T)))
        # Complete scalar EOM, including direct extra-matter force j.
        grad=rng.normal(size=5)*.2;tension=rng.normal(size=5)*.2;j=rng.normal(size=5)*.2
        Ephi=grad-tension+j
        lapU=float(np.trace(A@MX)+grad@tension)
        red=oldmet+np.trace(A@MX)/6+grad@grad/6+extra+grad@j/6
        total_error=float(abs(full-red-np.sum(H*T)+grad@Ephi/6+lapU/6))
        assert max(metric_error,total_error)<3e-14
        if case%3==1:assert np.max(abs(H))<1e-14
        rows.append(dict(case=case,stress_trace=t,extra_density=float(extra),
                         tracefree_shortcut_missing=float(-t*S/18-U*t/2+t*np.trace(A)/6+t*t/24),
                         metric_error=metric_error,complete_euler_error=total_error))
    assert max(abs(r['tracefree_shortcut_missing']) for r in rows)>.01
    return dict(rows=rows,max_error=max(max(r['metric_error'],r['complete_euler_error']) for r in rows),
                local_algebra_only=True,tracefree_limit_recovers_582=True)


def potential_gradient(phi):
    L,u,_=original.lattice.scalar.parameters()
    d=np.array([phi[:4]@phi[:4]-u[0],phi[4]**2-u[1]])
    v=L@d;f=float(original.F(phi));numerator=float(d@L@d)
    return phi*np.r_[np.full(4,v[0]),v[1]]/f**2+numerator*phi/(6*f**3)


def mass_directional(phi,v):
    # Mass numerators are linear; evaluate them through the saved 598 builder.
    f=float(original.F(phi));fv=float(original.F(v))
    hv,dv=fermion.mass_matrices(v)
    h,d=fermion.mass_matrices(phi)
    # M(phi)=N(phi)/sqrt(F(phi)); directional derivative of 1/sqrt(F).
    factor=float(phi@v)/(6*f)
    return hv*np.sqrt(fv/f)+h*factor,dv*np.sqrt(fv/f)+d*factor


def actual_yukawa_force():
    phi=np.array([.8,.3,-.2,.25,.7])
    covgrad=potential_gradient(phi)
    step=2e-5;finitegrad=[]
    for a in range(5):
        e=np.eye(5)[a]*step
        finitegrad.append((original.node_potential(phi+e)-original.node_potential(phi-e))/(2*step))
    gradient_error=float(np.max(abs(covgrad-finitegrad)))
    assert gradient_error<2e-10
    v=original.inverse(phi)@covgrad/6
    assert original.F(v)>0
    h,d=fermion.mass_matrices(phi);dh,dd=mass_directional(phi,v)
    assert np.linalg.norm(dh)+np.linalg.norm(dd)>1e-4
    rows=[]
    for lam in (.04,.02,.01,.005):
        hh,delta=fermion.mass_matrices(phi+lam*v)
        remainder=float(np.linalg.norm(hh-h-lam*dh)+np.linalg.norm(delta-d-lam*dd))
        rows.append(dict(eft_parameter=lam,mass_remainder=remainder,
                         remainder_over_parameter_squared=remainder/lam**2))
    ratios=[rows[i]['mass_remainder']/rows[i+1]['mass_remainder'] for i in range(3)]
    assert all(3.95<r<4.05 for r in ratios)
    return dict(original_phi=phi.tolist(),original_potential=float(original.node_potential(phi)),
                original_gradient_error=gradient_error,scalar_shift=v.tolist(),
                induced_dirac_coefficient_norm=float(np.linalg.norm(dh)),
                induced_majorana_coefficient_norm=float(np.linalg.norm(dd)),
                rows=rows,remainder_ratios=ratios,
                tests_bilinear_coefficients_not_a_prepared_fermion_state=True)


def curvature_mass_contact():
    # Same diagnostic Yukawas and all-left mass matrix as 599.
    M0,_=mass_loop.physical_mass_direction()
    q=.37;M=np.sqrt(6)*np.sinh(q/np.sqrt(6))*M0
    f=-float(np.trace(M.conj().T@M).real)/6  # Half-weight Majorana sector.
    rng=np.random.default_rng(6003);errors=[];contacts=[]
    for _ in range(16):
        X=rng.normal(size=(4,5))*.2;B=X@X.T;S=float(np.trace(B));U=.09
        tau=rng.normal(size=(4,4))*.2;tau=(tau+tau.T)/2;t=float(np.trace(tau))
        Ric=rng.normal(size=(4,4))*.2;Ric=(Ric+Ric.T)/2;R=float(np.trace(Ric))
        E=Ric-B-U*np.eye(4)-tau+t*np.eye(4)/2
        H=-(E-np.trace(E)*np.eye(4)/2)/2
        direct=f*R
        reduced=f*(S+4*U-t)
        residual=float(abs(direct-reduced-2*f*np.trace(H)))
        assert residual<1e-15
        errors.append(residual);contacts.append(float(-f*t))
    assert max(abs(c) for c in contacts)>.01
    return dict(mass_curvature_coefficient=f,max_factorization_error=max(errors),
                mass_stress_contacts=contacts,
                factorization_not_composite_operator_renormalization=True)


def run():
    tensor=tensor_factorization();force=actual_yukawa_force();contact=curvature_mass_contact()
    names=('research_note_351.md','research_note_352.md','research_note_581.md','research_note_582.md',
           'research_note_585.md','research_note_588.md','research_note_598.md','research_note_599.md',
           'joint_curved_quantum_source.py','joint_scalar_effective_geometry.py',
           'joint_fermion_gauss_completion.py','joint_fermion_scalar_loop_matching.py')
    return dict(round=600,tests_run=3,failures=0,errors=0,tensor=tensor,force=force,contact=contact,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope=dict(full_leading_action_used_for_redefinition=True,
                           extra_stress_trace_and_scalar_force_retained=True,
                           first_order_local_EFT_only=True,
                           known_bosonic_operator_sector_not_full_fermionic_loop_action=True,
                           actual_original_mass_shift_checked=True,
                           no_quantum_constraint_or_GR_derivation_claim=True,
                           no_factorization_of_quantum_stress_expectations=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=600,tests=3,all_passed=True,tensor_error=result['tensor']['max_error'],
                         induced_dirac=result['force']['induced_dirac_coefficient_norm'],
                         contact_error=result['contact']['max_factorization_error'])))
