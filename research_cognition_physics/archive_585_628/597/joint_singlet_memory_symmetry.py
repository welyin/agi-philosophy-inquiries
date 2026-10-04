"""597: original singlet wells, relational energy, and conditional memory.

First two groups use the original potential and hyperbolic target geometry.
The parity doublet test is an algebraic fixture, not original eigenvalues.
No autonomous writer or long lifetime is certified numerically.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_geometry_work_noise as common

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_singlet_memory_symmetry_results.json'


def wells_and_symmetry():
    L,u,_=original.lattice.scalar.parameters()
    a,b,d=L[0,0],L[0,1],L[1,1]
    A=6*original.M-u[0];v=u[1]
    x=u[0]+v*(b*A-d*v)/(a*A-b*v)
    barrier=9*v*v*np.linalg.det(L)/(a*A*A-2*b*v*A+d*v*v)
    assert 0<x<6*original.M and a*A-b*v>0
    cross=np.array([np.sqrt(x),0.,0.,0.,0.])
    assert abs(original.node_potential(cross)-barrier)<1e-15
    # Independent derivative and sampling; global minimality is analytic.
    diffs=[]
    for step in (1e-3,5e-4,2.5e-4):
        def f(y):return original.node_potential(np.array([np.sqrt(y),0,0,0,0.]))
        diffs.append(float(abs((f(x+step)-f(x-step))/(2*step))))
    assert diffs[-1]<diffs[0]/12
    vac=np.array([np.sqrt(u[0]),0,0,0,np.sqrt(u[1])])
    reflected=vac.copy();reflected[4]*=-1
    assert original.node_potential(vac)<1e-28
    assert original.node_potential(reflected)<1e-28
    rng=np.random.default_rng(597)
    phi=rng.normal(size=(48,5))*.3;chi=rng.normal(size=(48,5))*.4
    P=np.diag([1,1,1,1,-1]);pp=phi@P;cc=chi@P
    inv_error=float(np.max(abs(original.inverse(pp)-P@original.inverse(phi)@P)))
    metric_error=float(np.max(abs(original.metric(pp)-P@original.metric(phi)@P)))
    potential_error=float(np.max(abs(original.node_potential(pp)-original.node_potential(phi))))
    distance_error=float(np.max(abs(original.distance_squared(pp,cc)-original.distance_squared(phi,chi))))
    assert max(inv_error,metric_error,potential_error,distance_error)<1e-13
    fv=original.F(vac)
    edge=original.distance_squared(vac,reflected)/2
    assert abs(edge-3*np.arccosh(1+u[1]/(3*fv))**2)<1e-14
    return dict(original_M=original.M,L=L.tolist(),vacuum_squared=u.tolist(),
                vacuum_F=float(fv),crossing_h_squared=float(x),node_crossing_barrier=float(barrier),
                crossing_derivative_errors=diffs,opposite_vacuum_edge_energy_per_k=float(edge),
                symmetry_errors=dict(inverse=inv_error,metric=metric_error,
                                     potential=potential_error,distance=distance_error))


def radial_packet(n):
    # Gauge-invariant compact packet. Original measure h^3 sqrt(M) F^-3;
    # amplitude chi F^(3/2)/M^(1/4) leaves h^3 chi^2 in probabilities.
    _,u,_=original.lattice.scalar.parameters()
    z,w=np.polynomial.legendre.leggauss(n)
    h=np.sqrt(u[0])+.04*z[:,None]
    s=np.sqrt(u[1])+.03*z[None,:]
    h,s=np.broadcast_arrays(h,s)
    prob=w[:,None]*w[None,:]*h**3*np.exp(-2/(1-z[:,None]**2)-2/(1-z[None,:]**2))
    prob/=prob.sum()
    return h.ravel(),s.ravel(),prob.ravel()


def relational_packet_energy():
    rows=[]
    for n in (8,12,16):
        h,s,p=radial_packet(n);F=original.M-(h*h+s*s)/6
        na=3*n
        angle=np.arange(1,na+1)*np.pi/(na+1)
        cos=np.cos(angle);aw=2/(na+1)*np.sin(angle)**2
        assert abs(aw.sum()-1)<1e-14
        denom=np.sqrt(F[:,None]*F[None,:]);hh=h[:,None]*h[None,:];ss=s[:,None]*s[None,:]
        pair=p[:,None]*p[None,:]
        energies=[]
        for sign in (1,-1):
            total=0.
            for z,w in zip(cos,aw):
                excess=(original.M-(hh*z+sign*ss)/6)/denom-1
                assert excess.min()>-1e-13
                d2=24*np.arcsinh(np.sqrt(np.maximum(excess,0)/2))**2
                total+=w*np.sum(pair*d2)/2
            energies.append(float(total))
        delta=energies[1]-energies[0]
        assert delta>0
        rows.append(dict(nodes_per_radial_axis=n,angular_nodes=na,
                         same_sign_edge_energy_per_k=energies[0],
                         opposite_sign_edge_energy_per_k=energies[1],
                         relational_excess_per_k=delta,conformal_source_excess_per_k=2*delta))
    assert abs(rows[-1]['relational_excess_per_k']-rows[-2]['relational_excess_per_k'])<2e-5
    return dict(rows=rows,packet_halfwidths=dict(h=.04,s=.03),
                quantum_configuration_integral_not_propagation=True,
                isotropic_one_edge_coefficient_k_declared=True,
                point_well_and_packet_energy_distinguished=True)


def conditional_doublet():
    # Original gap and local contrast are unknown. Explicit diagnostic input.
    hbar=common.HBAR;gap=.07;c=.96;error=.05
    H=np.diag([.3,.3+gap]);P=np.diag([1.,-1.])
    Q=np.array([[0.,c],[c,0.]])
    G=np.diag([.8,1.1])
    plus=np.ones(2,dtype=complex)/np.sqrt(2);minus=P@plus
    T=hbar/gap*np.arccos((1-2*error)/c)
    rows=[]
    for t in (0.,T/2,T,2*T):
        U=np.diag(np.exp(-1j*np.diag(H)*t/hbar));v=U@plus;w=U@minus
        success=float(((v.conj()@((np.eye(2)+Q)/2)@v)+(w.conj()@((np.eye(2)-Q)/2)@w)).real/2)
        expected=(1+c*np.cos(gap*t/hbar))/2
        assert abs(success-expected)<1e-14
        assert abs(v.conj()@G@v-w.conj()@G@w)<1e-14
        rows.append(dict(time=float(t),success=success))
    assert rows[2]['success']>=1-error-1e-14 and rows[3]['success']<1-error
    # A symmetric blank remains symmetric; it cannot make a polarized state.
    blank=np.diag([.6,.4]);out=U@blank@U.conj().T
    assert np.linalg.norm(P@out@P-out)<1e-14
    assert abs(np.trace(Q@out))<1e-14
    return dict(diagnostic_only=True,gap_input=gap,contrast_input=c,
                allowed_error=error,conditional_window=float(T),rows=rows,
                original_gap_contrast_and_preparation_not_computed=True)


def run():
    wells=wells_and_symmetry();packet=relational_packet_energy();doublet=conditional_doublet()
    deps=('research_note_574.md','research_note_577.md','research_note_590.md','research_note_591.md',
          'research_note_596.md','joint_curved_quantum_source.py','joint_geometry_work_noise.py',
          'closed_time_path_bridge_review_586.md')
    return dict(round=597,tests_run=3,failures=0,errors=0,wells=wells,packet=packet,doublet=doublet,
                scope=dict(original_global_sign_symmetry=True,
                           geometry_blindness_only_to_global_sign_related_states=True,
                           relative_sign_records_can_cost_edge_energy=True,
                           long_lived_memory_and_autonomous_writing_not_established=True,
                           no_original_spectral_doublet_or_GR_derivation_claim=True),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=597,tests=3,all_passed=True,packet=result['packet']['rows'])))
