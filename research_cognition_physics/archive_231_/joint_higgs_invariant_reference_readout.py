"""556: radial domain obstruction, invariant records and instrument leakage.

Single-cell global SU(2) singlet mechanics with the inherited quartic potential.
No local gauge field, spacetime chart, autonomous apparatus or quantum gravity.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np
import joint_matter_energy_moment_control as p

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_higgs_invariant_reference_readout_results.json'


def radial_probe():
    path=HERE/'round556_drafts/radial_domain_probe.py'
    spec=importlib.util.spec_from_file_location('radial_probe_556',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    result=mod.run()
    assert result==json.loads((path.parent/'radial_domain_probe_results.json').read_text('utf8'))
    return result


def source_moments(k,n,omega=1.3,v=1.):
    """Normalized radial Laguerre state times a singlet Hermite state."""
    d=5;m=4;hbar=1.;a=1/(2*v);variance=1/(2*omega)
    q=[p.var(i,d) for i in range(d)]
    R=p.add(*(p.mul(x,x) for x in q[:4]));s2=p.mul(q[4],q[4])
    x=p.scale(R,omega)
    if k==0: radial=p.const(1,d)
    elif k==1: radial=p.scale(p.add(p.const(2,d),p.scale(x,-1)),1/math.sqrt(2))
    elif k==2: radial=p.scale(p.add(p.scale(p.mul(x,x),.5),p.scale(x,-3),p.const(3,d)),1/math.sqrt(3))
    else: raise ValueError(k)
    P=p.mul(radial,p.hermite_wave(n,q[4],omega,d))
    assert abs(p.norm2(P,variance)-1)<1e-10
    L,C,u=p.matter()
    delta=[p.add(R,p.const(-u[0],d)),p.add(s2,p.const(-u[1],d))]
    W=p.scale(p.quadratic(delta,L),v/4)
    grads=[];seconds=[]
    for i in range(d):
        first=p.deriv(P,i);g=p.scale(q[i],-omega)
        grads.append(p.add(first,p.mul(g,P)))
        seconds.append(p.add(p.deriv(first,i),p.scale(p.mul(g,first),2),
                       p.mul(p.add(p.mul(g,g),p.const(-omega,d)),P)))
    Tpsi=p.scale(p.add(*seconds),-a);Wpsi=p.mul(W,P);Hpsi=p.add(Tpsi,Wpsi)
    E1=p.real(p.expect(p.mul(p.conj(P),Hpsi),variance));E2=p.norm2(Hpsi,variance)
    lapW=p.add(*(p.deriv(p.deriv(W,i),i) for i in range(d)))
    wanted=p.scale(p.add(p.scale(R,6*L[0,0]+L[0,1]),
                         p.scale(s2,4*L[0,1]+3*L[1,1]),
                         p.const(-4*C[0]-C[1],d)),v)
    lap_error=max([abs(c) for c in p.add(lapW,p.scale(wanted,-1)).values()]+[0.])
    assert lap_error<1e-12
    identity=p.norm2(Tpsi,variance)+p.norm2(Wpsi,variance)
    identity+=2*a*sum(p.real(p.expect(p.mul(W,p.mul(p.conj(g),g)),variance)) for g in grads)
    identity-=a*p.real(p.expect(p.mul(p.mul(p.conj(P),P),lapW),variance))
    assert abs(identity-E2)<1e-8*max(1.,E2)
    ell=float(np.linalg.eigvalsh(L)[0]);r=1.
    alpha=max(0.,6*L[0,0]+L[0,1],4*L[0,1]+3*L[1,1])
    A=2*alpha/(ell*r);B=v*alpha*(sum(u)+r)+v*(4*abs(C[0])+abs(C[1]))
    G=E2+a*A*E1+a*B;S=sum(u)+r+2*E1/(v*ell*r)
    assert p.norm2(Tpsi,variance)+p.norm2(Wpsi,variance)<=G+1e-8
    operators=[p.mul(R,P),p.mul(q[4],P),p.scale(p.add(*seconds[:4]),-1),p.scale(seconds[4],-1)]
    means=[p.real(p.expect(p.mul(p.conj(P),f),variance)) for f in operators]
    second=[p.norm2(f,variance) for f in operators]
    tq=1/(2*omega);tp=omega/2
    record=[second[0]+4*tq*means[0]+2*m*tq*tq,
            second[1]+tq,
            second[2]+4*tp*means[2]+2*m*tp*tp,
            second[3]+4*tp*means[3]+2*tp*tp]
    bounds=[2*u[0]**2+8*E1/(v*ell)+4*tq*S+2*m*tq*tq,
            S+tq,4*v*v*G+8*v*tp*E1+2*m*tp*tp,
            4*v*v*G+8*v*tp*E1+2*tp*tp]
    assert all(x<=b+1e-8*max(1.,b) for x,b in zip(record,bounds))
    return dict(radial_excitation=k,singlet_oscillator_excitation=n,omega=omega,v=v,
                E1=E1,E2=E2,laplacian_polynomial_residual=float(lap_error),
                graph_identity_relative_residual=abs(identity-E2)/max(1.,E2),
                source_menu_means=means,source_menu_second_moments=second,
                actual_record_second_moments=record,
                record_bounds=[float(z) for z in bounds],
                potential_majorant=dict(A=float(A),B=float(B)),graph_budget=float(G))


def independent_bargmann_output(k,n,omega):
    """Positive coherent-state output density, integrated as a polynomial.

The radial creation polynomial is (sum alpha_i^2)^k with squared norm
4^k k! (2)_k. The singlet excitation contributes |alpha_s|^(2n)/n!.
All ten real integration variables are standard normal, unlike the source
wave-function calculation in five coordinates.
"""
    d=10;z=[p.var(i,d) for i in range(d)]
    alpha=[p.scale(p.add(z[i],p.scale(z[i+5],1j)),1/math.sqrt(2)) for i in range(5)]
    invariant=p.add(*(p.mul(x,x) for x in alpha[:4]))
    amplitude=p.const(1,d)
    for _ in range(k): amplitude=p.mul(amplitude,invariant)
    for _ in range(n): amplitude=p.mul(amplitude,alpha[4])
    denom=4**k*math.factorial(k)*math.factorial(k+1)*math.factorial(n)
    density=p.scale(p.mul(p.conj(amplitude),amplitude),1/denom)
    normalization=p.real(p.expect(density,1.))
    assert abs(normalization-1)<1e-9
    tq=1/(2*omega);tp=omega/2
    menu=[p.add(p.scale(p.add(*(p.mul(z[i],z[i]) for i in range(4))),1/omega),p.const(-4*tq,d)),
          p.scale(z[4],1/math.sqrt(omega)),
          p.add(p.scale(p.add(*(p.mul(z[i+5],z[i+5]) for i in range(4))),omega),p.const(-4*tp,d)),
          p.add(p.scale(p.mul(z[9],z[9]),omega),p.const(-tp,d))]
    means=[p.real(p.expect(p.mul(density,f),1.)) for f in menu]
    seconds=[p.real(p.expect(p.mul(density,p.mul(f,f)),1.)) for f in menu]
    return normalization,means,seconds


def compositions(n,d):
    if d==1: return [(n,)]
    return [(j,)+tail for j in range(n+1) for tail in compositions(n-j,d-1)]


def real_rep(matrix):
    result=np.zeros((4,4))
    for i in range(2):
        for j in range(2):
            a=matrix[i,j]
            result[2*i:2*i+2,2*j:2*j+2]=[[a.real,-a.imag],[a.imag,a.real]]
    return result


def su2_generators():
    pauli=(np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.]))
    return [1j*real_rep(-.5j*s) for s in pauli]


def fock_generators(n):
    basis=compositions(n,4);index={x:i for i,x in enumerate(basis)}
    result=[]
    for small in su2_generators():
        large=np.zeros((len(basis),len(basis)),complex)
        for j,occupation in enumerate(basis):
            for b in range(4):
                if not occupation[b]: continue
                middle=list(occupation);middle[b]-=1
                for a in range(4):
                    if small[a,b]==0: continue
                    final=middle.copy();final[a]+=1
                    large[index[tuple(final)],j]+=small[a,b]*math.sqrt(occupation[b]*(middle[a]+1))
        result.append(large)
    return result


def run():
    checks=[];old=radial_probe()
    checks.append('same_potential_radial_domain_counterexample_with_independent_log_integral')
    rows=[source_moments(k,n,v=v) for k,n,v in ((0,0,1.),(1,1,1.),(2,0,.7))]
    checks.append('full_five_component_interacting_energy_and_correct_laplacian_majorant')
    residuals=[]
    for row in rows:
        norm,means,seconds=independent_bargmann_output(row['radial_excitation'],
                                                     row['singlet_oscillator_excitation'],row['omega'])
        residual=max(np.max(abs(np.array(means)-row['source_menu_means'])),
                     np.max(abs(np.array(seconds)-row['actual_record_second_moments'])))
        assert residual<2e-8
        residuals.append(float(residual))
    checks.append('positive_invariant_POVM_moments_independent_Bargmann_output_integration')
    rng=np.random.default_rng(556);maxinv=0.
    for _ in range(40):
        x=rng.normal(size=4);x/=np.linalg.norm(x)
        a=x[0]+1j*x[1];b=x[2]+1j*x[3]
        U=np.array([[a,b],[-np.conj(b),np.conj(a)]])
        R=real_rep(U);assert np.max(abs(R.T@R-np.eye(4)))<1e-12
        X=rng.normal(size=4);P=rng.normal(size=4)
        maxinv=max(maxinv,abs((R@X)@(R@X)-X@X),abs((R@P)@(R@P)-P@P))
    assert maxinv<1e-12
    checks.append('SU2_real_representation_and_full_quadratic_output_invariance')
    sectors=[];total=Fraction(0);maxcomm=0.
    for n in range(9):
        J=fock_generators(n);casimir=sum(j@j for j in J)
        eig=np.linalg.eigvalsh(casimir)
        count=int(np.sum(abs(eig)<1e-9));assert count==int(n%2==0)
        comm=np.max(abs(J[0]@J[1]-J[1]@J[0]-1j*J[2]))
        maxcomm=max(maxcomm,float(comm));assert comm<1e-12
        total+=count*Fraction(1,2**(n+4))
        sectors.append(dict(excitation=n,dimension=len(eig),singlet_multiplicity=count,
                            smallest_Casimir_eigenvalue=float(eig[0])))
    tail=Fraction(1,12288);assert total+tail==Fraction(1,12)
    checks.append('finite_Fock_Casimirs_and_exact_thermal_singlet_leakage_with_analytic_tail')
    # The whole thermal-state normalization includes all Fock multiplicities.
    # Closed generating function: sum C(n+3,3)/2^(n+4)=1.
    assert Fraction(1,16)/(1-Fraction(1,2))**4==1
    probability=Fraction(1,12);leakage=1-probability
    odd=Fraction(40,81);assert leakage>=odd
    certificates=[]
    for row in rows:
        radius=math.sqrt(sum(row['record_bounds'])/.05)
        witnessed=sum(row['actual_record_second_moments'])/radius**2
        assert witnessed<=.05
        certificates.append(dict(radial_excitation=row['radial_excitation'],
                                 raw_invariant_record_radius=radius,
                                 failure_bound=.05,witness_Markov_bound=witnessed))
    checks.append('one_energy_budget_for_all_invariant_records_and_instrument_scope_comparison')
    deps=('joint_matter_energy_moment_control.py','joint_singlet_common_mass_rg_results.json',
          'round556_drafts/radial_domain_probe.py','round556_drafts/radial_domain_probe_results.json')
    return dict(round=556,tests_run=len(checks),failures=0,errors=0,checks=checks,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                radial_counterexample=old,witness_sources=rows,
                independent_record_moment_errors=residuals,SU2_menu_invariance_residual=float(maxinv),
                fock_sector_diagnostics=sectors,Casimir_commutator_residual=maxcomm,
                heterodyne_instrument=dict(input='invariant_reference_oscillator_vacuum_not_interacting_ground_state',
                    output='four_mode_thermal_mean_occupation_one_not_interacting_Gibbs_state',
                    singlet_probability=str(probability),leakage_probability=str(leakage),
                    checked_singlet_weight_through_n8=str(total),exact_singlet_tail=str(tail),
                    independent_center_parity_leakage_lower_bound=str(odd),
                    covariant_output_density_is_not_singlet_support=True),
                finite_record_certificates=certificates,
                scope=dict(single_cell_global_SU2_only=True,radial_quantization_is_inherited_math=True,
                           invariant_POVM_exists=True,finite_bin_square_root_instrument_is_abstract_input=True,
                           no_same_matter_autonomous_detector_proved=True,
                           no_post_measurement_energy_moment_bound_claimed=True,
                           no_local_Gauss_constraints_or_spatial_derivatives_included=True,
                           no_spacetime_coordinates_or_gravity_generation_claimed=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=556,tests=result['tests_run'],
                         max_record_residual=max(result['independent_record_moment_errors']),
                         singlet_leakage=result['heterodyne_instrument']['leakage_probability'])))
