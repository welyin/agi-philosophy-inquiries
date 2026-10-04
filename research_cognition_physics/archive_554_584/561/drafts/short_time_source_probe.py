"""561 preliminary: full-source polynomial domains and short-time certificates.

No evolution truncation. Norm upper bounds use the Gaussian monomial triangle
inequality. Their floating evaluations are not interval-arithmetic certificates.
"""
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
if str(BASE) not in sys.path:sys.path.insert(0,str(BASE))
import joint_matter_energy_moment_control as p
import joint_gauge_link_reference as link
import joint_finite_time_gauge_probe as old

TARGET=HERE/'short_time_source_probe_results.json'


def build():
    L,u,c=old.constants();d=10;a=.5;b=1.;omega=1.;variance=.5;k=.8;eta=.3
    q=[p.var(i,d) for i in range(d)]
    radii=[p.add(*(p.mul(q[j],q[j]) for j in range(o,o+4))) for o in (0,5)]
    W={}
    for j,o in enumerate((0,5)):
        delta=[p.add(radii[j],p.const(-u[0],d)),p.add(p.mul(q[o+4],q[o+4]),p.const(-u[1],d))]
        W=p.add(W,p.scale(p.quadratic(delta,L),.25))
    diffs=[p.add(q[i],p.scale(q[5+i],-1)) for i in range(4)]
    F=p.scale(p.add(*(p.mul(z,z) for z in diffs)),k/2)
    V=p.add(W,F)
    totalq2=p.add(*(p.mul(z,z) for z in q))
    hlocal=p.add(p.const(d*a*omega,d),p.scale(totalq2,-a*omega*omega),V)
    def G(P,j):
        return p.add(*(p.scale(p.mul(q[5+n],p.deriv(P,5+m)),link.J[j][m,n])
            for m in range(4) for n in range(4) if link.J[j][m,n]))
    df=[p.deriv(F,i) for i in range(d)]
    dg=[G(F,j) for j in range(3)]
    lapF=p.add(p.scale(p.add(*(p.deriv(df[i],i) for i in range(d))),a),
        p.scale(p.add(*(G(dg[j],j) for j in range(3))),b))
    Gamma=p.add(p.scale(p.add(*(p.mul(z,z) for z in df)),a),p.scale(p.add(*(p.mul(z,z) for z in dg)),b))
    def H(P):
        lap=p.add(*(p.deriv(p.deriv(P,i),i) for i in range(d)))
        drift=p.add(*(p.mul(q[i],p.deriv(P,i)) for i in range(d)))
        group=p.add(*(G(G(P,j),j) for j in range(3)))
        return p.add(p.mul(hlocal,P),p.scale(lap,-a),p.scale(drift,2*a*omega),p.scale(group,-b))
    def J(P):
        grad=p.add(*(p.mul(df[i],p.add(p.deriv(P,i),p.scale(p.mul(q[i],P),-omega))) for i in range(d)))
        group=p.add(*(p.mul(dg[j],G(P,j)) for j in range(3)))
        return p.scale(p.add(p.mul(lapF,P),p.scale(grad,2*a),p.scale(group,2*b)),1j)
    P=p.add(p.const(1,d),p.scale(p.add(F,p.const(-4*k*variance,d)),1j*eta))
    norm=1+eta*eta*8*k*k*variance*variance
    def upper(P):
        return sum(abs(value)*math.sqrt(math.prod(p.normal_moment(2*n,variance) for n in powers))
                   for powers,value in P.items())/math.sqrt(norm)
    return dict(H=H,J=J,G=G,F=F,Gamma=Gamma,P=P,upper=upper,norm=norm,c=c,
        omega=omega,eta=eta,variance=variance,k=k,a=a,b=b,hlocal=hlocal)


def run():
    model=build();H,J,P=model['H'],model['J'],model['P'];Gamma=model['Gamma'];up=model['upper']
    HP=H(P);HHP=H(HP);JP=J(P);JHP=J(HP)
    GP=p.mul(Gamma,P);GHP=p.mul(Gamma,HP)
    # Independent commutator expression for the first-order operator.
    by_comm=p.scale(p.add(p.mul(model['F'],HP),p.scale(H(p.mul(model['F'],P)),-1)),1j)
    diff=p.add(JP,p.scale(by_comm,-1))
    commutator_residual=max([abs(v) for v in diff.values()]+[0.])
    assert commutator_residual<2e-13
    E=p.real(p.expect(p.mul(p.conj(P),HP),model['variance']))/model['norm']
    source_norm=p.norm2(P,model['variance'])/model['norm']
    assert abs(source_norm-1)<1e-12
    saved=json.loads((BASE/'joint_probe_time_resolution_results.json').read_text('utf8'))['full_source_clock_witness']
    assert abs(E-saved['full_source_E1'])<1e-10
    names={'Hpsi':HP,'H2psi':HHP,'Jpsi':JP,'JHpsi':JHP,'Gamma_psi':GP,'Gamma_Hpsi':GHP}
    bounds={name:up(polynomial) for name,polynomial in names.items()}
    # g=1, sigma_Q=.5: mu2=1, mu4=3. All are true full-source norms.
    K=bounds['Hpsi']+.5*bounds['Jpsi']+math.sqrt(3)/3*bounds['Gamma_psi']
    KH=bounds['H2psi']+.5*bounds['JHpsi']+math.sqrt(3)/3*bounds['Gamma_Hpsi']
    response_constant=KH+bounds['Hpsi']*K
    theta=2*model['k']*model['variance'];D=model['k']*(4*model['a']+3*model['b']*model['variance']/4)
    sigmaQ=.5;eta=model['eta'];norm=model['norm'];g=1.
    lower=(2*eta*D*g/(norm*sigmaQ*math.sqrt(2*math.pi)))*(7*theta/3)*math.exp(-2-g*g*theta*theta/(2*sigmaQ*sigmaQ))
    tau=lower/(8*response_constant)
    slope_lower=lower-2*tau*response_constant
    c_mass=1/model['c']['B1'];M=c_mass*tau*tau;energy=1/(2*M)
    assert slope_lower>=.749*lower and tau*tau>M*model['c']['B1']/2
    nodes,w=np.polynomial.legendre.leggauss(160)
    f=8*theta*(nodes+1);weights=8*theta*w
    gaussian_derivative=np.exp(-(f-2*theta)**2/(2*sigmaQ*sigmaQ))/(sigmaQ*math.sqrt(2*math.pi))
    integrand=D*f*(f/theta**2)*np.exp(-f/theta)*gaussian_derivative
    diagnostic=2*eta/norm*float(weights@integrand)
    assert diagnostic>lower
    return dict(status='preliminary_561_not_frozen',source_E1=E,normalization=source_norm,
        polynomial_terms={name:len(polynomial) for name,polynomial in names.items()},
        norm_upper_bounds=bounds,commutator_coefficient_residual=float(commutator_residual),
        Kpsi_upper=K,KHpsi_upper=KH,time_derivative_error_coefficient=2*response_constant,
        instantaneous_threshold_slope_positive_bound=lower,
        instantaneous_threshold_slope_quadrature_diagnostic=diagnostic,
        sufficient_chosen_duration=tau,finite_duration_slope_lower_bound=slope_lower,
        sufficient_mass=M,initial_pointer_energy=energy,pointer_mass_coefficient=c_mass,
        source_state_norm_error_upper=tau*K,
        scope=dict(full_source_polynomial_norm_bounds_not_Galerkin=True,
            positive_event_slope_not_four_coordinates=True,
            extremely_conservative_sufficient_budget_not_optimal=True,
            actual_final_compensated_quadrature_is_input=True,
            floating_evaluations_not_machine_interval_proofs=True))


if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
