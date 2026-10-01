"""563 preliminary: direct radial pointer on the unchanged full polynomial source.

Weighted Gaussian monomial norms include inverse radii in J_R. These are
analytic triangle bounds with floating evaluations, not interval certificates.
"""
import importlib.util
import json
import math
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
if str(BASE) not in sys.path:sys.path.insert(0,str(BASE))
import joint_matter_energy_moment_control as p
spec=importlib.util.spec_from_file_location('source561',BASE/'round561_drafts/short_time_source_probe.py')
source=importlib.util.module_from_spec(spec);spec.loader.exec_module(source)
TARGET=HERE/'radial_pointer_source_probe_results.json'


def run():
    m=source.build();P=m['P'];H=m['H'];HP=H(P);a=m['a'];t=m['variance'];N=m['norm']
    d=10;q=[p.var(i,d) for i in range(d)]
    S=p.add(*(p.mul(q[i],q[i]) for i in (*range(4),*range(5,9))))
    def moment_norm(powers,alpha,beta):
        first=sum(powers[:4]);second=sum(powers[5:9])
        raw=math.prod(p.normal_moment(2*n,t) for n in powers)
        ratio=(2*t)**(alpha+beta)*math.gamma(2+first+alpha)/math.gamma(2+first)
        ratio*=math.gamma(2+second+beta)/math.gamma(2+second)
        return math.sqrt(raw*ratio/N)
    def up(blocks):
        return sum(abs(c)*moment_norm(power,alpha,beta)
            for poly,alpha,beta in blocks for power,c in poly.items())
    def radial_ops(Q):
        Dx=p.add(*(p.mul(q[i],p.deriv(Q,i)) for i in range(4)))
        Dy=p.add(*(p.mul(q[i],p.deriv(Q,i)) for i in range(5,9)))
        base=p.add(p.scale(Q,8),p.scale(Dx,2),p.scale(Dy,2),p.scale(p.mul(S,Q),-2))
        # J_R = i a [base -(r2/r1)(3Q+2Dx) -(r1/r2)(3Q+2Dy)+4r1r2 Q].
        blocks=[(p.scale(base,1j*a),0,0),
            (p.scale(p.add(p.scale(Q,3),p.scale(Dx,2)),-1j*a),-1,1),
            (p.scale(p.add(p.scale(Q,3),p.scale(Dy,2)),-1j*a),1,-1),
            (p.scale(Q,4j*a),1,1)]
        gamma=[(p.scale(p.mul(S,Q),2*a),0,0),(p.scale(Q,-4*a),1,1)]
        return up(blocks),up(gamma)
    Jpsi,Gpsi=radial_ops(P);JHp,GHp=radial_ops(HP)
    Hnorm=m['upper'](HP);H2norm=m['upper'](H(HP))
    K=Hnorm+Jpsi/2+math.sqrt(3)*Gpsi/3
    KH=H2norm+JHp/2+math.sqrt(3)*GHp/3
    C=KH+Hnorm*K
    # Radial Gaussian law: r^2 has Gamma(shape=2,scale=2t).
    r1,r2=p.var(0,2),p.var(1,2)
    r12=p.mul(r1,r1);r22=p.mul(r2,r2);rsum=p.add(r12,r22)
    R=p.add(p.scale(rsum,.5),p.scale(p.mul(r1,r2),-1))
    f2=p.scale(p.add(p.mul(r12,r12),p.mul(r22,r22),p.scale(p.mul(r12,r22),3),
        p.scale(rsum,-16*t),p.const(64*t*t,2)),m['k']**2/4)
    weight=p.add(p.const(1,2),p.scale(f2,m['eta']**2))
    def radial_expect(poly):
        return p.real(sum(c*math.prod((2*t)**(n/2)*math.gamma(2+n/2) for n in power)
                          for power,c in poly.items()))
    assert abs(radial_expect(weight)-N)<1e-12
    rmean=radial_expect(p.mul(R,weight))/N
    assert rmean>0
    # A positive source-independent rectangle contribution to E_phi[R G'(R)].
    sigmaQ=.5;y0=.4;Rmin=.125;Rmax=.78125
    radius_probability_lower=(.5*2*math.exp(-2.25))*(.25*2*.25**3*math.exp(-.25))
    derivative_lower=math.exp(-max(abs(Rmin-y0),abs(Rmax-y0))**2/(2*sigmaQ**2))/(sigmaQ*math.sqrt(2*math.pi))
    slope_lower=4*a*m['k']*m['eta']/N*Rmin*derivative_lower*radius_probability_lower
    tau=slope_lower/(8*C)
    D1=4/(m['c']['v']*m['c']['ell']);D0=2*m['c']['uh']**2
    M=tau*tau/D1
    assert tau>0 and M>0
    return dict(status='preliminary_563_not_frozen',same_source_normalization=N,
        same_source_E1=source.run()['source_E1'],source_radial_mean=rmean,
        instantaneous_nonselective_matter_energy_increase=a*rmean/sigmaQ**2,
        instantaneous_nonselective_electric_energy_increase=0.,
        full_source_norm_upper=dict(Hpsi=Hnorm,H2psi=H2norm,J_R_psi=Jpsi,
            J_R_Hpsi=JHp,Gamma_R_psi=Gpsi,Gamma_R_Hpsi=GHp),
        Kpsi=K,KHpsi=KH,response_error_coefficient=2*C,
        actual_threshold_derivative_positive_lower=slope_lower,
        sufficient_duration=tau,finite_duration_slope_lower=3*slope_lower/4,
        R_squared_constants=dict(D0=D0,D1=D1),pointer_mass=M,
        initial_pointer_energy=1/(2*M),total_H_lower=-D0/(2*D1),
        scope=dict(RP_coupling_and_compensated_reader_are_new_instrument_inputs=True,
            inverse_radial_factors_kept_in_full_Cartesian_norms=True,
            no_new_R_phase_source_or_frozen_angles=True,
            finite_time_electric_energy_not_claimed_unchanged=True,
            no_spacetime_or_unified_completion=True))


if __name__=='__main__':
    result=run()
    if TARGET.exists():assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
