"""Independent 1247 checks: exact early-time coefficients, full exchange,
and constraint-reduced Decimal evolution. No author's interval code imported.

The independent formal arithmetic is reused from our frozen 1246 checker.
Decimal step convergence is diagnostic, not an enclosure certificate.
Default replay is read-only; --write exclusively creates the result.
"""
from __future__ import annotations

import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
FORMAL_TOOL = HERE.parent/'1246'/'independent_check.py'
spec = importlib.util.spec_from_file_location('independent_formal_1246', FORMAL_TOOL)
old = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(old)
# These alter only this process's arithmetic configuration, never the old file.
old.RAD = F(111, 100)
old.ORDER = 12
Q, P = old.Q, old.P
ORDER = old.ORDER
RAD = old.RAD
T = P([0, 1])
W = 64*T**3-192*T**4+192*T**5-64*T**6
WP = W.derivative()
FNAMES = ('u','H','p','v','J','q','r','z','L','psi','X','V','Y','K')
RESULT = HERE/'independent_results.json'


def formal_energy(y, alpha):
    u,H,p,v,J,q,r,z,L,psi,X,V,Y,K = y
    mass = F(71,100)+J+(q*q+r*r+K*K)/2+alpha*W*q*p/4
    rods = F(3,20)*(1+p/5)*u*u
    dust = (mass+F(1,50)+p/100)*u**3
    rho = (v*v+p*p+V*V+X*X)/2+rods+dust
    rhoP = v*v+V*V+F(2,3)*rods+dust
    return mass,rho,rhoP


def formal_rhs(y, alpha):
    u,H,p,v,J,q,r,z,L,psi,X,V,Y,K = y
    mass,rho,rhoP = formal_energy(y,alpha)
    return [-H*u,-F(3,2)*rhoP,v,
            -3*H*v-p-F(3,100)*u*u-(F(1,100)+alpha*W*q/4)*u**3,
            -alpha*WP*q*p/4,r,-q-alpha*W*p/4,p/10,u,W*q*p/4,
            V,-3*H*V-X,K,P(0)]


def formal_branch(alpha):
    y=[P(a) for a in (1,Q(0,1),0,0,F(1,10),F(1,2),0,0,0,0,F(1,10),0,0,0)]
    for n in range(ORDER):
        rhs=formal_rhs(y,alpha)
        for a,b in zip(y,rhs):
            a.c[n+1]=b.c[n]/(n+1)
    mass,rho,rhoP=formal_energy(y,alpha)
    named=dict(zip(FNAMES,y))
    exchange=mass.derivative()-alpha*W*named['q']*named['v']/4
    continuity=rho.derivative()+3*named['H']*rhoP
    constraint=named['H']**2-rho
    assert all(c==0 for c in exchange.c[:ORDER])
    assert all(c==0 for c in continuity.c[:ORDER])
    assert all(c==0 for c in constraint.c)
    return named|{'mass':mass,'rho':rho}


def source_point(stage, alpha, theta):
    """Independent chain-rule evaluation, all arithmetic exactly rational."""
    u,H,p,v,J,q,r,X,V,Y,K = map(F,('4/5','7/10','-1/50','-1/40',
                                  '1/10','2/5','-1/5','1/10','-1/30','1/200','-1/300'))
    pulse=lambda s:(64*s**3*(1-s)**3,192*s**2*(1-s)**2*(1-2*s))
    w,wp=pulse(theta) if stage==0 else (F(0),F(0))
    f,fp=pulse(theta-1) if stage==1 else (F(0),F(0))
    mass=F(71,100)+J+(q*q+r*r+K*K)/2+alpha*w*q*p/4+f*X*Y/10
    rod=F(3,20)*(1+p/5)*u*u
    dust=(mass+F(1,50)+p/100)*u**3
    rhoP=v*v+V*V+F(2,3)*rod+dust
    up=-H*u
    Hp=-F(3,2)*rhoP
    vp=-3*H*v-p-F(3,100)*u*u-(F(1,100)+alpha*w*q/4)*u**3
    Vp=-3*H*V-X-f*Y*u**3/10
    Jp=-alpha*wp*q*p/4-fp*X*Y/10
    qp,rp,Yp,Kp=r,-q-alpha*w*p/4,K,-f*X/10
    mdot=Jp+q*qp+r*rp+K*Kp+alpha*(wp*q*p+w*qp*p+w*q*v)/4
    mdot+=(fp*X*Y+f*V*Y+f*X*Yp)/10
    assert mdot==alpha*w*q*v/4+f*Y*V/10
    rod_dot=F(3,100)*v*u*u-2*H*rod
    dust_dot=(mdot+v/100)*u**3-3*H*dust
    rho_dot=v*vp+p*v+V*Vp+X*V+rod_dot+dust_dot
    assert rho_dot+3*H*rhoP==0
    assert 2*H*Hp-rho_dot==0
    missing_X_force_defect=f*Y*V*u**3/10
    missing_clock_reaction_defect=fp*X*Y*u**3/10
    wrong_Vp=Vp+f*Y*u**3/10
    deleted_force_rho_dot=v*vp+p*v+V*wrong_Vp+X*V+rod_dot+dust_dot
    assert deleted_force_rho_dot+3*H*rhoP==missing_X_force_defect
    wrong_Jp=Jp+fp*X*Y/10
    deleted_control_mdot=wrong_Jp+q*qp+r*rp+K*Kp
    deleted_control_mdot+=alpha*(wp*q*p+w*qp*p+w*q*v)/4
    deleted_control_mdot+=(fp*X*Y+f*V*Y+f*X*Yp)/10
    deleted_control_dust_dot=(deleted_control_mdot+v/100)*u**3-3*H*dust
    deleted_control_rho_dot=v*vp+p*v+V*Vp+X*V+rod_dot+deleted_control_dust_dot
    assert deleted_control_rho_dot+3*H*rhoP==missing_clock_reaction_defect
    if stage==1:
        assert missing_X_force_defect!=0
        if theta!=F(3,2): assert missing_clock_reaction_defect!=0
    return dict(stage=stage,alpha=alpha,theta=str(theta),mass_rate=str(mdot),
                complete_continuity_defect='0',constraint_rate='0',
                missing_X_force_defect=str(missing_X_force_defect),
                missing_clock_reaction_defect=str(missing_clock_reaction_defect))


def formal_group():
    off,on=formal_branch(0),formal_branch(1)
    expected={'H':(6,Q(-F(21,1000))),
              'u':(7,Q(F(3,1000))),
              'X':(9,Q(-F(7,80000)))}
    leading={}
    for name,(order,coefficient) in expected.items():
        d=on[name]-off[name]
        assert all(c==0 for c in d.c[:order]),name
        assert d.c[order]==coefficient,name
        leading[name]=dict(power=order,coefficient={
            'rational':str(coefficient.a),
            'sqrt_111_100_coefficient':str(coefficient.b)})
    assert off['rho'].c[0]==Q(RAD)==on['rho'].c[0]
    assert off['mass'].c[0]==Q(F(187,200))==on['mass'].c[0]
    assert all(c==0 for c in off['Y'].c+off['K'].c+on['Y'].c+on['K'].c)
    rows=[source_point(stage,a,t) for a in (0,1)
          for stage in (0,1) for t in (F(1,4)+stage,F(1,2)+stage,F(3,4)+stage)]
    # Rational guards for the off-pulse energy bootstrap. The differential
    # argument (positivity continuation and rho'<=0) is reviewed analytically.
    assert 2*RAD<F(3,2)**2
    assert RAD<F(11,10)**2
    assert 1-F(3,2)/5==F(7,10)
    assert F(1,5)-F(3,2)/10==F(1,20)
    assert 1-F(11,10)*F(1,4)==F(29,40)
    return dict(status='PASS',method='exact Q(sqrt(111/100)) series plus full chain rule',
                first_nonzero_differences=leading,initial_rho='111/100',
                initial_clock_mass='187/200',
                complete_exchange_and_continuity_coefficients_through_order11='0',
                constraint_coefficients_through_order12='0',
                memory_Y_K_identically_zero_before_write=True,source_point_checks=rows,
                hold_bootstrap_guards=dict(rho_upper='111/100',abs_p_X_upper='3/2',
                    H_upper='11/10',F_lower='7/10',M_r_lower='1/20',
                    u_retention_lower_factor='29/40',duration='1/4',
                    requires_positive_end_of_write_domain=True),
                scope='Early series is exact; point samples only check the analytic source identities.')


DNAMES=('a','p','v','J','q','r','z','L','psi','X','V','Y','K')


def wave(t):
    def pulse(s):
        return 64*s**3*(1-s)**3,192*s**2*(1-s)**2*(1-2*s)
    zero=D(0)
    w,wp=pulse(t) if zero<=t<=1 else (zero,zero)
    f,fp=pulse(t-1) if D(1)<=t<=2 else (zero,zero)
    return w,wp,f,fp


def decimal_energy(t,y,alpha):
    a,p,v,J,q,r,z,L,psi,X,V,Y,K=y
    w,wp,f,fp=wave(t)
    mass=D('.71')+J+(q*q+r*r+K*K)/2+D(alpha)*w*q*p/4+f*X*Y/10
    rods=D('.15')*(1+p/5)/a**2
    dust=(mass+D('.02')+p/100)/a**3
    rho=(v*v+p*p+V*V+X*X)/2+rods+dust
    return mass,rods,dust,rho,rho.sqrt()


def decimal_rhs(t,y,alpha):
    a,p,v,J,q,r,z,L,psi,X,V,Y,K=y
    w,wp,f,fp=wave(t)
    mass,rods,dust,rho,H=decimal_energy(t,y,alpha)
    return [a*H,v,-3*H*v-p-D('.03')/a**2-(D('.01')+D(alpha)*w*q/4)/a**3,
            -D(alpha)*wp*q*p/4-fp*X*Y/10,r,-q-D(alpha)*w*p/4,
            p/10,1/a,w*q*p/4,V,-3*H*V-X-f*Y/(10*a**3),K,-f*X/10]


def integrate_decimal(alpha,steps_per_unit):
    y=list(map(D,('1','0','0','.1','.5','0','0','0','0','.1','0','0','0')))
    h=D(1)/steps_per_unit
    minima={key:D(10) for key in ('J','M_c','F','M_r')}
    defect=D(0)
    at_write_end=None
    # Two pulse windows, followed by a quarter-unit retention diagnostic.
    total=9*steps_per_unit//4
    for n in range(total):
        t=D(n)/steps_per_unit
        f1=decimal_rhs(t,y,alpha)
        f2=decimal_rhs(t+h/2,[a+h*b/2 for a,b in zip(y,f1)],alpha)
        f3=decimal_rhs(t+h/2,[a+h*b/2 for a,b in zip(y,f2)],alpha)
        f4=decimal_rhs(t+h,[a+h*b for a,b in zip(y,f3)],alpha)
        a,p,v,J,q,r,z,L,psi,X,V,Y,K=y
        ap,pp,vp,Jp,qp,rp,zp,Lp,psip,Xp,Vp,Yp,Kp=f1
        mass,rods,dust,rho,H=decimal_energy(t,y,alpha)
        w,wp,f,fp=wave(t)
        mdot=Jp+q*qp+r*rp+K*Kp+D(alpha)*(wp*q*p+w*qp*p+w*q*pp)/4
        mdot+=(fp*X*Y+f*Xp*Y+f*X*Yp)/10
        rhop=v*vp+p*pp+V*Vp+X*Xp+D('.03')*pp/a**2-2*H*rods
        rhop+=(mdot+pp/100)/a**3-3*H*dust
        identity=rhop/(2*H)+D('1.5')*(v*v+V*V+2*rods/3+dust)
        defect=max(defect,abs(identity))
        for key,val in {'J':J,'M_c':mass,'F':1+p/5,'M_r':D('.2')+p/10}.items():
            minima[key]=min(minima[key],val)
        y=[a+h*(b+2*c+2*d+e)/6 for a,b,c,d,e in zip(y,f1,f2,f3,f4)]
        if n+1==2*steps_per_unit:
            mass,rods,dust,rho,H=decimal_energy(D(2),y,alpha)
            at_write_end=dict(zip(DNAMES,y))|{'H':H,'u':1/y[0],'Brod':D('.1')/y[0]**2,'M_c':mass}
    assert at_write_end is not None
    assert y[-1]==at_write_end['K'] # Every subsequent K RHS is exactly zero.
    assert abs(y[-2]-at_write_end['Y']-at_write_end['K']/4)<D('1e-55')
    assert defect<D('1e-50')
    assert min(minima.values())>0
    return at_write_end,minima,defect,dict(hold_until='9/4',K_exactly_constant=True,
                                        Y_linear_drift_is_not_record_loss=True)


def decimal_group():
    rows=[];deltas=[]
    with localcontext() as ctx:
        ctx.prec=60
        for steps in (256,512,1024):
            off,min0,res0,hold0=integrate_decimal(0,steps)
            on,min1,res1,hold1=integrate_decimal(1,steps)
            delta={k:on[k]-off[k] for k in ('K','Y','X','Brod','L','z','u','H')}
            rows.append(dict(steps_per_unit=steps,
                             off={k:str(v) for k,v in off.items()},
                             on={k:str(v) for k,v in on.items()},
                             difference={k:str(v) for k,v in delta.items()},
                             diagnostic_sampled_minima_off={k:str(v) for k,v in min0.items()},
                             diagnostic_sampled_minima_on={k:str(v) for k,v in min1.items()},
                             chain_rule_residual_max=str(max(res0,res1)),
                             postwrite_hold_off=hold0,postwrite_hold_on=hold1))
            deltas.append(delta)
        convergence={}
        for name in deltas[-1]:
            coarse=abs(deltas[1][name]-deltas[0][name])
            fine=abs(deltas[2][name]-deltas[1][name])
            assert fine<coarse/7 and fine<D('1e-10'),name
            assert abs(deltas[-1][name])>100*fine,name
            convergence[name]=dict(coarse_change=str(coarse),fine_change=str(fine),ratio=str(coarse/fine))
        assert deltas[-1]['K']>2*D('1e-10')
    return dict(status='PASS',method='60-digit Decimal RK4, a variable and H=sqrt(full rho)',
                rows=rows,convergence_diagnostics=convergence,
                no_rigorous_integration_error_bound_claimed=True,
                postwrite_hold_proof_is_Kdot_zero_not_numerical_stability=True,
                terminal_error_1e_10_is_an_unimplemented_classical_readout_contract=True)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(diagnose=False):
    result=dict(round=1247,status='PASS',scientific_groups=2,
                formal_group=formal_group(),decimal_group=decimal_group(),
                independent_formal_tool_sha256=sha(FORMAL_TOOL),
                CT_adopted=False,actual_quantum_readout_certified=False,whole_GR_proved=False)
    if not diagnose:
        files=('admission.md','proof.md','check.py','results.json','NEXT.md')
        result['author_sha256']={p:sha(HERE/p) for p in files}
        result['independent_code_sha256']=sha(Path(__file__))
        author=json.loads((HERE/'results.json').read_text('utf8'))
        assert author['status']=='PASS'
        result['author_status']='PASS'
        with localcontext() as ctx:
            ctx.prec=60
            lo,hi=map(D,author['group_finite_memory']['delta_K_final'])
            value=D(result['decimal_group']['rows'][-1]['difference']['K'])
            mesh=D(result['decimal_group']['convergence_diagnostics']['K']['fine_change'])
            distance=max(lo-value,value-hi,D(0))
            assert distance<2*mesh
            result['finite_certificate_comparison']=dict(author_interval=[str(lo),str(hi)],
                independent_distance_from_interval=str(distance),last_mesh_change=str(mesh),
                consistency_only_not_rigorous_error_estimate=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--diagnose',action='store_true')
    args=parser.parse_args()
    assert not(args.write and args.diagnose)
    value=run(args.diagnose)
    if args.write:
        with RESULT.open('x',encoding='utf8') as f:
            json.dump(value,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='INITIAL_SAVE'
    elif args.diagnose:
        mode='DIAGNOSTIC_ONLY'
    else:
        assert json.loads(RESULT.read_text('utf8'))==value,'read-only result mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status=value['status'],mode=mode,scientific_groups=2,
                         early=value['formal_group']['first_nonzero_differences'],
                         finite=value['decimal_group']['rows'][-1]['difference']),ensure_ascii=False))
