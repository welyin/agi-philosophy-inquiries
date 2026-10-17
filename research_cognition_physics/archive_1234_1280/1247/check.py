"""1247: fixed autonomous read pulse, complete source and retained memory.

The exact integer outward interval arithmetic is reused from 1242 through the
frozen 1246 module.  Each C2 pulse is integrated on its own polynomial side;
no ninth derivative is taken across theta=1 or theta=2.  Default is read-only.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / '1246' / 'check.py'
TOOLS = HERE.parent / '1242' / 'check.py'
spec = importlib.util.spec_from_file_location('round1246_interval_reuse', PREVIOUS)
old = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(old)
I, Jet, S, sqrt_interval, hull = old.I, old.Jet, old.S, old.sqrt_interval, old.hull
NAMES = old.NAMES + ['X', 'V', 'Y', 'K']


def pulses(theta, phase):
    if phase == 0:
        w, wp = old.pulse(theta)
        return w, wp, 0, 0
    if phase == 1:
        f, fp = old.pulse(theta-1)
        return 0, 0, f, fp
    assert phase == 2
    return 0, 0, 0, 0


def energy(y, alpha, phase):
    u,H,p,v,z,light,theta,J,q,r,psi,X,V,Y,K = y
    w,wp,f,fp = pulses(theta,phase)
    osc = (q*q+r*r)*F(1,2)
    memory = F(1,100)+K*K*F(1,2)+F(1,10)*f*X*Y
    Mc = F(7,10)+J+osc+alpha*F(1,4)*w*q*p+memory
    D = Mc+F(1,50)
    rod = F(3,20)*(1+F(1,5)*p)*u*u
    dust = (D+F(1,100)*p)*u**3
    scalar = (v*v+p*p+V*V+X*X)*F(1,2)
    rho = scalar+rod+dust
    pressure = (v*v-p*p+V*V-X*X)*F(1,2)-rod*F(1,3)
    return Mc,D,rod,dust,rho,pressure,memory,osc


def rhs(y, alpha, phase):
    u,H,p,v,z,light,theta,J,q,r,psi,X,V,Y,K = y
    w,wp,f,fp = pulses(theta,phase)
    Mc,D,rod,dust,rho,pressure,memory,osc = energy(y,alpha,phase)
    coupling = alpha*F(1,4)*w
    return [-H*u, -F(3,2)*(v*v+V*V+F(2,3)*rod+dust), v,
            -3*H*v-p-F(3,100)*u*u-(F(1,100)+coupling*q)*u**3,
            F(1,10)*p, u, 1,
            -alpha*F(1,4)*wp*q*p-F(1,10)*fp*X*Y,
            r, -q-coupling*p, F(1,4)*w*q*p,
            V, -3*H*V-X-F(1,10)*f*Y*u**3, K, -F(1,10)*f*X]


def coefficients(y, alpha, phase, order):
    c = [[I(a)] for a in y]
    for k in range(order):
        f = rhs([Jet(a) for a in c],alpha,phase)
        for j in range(len(y)):
            term = f[j].c[k] if isinstance(f[j],Jet) else I(f[j] if k==0 else 0)
            c[j].append(term.div(k+1))
    return c


def initial_interval():
    y = old.initial_interval()
    y[1] = sqrt_interval(F(111,100))
    return y+[I(F(1,10)),I(0),I(0),I(0)]


def endpoint_data(y, alpha, phase):
    Mc,D,rod,dust,rho,pressure,memory,osc = energy(y,alpha,phase)
    return dict(state={n:a.out() for n,a in zip(NAMES,y)},
                energy=dict(M_c=Mc.out(),oscillator=osc.out(),
                  memory_complete=memory.out(),memory_kinetic=(y[14]**2*F(1,2)).out(),
                  rho=rho.out(),pressure=pressure.out()),
                Brod=(y[0]**2*F(1,10)).out())


def integrate(alpha, steps_per_phase=200, order=8):
    dt = F(1,steps_per_phase)
    h,time = I(dt),I(0,dt)
    y = initial_interval()
    radii = [2*dt,3*dt,F(1,5)*dt,F(1,2)*dt,F(1,30)*dt,
             F(11,10)*dt,0,F(1,4)*dt,F(3,5)*dt,F(3,5)*dt,F(1,20)*dt,
             F(1,5)*dt,F(1,5)*dt,F(1,10)*dt,F(1,10)*dt]
    domain = [I('.2','1.02'),I('.1','1.08'),I('-.15','.03'),I('-.15','.05'),
              I('-.03','.002'),I('-.01','2.02'),I(0,2),I('.04','.17'),
              I('-.3','.52'),I('-.52','.03'),I('-.02','.01'),
              I('-.1','.12'),I('-.15','.03'),I('-.1','.03'),I('-.1','.03')]
    hp=[I(1)]
    for _ in range(order+1): hp.append(hp[-1]*h)
    all_tube=None;positivity=None;checkpoints={}
    max_width=max_remainder=max_constraint=F(0)
    for phase in [0,1]:
        assert y[6].lo==phase*S and y[6].hi==phase*S
        for step in range(steps_per_phase):
            box=[a+I(-r,r) for a,r in zip(y,radii)]
            box[6]=I(phase+F(step,steps_per_phase),phase+F(step+1,steps_per_phase))
            assert all(a.subset(b) for a,b in zip(box,domain)),('common domain',alpha,phase,step)
            proposed=[a+time*b for a,b in zip(y,rhs(box,alpha,phase))]
            assert all(a.subset(b) for a,b in zip(proposed,box)),('Picard',alpha,phase,step)
            c=coefficients(y,alpha,phase,order)
            tail=coefficients(box,alpha,phase,order+1)
            out=[]
            for j in range(len(y)):
                rem=tail[j][-1]*hp[order+1]
                out.append(sum((c[j][k]*hp[k] for k in range(order+1)),I(0))+rem)
                max_remainder=max(max_remainder,abs(F(rem.lo,S)),abs(F(rem.hi,S)))
            assert all(a.subset(b) for a,b in zip(out,box)),('Taylor endpoint',alpha,phase,step)
            y=out
            Mc,D,rod,dust,rho,pressure,memory,osc=energy(box,alpha,phase)
            pos=[box[0],box[1],box[7],Mc,1+F(1,5)*box[2],
                 F(1,5)+F(1,10)*box[2],memory,D+F(1,100)*box[2]]
            assert all(a.lo>0 for a in pos),('positive resources',alpha,phase,step)
            positivity=pos if positivity is None else [hull(a,b) for a,b in zip(positivity,pos)]
            con=y[1]**2-energy(y,alpha,phase)[4]
            assert con.contains_zero(),('constraint enclosure',alpha,phase,step,con.out())
            max_constraint=max(max_constraint,abs(F(con.lo,S)),abs(F(con.hi,S)))
            max_width=max(max_width,max(a.width() for a in y))
            all_tube=box if all_tube is None else [hull(a,b) for a,b in zip(all_tube,box)]
        assert y[6].lo==(phase+1)*S and y[6].hi==(phase+1)*S
        checkpoints[str(phase+1)]=endpoint_data(y,alpha,phase)
    # Retention uses the actual off-pulse Hamiltonian, not a reset or a new map.
    off=rhs(y,alpha,2)
    assert off[14].lo==off[14].hi==0
    assert off[7].lo==off[7].hi==0
    assert off[13].lo==y[14].lo and off[13].hi==y[14].hi
    data=dict(alpha=alpha,T='2',steps_per_phase=steps_per_phase,order=order,
              checkpoints=checkpoints,
              all_time_tube={n:a.out() for n,a in zip(NAMES,all_tube)},
              positive_all_time={n:a.out() for n,a in zip(
                  ['u','H','J','M_c','F','M_r','memory_complete','D_plus_Qp'],positivity)},
              max_endpoint_width=str(max_width),max_step_remainder=str(max_remainder),
              max_constraint_enclosure_abs=str(max_constraint),
              exact_constraint_uses_analytic_continuity=True,
              pulse_segments_separated_at_exact_theta_1=True,
              memory_K_off_pulse_derivative='0',memory_Y_off_pulse_derivative='K')
    return data,y


def source_checks():
    cases=[]
    for alpha in [0,1]:
        for phase in [0,1,2]:
            for fraction in [F(1,4),F(1,2),F(3,4)]:
                theta=phase+fraction
                y=[F(4,5),F(7,10),F(-1,50),F(-1,40),F(-1,1000),F(1,3),
                   theta,F(1,10),F(2,5),F(-1,5),F(0),
                   F(2,25),F(-1,50),F(-1,1000),F(-1,500)]
                u,H,p,v,z,light,theta,J,q,r,psi,X,V,Y,K=y
                dy=rhs(y,alpha,phase)
                up,Hp,pp,vp,zp,Lp,tp,Jp,qp,rp,psip,Xp,Vp,Yp,Kp=dy
                w,wp,f,fp=pulses(theta,phase)
                Mc,D,rod,dust,rho,pressure,memory,osc=energy(y,alpha,phase)
                Mdot=(Jp+q*qp+r*rp+K*Kp
                     +alpha*F(1,4)*(wp*q*p+w*qp*p+w*q*pp)
                     +F(1,10)*(fp*X*Y+f*Xp*Y+f*X*Yp))
                expected=alpha*F(1,4)*w*q*v+F(1,10)*f*Y*V
                assert Mdot==expected
                rod_dot=F(3,100)*v*u*u-2*H*rod
                dust_dot=(Mdot+F(1,100)*v)*u**3-3*H*dust
                scalar_dot=v*vp+p*v+V*Vp+X*V
                rhodot=rod_dot+dust_dot+scalar_dot
                assert rhodot+3*H*(rho+pressure)==0
                assert 2*H*Hp-rhodot==0
                # Keep the full source but omit the X reaction force.
                wrong=F(1,10)*f*Y*V*u**3
                if phase==1: assert wrong!=0
                # Keep the read interaction but omit the internal clock work.
                missing_control=F(1,10)*fp*X*Y*u**3
                if phase==1 and fraction!=F(1,2): assert missing_control!=0
                if phase==2:
                    assert Kp==0 and Yp==K and Jp==0
                cases.append(dict(alpha=alpha,phase=phase,theta=str(theta),
                    mass_exchange_rate=str(Mdot),continuity_residual='0',constraint_derivative='0',
                    missing_X_reaction_defect=str(wrong),
                    missing_clock_work_defect=str(missing_control)))
    y0=[F(1),F(0),F(0),F(0),F(0),F(0),F(0),F(1,10),F(1,2),F(0),F(0),
        F(1,10),F(0),F(0),F(0)]
    assert energy(y0,0,0)[4]==energy(y0,1,0)[4]==F(111,100)
    assert energy(y0,0,0)[0]==energy(y0,1,0)[0]==F(187,200)
    return dict(status='PASS',initial_rho_both='111/100',initial_Mc_both='187/200',
                finite_rational_diagnostics=cases,
                finite_points_do_not_replace_analytic_source_proof=True,
                shadow_equality_is_analytic_ODE_uniqueness_on_fixed_u_H=True,
                read_subsystem_has_no_direct_alpha_p_q_J_coupling=True,
                retained_quantity='K; Y is allowed to drift; no reset')


def rk_diagnostic(steps_per_phase=4000):
    import numpy as np
    y0=np.array([1,(111/100)**.5,0,0,0,0,0,.1,.5,0,0,.1,0,0,0],dtype=float)
    yy=np.concatenate([y0,y0]);dt=1/steps_per_phase;checkpoints={}
    for phase in [0,1]:
        def f(y):
            a,b=np.split(y,2)
            return np.concatenate([np.asarray(rhs(a,0,phase),dtype=float),
                                   np.asarray(rhs(b,1,phase),dtype=float)])
        for step in range(steps_per_phase):
            k1=f(yy);k2=f(yy+dt*k1/2);k3=f(yy+dt*k2/2);k4=f(yy+dt*k3)
            yy+=dt*(k1+2*k2+2*k3+k4)/6
            yy[6]=yy[21]=phase+(step+1)/steps_per_phase
        a,b=np.split(yy,2)
        checkpoints[str(phase+1)]=dict(alpha0=dict(zip(NAMES,map(float,a))),
             alpha1=dict(zip(NAMES,map(float,b))),delta_X=float(b[11]-a[11]),delta_K=float(b[14]-a[14]))
    return dict(method='binary64 piecewise RK4, diagnostic only',steps_per_phase=steps_per_phase,
                checkpoints=checkpoints,shadow_not_numerically_constructed=True)


def run():
    a,ya=integrate(0);b,yb=integrate(1)
    dK=yb[14]-ya[14]
    eta=F(1,10**10)  # Locked before any diagnostic: not fitted to the signal.
    absK=dK if dK.lo>0 else -dK
    margin=absK-2*eta
    assert not dK.contains_zero() and margin.lo>0,('fixed memory read contract failed',dK.out())
    xa=I(*a['checkpoints']['1']['state']['X'])
    xb=I(*b['checkpoints']['1']['state']['X'])
    return dict(round=1247,status='PASS',science_groups=2,
                arithmetic=dict(grid_denominator=str(S),outward_rounding='integer floor/ceil',
                   reused='1242 I/Jet/sqrt via frozen 1246; 1246 hull and initial-state helper',
                   certificate='per-side Picard tube and Taylor 8 + interval 9 remainder',
                   junctions='theta=1,2 exact; no high derivative across C2 pulse endpoints'),
                group_source_and_retention=source_checks(),
                group_finite_memory=dict(status='PASS',alpha0=a,alpha1=b,
                   delta_X_before_read=(xb-xa).out(),delta_K_final=dK.out(),
                   read_error_per_history_K=str(eta),margin_after_two_errors=margin.out(),
                   readout_accuracy_is_unimplemented_contract=True,
                   K_retains_after_read_while_Y_drifts=True),
                diagnostic=rk_diagnostic(),
                source_sha256={str(n):hashlib.sha256(p.read_bytes()).hexdigest() for n,p in
                   [('1242/check.py',TOOLS),('1246/check.py',PREVIOUS),('1247/check.py',Path(__file__))]},
                new_adopted_axioms=0,whole_GR_proved=False,
                limitation='Conditional classical Einstein branch; fixed-geometry deletion is not a second Einstein solution. No quantum instrument or empirical accuracy derived.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    parser.add_argument('--diagnose',action='store_true')
    args=parser.parse_args()
    if args.diagnose:
        assert not args.write
        print(json.dumps(rk_diagnostic(),ensure_ascii=False,indent=2))
    else:
        result=run();path=HERE/'results.json'
        if args.write:
            with path.open('x',encoding='utf-8') as f:
                json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
            mode='CREATED_ONCE'
        else:
            assert json.loads(path.read_text('utf-8'))==result,'saved result mismatch'
            mode='READ_ONLY_MATCH'
        print(json.dumps(dict(status='PASS',mode=mode,groups=result['science_groups'],
            delta_K=result['group_finite_memory']['delta_K_final'],
            margin=result['group_finite_memory']['margin_after_two_errors']),ensure_ascii=False))
