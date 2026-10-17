"""1246: one fixed autonomous pulse, complete sources and finite geometry.

Reuse the frozen 1242 integer outward I/Jet arithmetic.  Default replay is
read-only; --write creates results exactly once.  Floating point RK4 is a
separate diagnostic and is never used as an enclosure or an initial condition.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parent / '1242' / 'check.py'
spec = importlib.util.spec_from_file_location('round1242_interval_tools', TOOLS)
old = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(old)
I, Jet, S, sqrt_interval = old.I, old.Jet, old.S, old.sqrt_interval
NAMES = ['u', 'H', 'p', 'v', 'z', 'L', 'theta', 'J', 'q', 'r', 'psi']


def pulse(theta):
    a = theta * (1-theta)
    return 64*a**3, 192*a**2*(1-2*theta)


def energy(y, alpha):
    u, H, p, v, z, light, theta, J, q, r, psi = y
    w, wp = pulse(theta)
    Mc = F(7,10)+J+(q*q+r*r)*F(1,2)+alpha*F(1,4)*w*q*p
    D = Mc+F(1,50)
    rod = F(3,20)*(1+F(1,5)*p)*u*u
    dust = (D+F(1,100)*p)*u**3
    rho = (v*v+p*p)*F(1,2)+rod+dust
    pressure = (v*v-p*p)*F(1,2)-rod*F(1,3)
    return Mc, D, rod, dust, rho, pressure


def rhs(y, alpha):
    u, H, p, v, z, light, theta, J, q, r, psi = y
    w, wp = pulse(theta)
    Mc, D, rod, dust, rho, pressure = energy(y, alpha)
    coupling = alpha*F(1,4)*w
    return [-H*u, -F(3,2)*(v*v+F(2,3)*rod+dust), v,
            -3*H*v-p-F(3,100)*u*u-(F(1,100)+coupling*q)*u**3,
            F(1,10)*p, u, 1, -alpha*F(1,4)*wp*q*p,
            r, -q-coupling*p, F(1,4)*w*q*p]


def coefficients(y, alpha, order):
    c = [[I(a)] for a in y]
    for k in range(order):
        f = rhs([Jet(a) for a in c], alpha)
        for j in range(len(y)):
            term = f[j].c[k] if isinstance(f[j], Jet) else I(f[j] if k == 0 else 0)
            c[j].append(term.div(k+1))
    return c


def initial_interval():
    return [I(1), sqrt_interval(F(219,200)), I(0), I(0), I(0), I(0),
            I(0), I(F(1,10)), I(F(1,2)), I(0), I(0)]


def hull(a, b):
    return I(min(a.lo,b.lo), max(a.hi,b.hi), raw=True)


def integrate(alpha, steps=200, order=8):
    dt = F(1,steps)
    h, time = I(dt), I(0,dt)
    y = initial_interval()
    # These are a priori tube paddings, not values imported from RK4.
    radii = [F(2)*dt, F(3)*dt, F(1,5)*dt, F(1,2)*dt,
             F(1,30)*dt, F(11,10)*dt, F(0), F(1,4)*dt,
             F(3,5)*dt, F(3,5)*dt, F(1,20)*dt]
    domain = [I('.3','1.02'), I('.2','1.08'), I('-.15','.03'),
              I('-.15','.05'), I('-.02','.002'), I('-.01','1.02'),
              I(0,1), I('.04','.17'), I('.2','.52'), I('-.5','.03'),
              I('-.02','.01')]
    hp = [I(1)]
    for k in range(order+1):
        hp.append(hp[-1]*h)
    all_tube = None
    max_width, max_remainder = F(0), F(0)
    positivity = None
    max_constraint = F(0)
    for step in range(steps):
        box = [a+I(-r,r) for a,r in zip(y,radii)]
        # theta=t exactly; polynomial continuation outside [0,1] is not used.
        box[6] = I(F(step,steps), F(step+1,steps))
        assert all(a.subset(b) for a,b in zip(box,domain)), ('common domain',alpha,step)
        proposed = [a+time*b for a,b in zip(y,rhs(box,alpha))]
        assert all(a.subset(b) for a,b in zip(proposed,box)), ('Picard inclusion',alpha,step)
        c, tail = coefficients(y,alpha,order), coefficients(box,alpha,order+1)
        out = []
        for j in range(len(y)):
            rem = tail[j][-1]*hp[order+1]
            out.append(sum((c[j][k]*hp[k] for k in range(order+1)),I(0))+rem)
            max_remainder = max(max_remainder, abs(F(rem.lo,S)), abs(F(rem.hi,S)))
        assert all(a.subset(b) for a,b in zip(out,box)), ('Taylor endpoint',alpha,step)
        y = out
        Mc,D,rod,dust,rho,pressure = energy(box,alpha)
        pos = [box[7],Mc,1+F(1,5)*box[2],F(1,5)+F(1,10)*box[2],D+F(1,100)*box[2]]
        assert all(a.lo>0 for a in pos), ('positive resource domain',alpha,step)
        positivity = pos if positivity is None else [hull(a,b) for a,b in zip(positivity,pos)]
        con = y[1]*y[1]-energy(y,alpha)[4]
        assert con.contains_zero(), ('constraint zero inclusion',alpha,step,con.out())
        max_constraint = max(max_constraint,abs(F(con.lo,S)),abs(F(con.hi,S)))
        max_width = max(max_width,max(a.width() for a in y))
        all_tube = box if all_tube is None else [hull(a,b) for a,b in zip(all_tube,box)]
    assert y[6].lo==S and y[6].hi==S
    assert pulse(y[6])[0].lo==0 and pulse(y[6])[0].hi==0
    data = dict(alpha=alpha,T='1',steps=steps,order=order,
                endpoint={n:a.out() for n,a in zip(NAMES,y)},
                Brod=(F(1,10)*y[0]**2).out(),
                all_time_tube={n:a.out() for n,a in zip(NAMES,all_tube)},
                positive_all_time={n:a.out() for n,a in zip(['J','M_c','F','M_r','D_plus_Qp'],positivity)},
                max_endpoint_width=str(max_width),max_step_remainder=str(max_remainder),
                max_constraint_enclosure_abs=str(max_constraint),
                constraint_exactness_uses_analytic_continuity_identity=True,
                pulse_finished_exact=True)
    return data,y


def source_checks():
    cases=[]
    # Exact rational point diagnostics check the implementation of the proof's
    # universal polynomial identities; the finite points do not prove them.
    for alpha in [0,1]:
        for theta in [F(1,4),F(1,2),F(3,4)]:
            y=[F(4,5),F(7,10),F(-1,50),F(-1,40),F(-1,1000),F(1,3),
               theta,F(1,10),F(2,5),F(-1,5),F(0)]
            u,H,p,v,z,light,theta,J,q,r,psi=y
            dy=rhs(y,alpha);up,Hp,pp,vp,zp,Lp,tp,Jp,qp,rp,psip=dy
            w,wp=pulse(theta)
            Mc,D,rod,dust,rho,pressure=energy(y,alpha)
            Mdot=Jp+q*qp+r*rp+alpha*F(1,4)*(wp*tp*q*p+w*qp*p+w*q*pp)
            expected=alpha*F(1,4)*w*q*v
            assert Mdot==expected
            rod_dot=F(3,100)*v*u*u-2*H*rod
            dust_dot=(Mdot+F(1,100)*v)*u**3-3*H*dust
            scalar_dot=v*vp+p*v
            rhodot=rod_dot+dust_dot+scalar_dot
            assert rhodot+3*H*(rho+pressure)==0
            assert 2*H*Hp-rhodot==0
            # Delete only pulse force in p, retaining the full matter source.
            wrong_scalar=scalar_dot+alpha*F(1,4)*w*q*v*u**3
            wrong=rod_dot+dust_dot+wrong_scalar+3*H*(rho+pressure)
            assert wrong==alpha*F(1,4)*w*q*v*u**3
            if alpha: assert wrong!=0
            assert psip==F(1,4)*w*q*p  # B, not alpha, is the canonical program.
            cases.append(dict(alpha=alpha,theta=str(theta),mass_rate=str(Mdot),
                              continuity_residual='0',constraint_derivative='0',
                              missing_pulse_force_defect=str(wrong),program_phase_rate=str(psip)))
    y0=[F(1),F(0),F(0),F(0),F(0),F(0),F(0),F(1,10),F(1,2),F(0),F(0)]
    assert energy(y0,0)[4]==energy(y0,1)[4]==F(219,200)
    return dict(status='PASS',initial_rho_both='219/200',initial_Mc_both='37/40',
                finite_rational_diagnostics=cases,
                universal_identity_proof_is_in_proof_md=True,
                program_B_is_conserved_internal_initial_data=True)


def rk_diagnostic(steps=4000):
    import numpy as np
    # Independent from enclosures; stable ordinary RK4 and a joint shadow state.
    y0=np.array([1,(219/200)**.5,0,0,0,0,0,.1,.5,0,0],dtype=float)
    y=np.concatenate([y0,y0,y0])
    def f(yy):
        a,b,c=np.split(yy,3)
        da=np.asarray(rhs(a,0),dtype=float)
        db=np.asarray(rhs(b,1),dtype=float)
        sc=c.copy();sc[0]=a[0];sc[1]=a[1]
        dc=np.asarray(rhs(sc,1),dtype=float)
        # Prescribed geometry exactly copies alpha=0. No claim of Einstein.
        dc[0]=da[0];dc[1]=da[1];dc[5]=da[5]
        return np.concatenate([da,db,dc])
    dt=1/steps
    for _ in range(steps):
        k1=f(y);k2=f(y+dt*k1/2);k3=f(y+dt*k2/2);k4=f(y+dt*k3)
        y+=dt*(k1+2*k2+2*k3+k4)/6
    a,b,c=np.split(y,3)
    assert a[0]==c[0] and a[1]==c[1] and a[5]==c[5]
    return dict(method='ordinary binary64 RK4; diagnostic only',steps=steps,
                alpha0=dict(zip(NAMES,map(float,a))),alpha1=dict(zip(NAMES,map(float,b))),
                shadow_alpha1=dict(zip(NAMES,map(float,c))),
                true_delta_B=float(.1*(b[0]**2-a[0]**2)),true_delta_L=float(b[5]-a[5]),
                shadow_delta_z=float(c[4]-a[4]),shadow_delta_B='0 exactly by construction',
                shadow_delta_L='0 exactly by construction',
                shadow_is_not_second_Einstein_solution=True,
                shadow_z_nonzero_is_float_diagnostic_not_certified=True)


def run():
    a,ya=integrate(0);b,yb=integrate(1)
    dB=F(1,10)*(yb[0]**2-ya[0]**2)
    dL=yb[5]-ya[5]
    dz=yb[4]-ya[4]
    assert not dB.contains_zero() and not dL.contains_zero()
    absB=dB if dB.lo>0 else -dB
    absL=dL if dL.lo>0 else -dL
    # These numbers will be fixed after the diagnostic, before the first save.
    eta_B=F(1,10**8);eta_L=F(1,10**8)
    marginB=absB-2*eta_B;marginL=absL-2*eta_L
    assert marginB.lo>0 and marginL.lo>0
    return dict(round=1246,status='PASS',science_groups=2,
                arithmetic=dict(grid_denominator=str(S),outward_rounding='integer floor/ceil',
                   inherited_tools='1242 I, Jet, sqrt_interval; no scipy or new interval library',
                   time='theta=t exactly; polynomial pulse only on 0<=theta<=1',
                   certificate='Picard inclusion and Taylor order 8 + interval order 9 remainder'),
                group_complete_sources=source_checks(),
                group_finite_geometry=dict(status='PASS',alpha0=a,alpha1=b,
                   delta_Brod=dB.out(),delta_L=dL.out(),delta_z=dz.out(),
                   read_error_per_history_B=str(eta_B),read_error_per_history_L=str(eta_L),
                   margin_B=marginB.out(),margin_L=marginL.out(),
                   readout_budgets_are_unimplemented_contracts=True),
                shadow_diagnostic=rk_diagnostic(),
                source_sha256={'1242/check.py':hashlib.sha256(TOOLS.read_bytes()).hexdigest(),
                               '1246/check.py':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                new_adopted_axioms=0,whole_GR_proved=False,
                limitation='Conditional Einstein effective model and classical internal program. No quantum preparation/readout or source-to-geometry selection derived.')


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
                 delta_Brod=result['group_finite_geometry']['delta_Brod'],
                 delta_L=result['group_finite_geometry']['delta_L']),ensure_ascii=False))
