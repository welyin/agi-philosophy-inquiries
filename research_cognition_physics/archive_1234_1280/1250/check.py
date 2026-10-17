"""1250: the fixed 1247 hold history, null area, complete boost and omission.

Exact integer outward intervals; no new integration of the frozen [0,2] run.
Default verifies saved results without writing. --write creates results once.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNTIME = HERE.parent/'1242'/'check.py'
PREVIOUS = HERE.parent/'1247'/'results.json'
EXPECTED = {
    '1242/check.py':'2a25be6a9b371c0689f03458c0151dcd534339d983989ca8c914056749aabc16',
    '1247/results.json':'dbcd590e7816bc423c49a4e1bf83bab6fa50a84dd5578249cc493a15cd72bf06',
}
for name, digest in EXPECTED.items():
    assert hashlib.sha256((HERE.parent/name).read_bytes()).hexdigest()==digest
spec=importlib.util.spec_from_file_location('round1242_interval_1250',RUNTIME)
old=importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(old)
I,Jet,S=old.I,old.Jet,old.S
BG=['u','H','p','v','X','V','a']
NULL=['x','y','G0','G1','Z0','Z1','B','Dmiss']
ETA=F(1,10**9)
WIDTH=F(1,64)
TSTAR=F(9,4)


def reciprocal(x):
    assert x.lo>0
    return I(S*S//x.hi,-((-S*S)//x.lo),raw=True)


def hull(a,b):
    return I(min(a.lo,b.lo),max(a.hi,b.hi),raw=True)


def energy(z,mc):
    u,H,p,v,X,V,a=z[:7]
    rod=F(3,20)*(1+F(1,5)*p)*u*u
    dust=(F(1,100)*p+mc+F(1,50))*u**3
    rho=(v*v+p*p+V*V+X*X)*F(1,2)+rod+dust
    pressure=(v*v-p*p+V*V-X*X)*F(1,2)-rod*F(1,3)
    rhoP=v*v+V*V+rod*F(2,3)+dust
    return rho,pressure,rhoP,rod,dust


def bg_rhs(z,mc):
    u,H,p,v,X,V,a=z[:7]
    # The original full-material evolution, including the original clock J.
    rho,pressure,rhoP,rod,dust=energy(z,mc)
    return [-H*u,-F(3,2)*rhoP,v,
            -3*H*v-p-F(3,100)*u*u-F(1,100)*u**3,
            V,-3*H*V-X,H*a]


def screen_rhs(z,mc,mk,us,hs,ius):
    u,H,p,v,X,V,a=z[:7]
    x,y,g0,g1,z0,z1,b,d=z[7:]
    forward=bg_rhs(z,mc)
    rho,pressure,rhoP,rod,dust=energy(z,mc)
    kt=u*ius
    # Ricci from the FLRW metric derivative: Rkk=-2 Hdot (kt)^2.
    # No definition of geometric correction as area minus matter heat.
    geometric=-2*forward[1]*y*y*kt
    thermal=3*rhoP*y*y*kt
    expansion=2*kt*(H*y-hs)**2
    omission=3*x*y*y*mk*u**3*kt
    return [-f for f in forward]+[a*us,-H*y+hs,
        geometric,x*geometric,expansion,x*expansion,x*thermal,omission]


def coefficients(z,fun,order):
    c=[[I(a)] for a in z]
    for k in range(order):
        d=fun([Jet(a) for a in c])
        for j in range(len(z)):
            term=d[j].c[k] if isinstance(d[j],Jet) else I(d[j] if k==0 else 0)
            c[j].append(term.div(k+1))
    return c


def integrate(z,fun,duration,steps,speed_bounds,order=8):
    dt=duration/steps
    h,dtbox=I(dt),I(0,dt)
    hp=[I(1)]
    for _ in range(order+1):hp.append(hp[-1]*h)
    maxrem=F(0);maxwidth=F(0);allbox=None
    for step in range(steps):
        box=[a+I(-dt*v,dt*v) for a,v in zip(z,speed_bounds)]
        picard=[a+dtbox*f for a,f in zip(z,fun(box))]
        assert all(a.subset(b) for a,b in zip(picard,box)),('Picard inclusion',step)
        c=coefficients(z,fun,order)
        tail=coefficients(box,fun,order+1)
        new=[]
        for j in range(len(z)):
            rem=tail[j][-1]*hp[-1]
            new.append(sum((c[j][k]*hp[k] for k in range(order+1)),I(0))+rem)
            maxrem=max(maxrem,abs(F(rem.lo,S)),abs(F(rem.hi,S)))
        assert all(a.subset(b) for a,b in zip(new,box)),('endpoint tube',step)
        z=new
        allbox=box if allbox is None else [hull(a,b) for a,b in zip(allbox,box)]
        maxwidth=max(maxwidth,max(a.width() for a in z))
    return z,dict(steps=steps,duration=str(duration),order=order,remainder_order=order+1,
        picard_tubes=steps,all_time_tube=[a.out() for a in allbox],
        max_endpoint_width=str(maxwidth),max_step_remainder=str(maxrem)),allbox


def branch(alpha,previous):
    start=previous['group_finite_memory'][f'alpha{alpha}']['checkpoints']['2']
    state=start['state'];eng=start['energy']
    z=[I(*state[n]) for n in BG[:-1]]
    z.append(reciprocal(z[0]))
    mc=I(*eng['M_c'])
    # Same physical clock source and J; omit only this stored memory component.
    mk=F(1,100)+I(*state['K'])**2*F(1,2)
    assert mk.lo>I(F(1,100)).lo and mc.lo>I(F(9,10)).lo
    speeds=[F(1,5),F(1,5),F(1,10),F(1,10),F(1,10),F(1,10),F(1)]
    zs,hold,holdbox=integrate(z,lambda zz:bg_rhs(zz,mc),F(1,4),64,speeds)
    us,hs=zs[:2];ius=reciprocal(us)
    rho,pressure,rhoP,_,_=energy(zs,mc)
    assert (hs*hs-rho).contains_zero()
    assert (zs[0]*zs[6]-1).contains_zero()
    assert holdbox[0].lo>0 and holdbox[1].lo>0
    assert (1+holdbox[2]*F(1,5)).lo>0
    assert (F(1,5)+holdbox[2]*F(1,10)).lo>0
    nz=zs+[I(0),I(1)]+[I(0)]*6
    nspeeds=speeds+[F(2),F(1,20),F(1),F(1,20),F(1,20),F(1,100),F(1,20),F(1,100)]
    fun=lambda zz:screen_rhs(zz,mc,mk,us,hs,ius)
    end,null,box=integrate(nz,fun,WIDTH,16,nspeeds)
    ell,y,g0,g1,z0,z1,b,d=end[7:]
    area=1-y*y
    area_geometric=ell*g0-g1-ell*z0+z1
    optical_residual=area-area_geometric
    ew=ell*g0-2*g1
    ex=-ell*z0+z1
    correction=ew+ex
    corrected=area-b-correction
    omitted=area-(b-d)-correction
    raw=area-b
    raw_omitted=area-(b-d)
    assert optical_residual.contains_zero()
    assert corrected.contains_zero()
    assert (g1-b).contains_zero()
    assert box[8].lo>0 and ell.lo>0 and d.lo>0
    assert all(box[j].lo>0 for j in (0,1,6))
    assert ex.hi<0 and area.lo>0
    # Whole sphere fits the original L=10 torus; a cap requires less.
    comoving_radius=box[8]*reciprocal(hs)*box[0]
    assert comoving_radius.hi<5*S
    assert omitted.lo>7*I(ETA).hi
    assert omitted.width()<F(1,10**12)
    # Menu observation errors are extra permissions, not numerical errors.
    return dict(status='PASS',alpha=alpha,
        frozen_t2=dict(M_c=mc.out(),memory_omitted=mk.out(),J=state['J']),
        hold=hold,terminal_geometry={n:a.out() for n,a in zip(BG,zs)},
        terminal_constraint=(hs*hs-rho).out(),terminal_ua=(us*zs[6]-1).out(),
        null_tube=null,null_names=BG+NULL,
        screen=dict(affine_length=ell.out(),area_ratio_past=(y*y).out(),
          geometric_area_delta=area.out(),area_from_metric_optics=area_geometric.out(),
          optical_residual=optical_residual.out(),boost_full=b.out(),
          weight_correction=ew.out(),expansion_correction=ex.out(),
          geometry_correction=correction.out(),curvature_boost_minus_source=(g1-b).out(),
          raw_balance=raw.out(),raw_omitted_balance=raw_omitted.out(),
          corrected_balance=corrected.out(),corrected_omitted_balance=omitted.out(),
          positive_memory_gap=d.out(),comoving_radius_tube=comoving_radius.out(),
          positive_area_certified=True,terminal_expansion_exact='0',terminal_kt_exact='1'),
        read_contract=dict(error_per_normalized_menu=str(ETA),
          full_corrected_error=str(3*ETA),omitted_corrected_error=str(4*ETA),
          corrected_separation_margin=(d-7*ETA).out(),physical_instrument_not_derived=True))


def rational_diagnostics():
    out=[]
    for u,H,p,v,X,V,mc in [tuple(map(F,['1/3','1/4','-1/100','1/1000','1/40','-1/25','93/100'])),
                            tuple(map(F,['7/20','13/50','1/200','-1/300','1/30','1/50','94/100']))]:
        a=1/u;z=[u,H,p,v,X,V,a]
        du,dH,dp,dv,dX,dV,da=bg_rhs(z,mc)
        rho,P,rhoP,rod,dust=energy(z,mc)
        rod_dot=F(3,100)*v*u*u-2*H*rod
        dust_dot=F(1,100)*v*u**3-3*H*dust
        rho_dot=v*dv+p*v+V*dV+X*V+rod_dot+dust_dot
        continuity=rho_dot+3*H*rhoP
        constraint_dot=2*H*dH-rho_dot
        # At a general positive screen point use independent metric derivative.
        y=F(999,1000);hs=F(1,4);us=F(1,3);kt=u/us
        # Future-time ydot=Hy-hs, theta=2kt(H-hs/y).
        dydt=H*y-hs;ktdot=-H*kt
        theta=2*kt*(H-hs/y)
        theta_dot=2*ktdot*(H-hs/y)+2*kt*(dH+hs*dydt/(y*y))
        ray=kt*theta_dot+theta*theta/2-2*dH*kt*kt
        assert continuity==constraint_dot==ray==0
        out.append(dict(continuity=str(continuity),constraint_derivative=str(constraint_dot),
                        affine_raychaudhuri=str(ray)))
    return out


def calculate():
    previous=json.loads(PREVIOUS.read_text('utf-8'))
    branches=[branch(i,previous) for i in (0,1)]
    ec=[I(*b['screen']['geometry_correction']) for b in branches]
    raw_om=[I(*b['screen']['raw_omitted_balance']) for b in branches]
    maxrem=max(abs(q) for e in ec for q in (F(e.lo,S),F(e.hi,S)))
    gamma=maxrem+2*ETA
    # Conservative common raw acceptance |area-boost| <= max |Egeom|+2 eta.
    # The omitted heat uses full heat minus the separately bounded omission.
    raw_separated=all((r-3*ETA).lo>I(gamma).hi for r in raw_om)
    raw_om_inside=all(max(abs(F(r.lo,S)),abs(F(r.hi,S)))+3*ETA<gamma for r in raw_om)
    assert raw_om_inside and not raw_separated,'fixed raw budget must retain its certified failure'
    return dict(round=1250,status='PASS',science_groups=2,
        fixed=dict(tstar=str(TSTAR),time_width=str(WIDTH),error_per_menu=str(ETA),
                   cap_solid_angle='1',Wald_s_times_hbar='2*pi/3',eight_pi_G='3'),
        arithmetic=dict(runtime='frozen 1242 integer I/Jet',grid_denominator=str(S),
          method='Picard inclusion + Taylor order8 and interval order9 remainder',
          initial_data='frozen 1247 t=2 endpoint intervals; no [0,2] replay'),
        group_same_history_null=dict(status='PASS',rational_diagnostics=rational_diagnostics(),
                                    branches=branches),
        group_finite_balance=dict(status='PASS',
          raw_common_remainder_bound=str(maxrem),raw_common_acceptance=str(gamma),
          raw_deletion_certified=raw_separated,corrected_deletion_certified=True,
          omitted_raw_entire_error_interval_inside_acceptance=raw_om_inside,
          raw_claim=('SEPARATED_WITH_DECLARED_REMAINDER_AND_READ_BUDGET' if raw_separated else
                     'NOT_SEPARATED_BY_THIS_RAW_BUDGET_USE_CORRECTED_BALANCE'),
          corrected_gap_intervals=[b['screen']['positive_memory_gap'] for b in branches],
          corrected_margin_intervals=[b['read_contract']['corrected_separation_margin'] for b in branches]),
        source_sha256={n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in
                       [('1242/check.py',RUNTIME),('1247/results.json',PREVIOUS),('1250/check.py',Path(__file__))]},
        new_adopted_axioms=0,CT_adopted=False,actual_Unruh_reference_delivered=False,
        Wald_equals_bare_region_vN_proved=False,whole_GR_proved=False,
        scope='Conditional EH effective history, local null cap and declared normalized read contract. Omitted source is an incomplete ledger on the same geometry, not a new Einstein solution.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    result=calculate();path=HERE/'results.json'
    if args.write:
        with path.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        assert json.loads(path.read_text('utf-8'))==result,'saved result mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status='PASS',mode=mode,groups=result['science_groups'],
        balance=result['group_finite_balance']),ensure_ascii=False))
