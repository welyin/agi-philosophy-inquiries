"""1250 independent null-sheet audit; symbolic identities and Decimal diagnostic.

No author checker is imported. The numerical part starts at the frozen t=2
endpoint, uses a rather than u, and the canonical geometric H equation.
It is a convergence diagnostic, not a substitute for the author's interval tube.
Default replay is read-only; --write creates the first saved result exclusively.
"""
from __future__ import annotations
import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
OUT=HERE/'independent_results.json'
ETA=F(1,10**9)
NAMES=('a','H','p','v','X','V','Mc','y','x','us','Hs','mk')
INDEX={n:i for i,n in enumerate(NAMES)}
ZERO=(0,)*len(NAMES)


class P:
    def __init__(self,x=0):
        if isinstance(x,P):self.c=x.c.copy()
        elif isinstance(x,dict):self.c={k:F(v) for k,v in x.items() if v}
        else:self.c={ZERO:F(x)} if x else {}
    @staticmethod
    def var(n):
        e=list(ZERO);e[INDEX[n]]=1
        return P({tuple(e):F(1)})
    def __add__(self,b):
        b=P(b);out=self.c.copy()
        for k,v in b.c.items():out[k]=out.get(k,0)+v
        return P(out)
    __radd__=__add__
    def __neg__(self):return P({k:-v for k,v in self.c.items()})
    def __sub__(self,b):return self+-P(b)
    def __rsub__(self,b):return P(b)+-self
    def __mul__(self,b):
        b=P(b);out={}
        for k,v in self.c.items():
            for l,w in b.c.items():
                e=tuple(i+j for i,j in zip(k,l));out[e]=out.get(e,0)+v*w
        return P(out)
    __rmul__=__mul__
    def __truediv__(self,b):return self*F(1,b)
    def __pow__(self,n):
        if n<0:
            assert len(self.c)==1
            e,v=next(iter(self.c.items()))
            return P({tuple(i*n for i in e):v**n})
        out=P(1)
        for _ in range(n):out=out*self
        return out
    def d(self,n):
        i=INDEX[n];out={}
        for e,c in self.c.items():
            if e[i]:
                q=list(e);q[i]-=1;out[tuple(q)]=c*e[i]
        return P(out)
    def __bool__(self):return bool(self.c)


def symbolic_group():
    a,H,p,v,X,V,Mc,y,x,us,Hs,mk=[P.var(n) for n in NAMES]
    Frod=1+p/5
    rho=(v*v+V*V+p*p+X*X)/2+F(3,20)*Frod*a**-2+(Mc+F(1,50)+p/100)*a**-3
    pressure=(v*v+V*V-p*p-X*X)/2-F(1,20)*Frod*a**-2
    forward={'a':a*H,'H':-F(3,2)*(H*H+pressure),'p':v,
             'v':-3*H*v-p-F(3,100)*a**-2-F(1,100)*a**-3,
             'X':V,'V':-3*H*V-X}
    derivative=lambda f,flow:sum((f.d(n)*z for n,z in flow.items()),P())
    C=H*H-rho
    reverse={n:-z for n,z in forward.items()}
    reverse.update(y=Hs-H*y,x=us*a)
    kt=a**-1*us**-1
    theta=2*kt*(H-Hs*y**-1)
    Rgeom=3*(H*H+pressure)*kt*kt
    Tnull=(rho+pressure)*kt*kt
    return _finish_symbolic(a,H,rho,pressure,C,forward,reverse,derivative,
                            y,us,kt,theta,Rgeom,Tnull,mk)


def _finish_symbolic(a,H,rho,pressure,C,forward,reverse,derivative,y,us,kt,theta,Rgeom,Tnull,mk):
    I0prime=Rgeom*y*y*a*us
    Z0prime=theta*theta*y*y*a*us/2
    checks={
        'full_source_continuity':derivative(rho,forward)+3*H*(rho+pressure),
        'canonical_constraint_transport':derivative(C,forward)+3*H*C,
        'independent_curvature_source_difference':Rgeom-3*Tnull-3*C*kt*kt,
        'affine_area_velocity':-kt*derivative(y*y,reverse)-theta*y*y,
        'Raychaudhuri_shear_free':-kt*derivative(theta,reverse)+theta*theta/2+Rgeom,
        'area_derivative_integrated_kernel':derivative(y*y*theta,reverse)-I0prime+Z0prime,
        'memory_null_omission':3*(Tnull-(Tnull-mk*a**-3*kt*kt))-3*mk*a**-3*kt*kt,
    }
    for n,z in checks.items():assert not z,(n,z.c)
    # Whole sphere remains inside one injectivity chart of the original T3.
    assert F(93,100)/27>F(9,50)**2
    comoving_radius_upper=F(3,8)/F(9,50)+F(3,8)/64
    assert comoving_radius_upper<5
    ellmin=F(1,72);ellmax=F(1,64)
    area_ratio_min=F(64,81)
    deletion_lower=3*F(1,100)*F(1,27)*area_ratio_min*ellmin**2/2
    assert deletion_lower>7*ETA
    return {'status':'PASS','method':'Exact Laurent-polynomial coefficient identities; geometric H equation and material source computed separately.',
            'zero_identities':list(checks),
            'finite_area_identity':'DeltaA/Astar = ell*I0-I1-ell*Z0+Z1 (theta0=0)',
            'raw_Clausius_geometric_remainder':'ell*I0-2*I1-ell*Z0+Z1',
            'independent_positive_bounds':{'Hstar_lower':'9/50','comoving_radius_upper':str(comoving_radius_upper),
                'T3_injectivity_radius_comoving':5,'affine_length_lower':str(ellmin),'affine_length_upper':str(ellmax),
                'area_ratio_lower':str(area_ratio_min),'normalized_memory_deletion_lower':str(deletion_lower),
                'independent_corrected_read_margin_lower':str(deletion_lower-7*ETA)},
            'Wald_normalization':'8*pi*G=3, s*hbar=2*pi/3, T_U*s=kappa_B/3; temperature reference remains conditional.'}


def midpoint(pair):return (D(pair[0])+D(pair[1]))/2


def source(y,Mc):
    a,H,p,v,X,V=y[:6]
    rho=(v*v+V*V+p*p+X*X)/2+D('.15')*(1+p/5)/a**2+(Mc+D('.02')+p/100)/a**3
    pressure=(v*v+V*V-p*p-X*X)/2-D('.05')*(1+p/5)/a**2
    return rho,pressure


def background(y,Mc):
    a,H,p,v,X,V=y[:6]
    _,pressure=source(y,Mc)
    return [a*H,-D('1.5')*(H*H+pressure),v,
            -3*H*v-p-D('.03')/a**2-D('.01')/a**3,V,-3*H*V-X]


def rk4(initial,span,steps,rhs):
    y=list(initial);dt=span/steps
    for _ in range(steps):
        k1=rhs(y)
        k2=rhs([a+dt*b/2 for a,b in zip(y,k1)])
        k3=rhs([a+dt*b/2 for a,b in zip(y,k2)])
        k4=rhs([a+dt*b for a,b in zip(y,k3)])
        y=[a+dt*(b+2*c+2*d+e)/6 for a,b,c,d,e in zip(y,k1,k2,k3,k4)]
    return y


def numeric_branch(branch,steps):
    endpoint=branch['checkpoints']['2'];st=endpoint['state']
    assert F(endpoint['energy']['M_c'][0])>F(93,100)
    Mc=midpoint(endpoint['energy']['M_c']);mk=midpoint(endpoint['energy']['memory_complete'])
    initial=[1/midpoint(st['u'])]+[midpoint(st[n]) for n in ('H','p','v','X','V')]
    end=rk4(initial,D('.25'),steps,lambda z:background(z,Mc))
    astar,Hstar=end[:2];us=1/astar
    # y=R/Rstar,x=-lambda,Ig0,Ig1,Z0,Z1,It1,Dmem.
    start=end+[D(1)]+[D(0)]*7
    def rhs(z):
        a,H=z[:2];ys,x,Ig0,Ig1,Z0,Z1,It1,Dmem=z[6:]
        rho,pr=source(z,Mc);kt=1/(a*us)
        G0=3*(H*H+pr)*ys*ys*kt
        T0=3*(rho+pr)*ys*ys*kt
        Z=2*kt*(H*ys-Hstar)**2
        return [-v for v in background(z,Mc)]+[Hstar-H*ys,us*a,G0,x*G0,Z,x*Z,x*T0,3*x*ys*ys*mk*kt/a**3]
    final=rk4(start,D(1)/64,max(32,steps//8),rhs)
    a,H=final[:2];ys,x,Ig0,Ig1,Z0,Z1,It1,Dmem=final[6:]
    delta=1-ys*ys
    exact=x*Ig0-Ig1-x*Z0+Z1
    correction=x*Ig0-2*Ig1-x*Z0+Z1
    rho,_=source(final,Mc)
    return {'a_star':astar,'H_star':Hstar,'affine_length':x,'area_change':delta,
            'geometric_boost_integral':Ig1,'material_boost_integral':It1,
            'two_weight_geometric_correction':correction,
            'raw_material_balance':delta-It1,'omitted_memory_shift':Dmem,
            'corrected_material_balance':delta-It1-correction,
            'exact_area_identity_residual':delta-exact,
            'terminal_constraint':end[1]**2-source(end,Mc)[0],
            'past_constraint':H*H-rho}


def numeric_group(compare_author=False):
    old=json.loads((HERE.parent/'1247/results.json').read_text(encoding='utf-8'))
    rows={}
    with localcontext() as ctx:
        ctx.prec=48
        for label in ('alpha0','alpha1'):
            coarse=numeric_branch(old['group_finite_memory'][label],1024)
            fine=numeric_branch(old['group_finite_memory'][label],2048)
            changes={k:abs(fine[k]-coarse[k]) for k in fine}
            assert max(changes.values())<D('1e-12')
            assert abs(fine['exact_area_identity_residual'])<D('1e-15')
            assert abs(fine['corrected_material_balance'])<D('1e-15')
            assert fine['area_change']>0 and fine['omitted_memory_shift']>D('8e-8')
            rows[label]={'fine':{k:str(v) for k,v in fine.items()},
                         'coarse_fine_absolute_change':{k:str(v) for k,v in changes.items()}}
    crosscheck={}
    if compare_author:
        saved=json.loads((HERE/'results.json').read_text(encoding='utf-8'))
        balance=saved['group_finite_balance'];gamma=F(balance['raw_common_acceptance'])
        assert not balance['raw_deletion_certified'] and balance['corrected_deletion_certified']
        assert balance['omitted_raw_entire_error_interval_inside_acceptance']
        pairs={'affine_length':'affine_length','area_change':'geometric_area_delta',
               'material_boost_integral':'boost_full','two_weight_geometric_correction':'geometry_correction',
               'raw_material_balance':'raw_balance','omitted_memory_shift':'positive_memory_gap'}
        for i,branch in enumerate(saved['group_same_history_null']['branches']):
            label=f'alpha{i}';screen=branch['screen']
            lo,hi=map(F,screen['raw_omitted_balance'])
            assert -gamma<lo-3*ETA<hi+3*ETA<gamma
            dl,dh=map(F,screen['positive_memory_gap'])
            assert dl>F(1,11809800)>7*ETA
            assert F(branch['read_contract']['full_corrected_error'])==3*ETA
            assert F(branch['read_contract']['omitted_corrected_error'])==4*ETA
            distance={}
            for own,their in pairs.items():
                lo,hi=map(F,screen[their]);point=F(rows[label]['fine'][own])
                gap=max(lo-point,point-hi,F(0));distance[own]=str(gap)
                # A diagnostic comparison only; RK step-size agreement is not
                # converted into a rigorous enclosure of the exact trajectory.
                assert gap<F(1,10**16)
            crosscheck[label]={'RK_point_distance_to_author_interval':distance,
                              'omitted_raw_entire_read_error_inside_acceptance':True,
                              'corrected_gap_exceeds_independent_lower':True,
                              'strict_author_gap':list(map(str,(dl,dh)))}
    return {'status':'PASS' if compare_author else 'DIAGNOSTIC_PASS',
            'method':'48-digit Decimal RK4, a/H canonical geometric equation, 1024/2048 hold steps and separate null integration; exact rational error-contract checks.',
            'scope':'The ODE calculation is a diagnostic, not a rigorous enclosure or an assertion that the midpoint is the original exact state. The rational budget comparisons are exact.',
            'branches':rows,'author_interval_and_budget_crosscheck':crosscheck}


def run(draft=False):
    result={'round':1250,'status':'DRAFT_NOT_REVIEWED' if draft else 'PASS','scientific_groups':2,
            'symbolic_group':symbolic_group(),'numeric_group':numeric_group(compare_author=not draft)}
    if not draft:
        paths={
            'independent_check.py':Path(__file__),
            **{name:HERE/name for name in ('admission.md','proof.md','check.py','results.json','NEXT.md')},
            '../research_note_1250.md':HERE.parent/'research_note_1250.md',
            '../1247/proof.md':HERE.parent/'1247/proof.md',
            '../1247/results.json':HERE.parent/'1247/results.json',
            '../1248/proof.md':HERE.parent/'1248/proof.md',
            '../1242/check.py':HERE.parent/'1242/check.py',
        }
        result['bindings']={name:hashlib.sha256(path.read_bytes()).hexdigest() for name,path in paths.items()}
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--draft',action='store_true');parser.add_argument('--write',action='store_true');args=parser.parse_args()
    assert not(args.draft and args.write)
    result=run(draft=args.draft)
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
        print('FIRST_EXCLUSIVE_WRITE PASS; scientific_groups=2')
    elif not args.draft:
        assert json.loads(OUT.read_text(encoding='utf-8'))==result,'READ_ONLY_MISMATCH'
        print('READ_ONLY_MATCH PASS; scientific_groups=2')
    print(json.dumps({'status':result['status'],'symbolic_identities':len(result['symbolic_group']['zero_identities']),
                      'deletion_lower':result['symbolic_group']['independent_positive_bounds']['normalized_memory_deletion_lower'],
                      'raw_deletion':False,'corrected_deletion':True},ensure_ascii=False))
