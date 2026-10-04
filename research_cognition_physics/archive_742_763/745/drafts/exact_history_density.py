"""Exact rational Gaussian-integer jets of the inherited native onsite H.

Frozen matching-table decimals and diagnostic Yukawa decimals are interpreted
as exact rational model inputs. The linear solve for u is exact. No coefficient
threshold or floating polynomial operation is used. Quadrature is not certified.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import lcm,comb
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'round744_drafts'))
import history_jet_probe as old
ZERO=(0,)*5

def parameters():
    row=json.loads((ROOT/'joint_singlet_common_mass_rg_results.json').read_text('utf8'),
                   parse_float=Q)['examples'][2]['state']
    L=[[row['lambda_H'],row['p']],[row['p'],row['lambda_s']]]
    C=[row['x']/4,row['y']/4]
    determinant=L[0][0]*L[1][1]-L[0][1]*L[1][0]
    u=[(C[0]*L[1][1]-C[1]*L[0][1])/determinant,
       (L[0][0]*C[1]-L[1][0]*C[0])/determinant]
    P=[[Q(int(i==j))-u[i]/12 for j in range(2)] for i in range(2)]
    quad=[[sum(P[k][i]*L[k][ell]*P[ell][j] for k in range(2) for ell in range(2))/4
           for j in range(2)] for i in range(2)]
    lin=[-sum(P[k][i]*L[k][ell]*u[ell] for k in range(2) for ell in range(2))/4
         for i in range(2)]
    const=sum(u[k]*L[k][ell]*u[ell] for k in range(2) for ell in range(2))/16
    return dict(L=L,u=u,constant=const,linear=lin,quadratic=quad)

PAR=parameters()
# Kinetic denominators 12, Yukawa diagnostic decimal denominators 100.
DEN=lcm(12,100,PAR['constant'].denominator,
        *(p.denominator for p in PAR['linear']),
        *(p.denominator for row in PAR['quadratic'] for p in row))
def integer(q):
    q=Q(q)*DEN
    assert q.denominator==1
    return q.numerator
U0=integer(PAR['constant']);U1=list(map(integer,PAR['linear']))
U2=[[integer(p) for p in row] for row in PAR['quadratic']]

def plus(out,key,a,b=0):
    c,d=out.get(key,(0,0));out[key]=(c+a,d+b)

def car_action(state,index,create):
    occupied=(state>>index)&1
    if occupied==create:return None
    sign=-1 if (state&((1<<index)-1)).bit_count()%2 else 1
    return state^(1<<index),sign

def exact_mass_terms(axis):
    # Coefficients in units of 1/100; block ordering is the original598 one.
    if axis==4:
        return [(((31,True),(30,True)),31,9),
                (((30,False),(31,False)),31,-9)]
    X=[(int(axis==j),int(axis==j+2)) for j in range(2)]
    tilde=[(X[1][0],-X[1][1]),(-X[0][0],X[0][1])]
    ans=[]
    for left,right,col,colors,y in ((0,12,tilde,3,(120,30)),(0,18,X,3,(50,-10)),
                                  (24,30,tilde,1,(40,-12)),(24,28,X,1,(20,7))):
        for color in range(colors):
            for weak in range(2):
                ar,ai=col[weak];yr,yi=y;re=ar*yr-ai*yi;im=ar*yi+ai*yr
                if not(re or im):continue
                for spin in range(2):
                    i=left+4*color+2*weak+spin;j=right+2*color+spin
                    ans.extend([(((j,False),(i,True)),re,im),
                                (((i,False),(j,True)),re,-im)])
    return ans

MASS_TERMS=[exact_mass_terms(i) for i in range(5)]
@lru_cache(None)
def mass_steps(axis,state):
    accumulated={}
    for actions,re,im in MASS_TERMS[axis]:
        target=state;sign=1
        for index,create in actions:
            step=car_action(target,index,create)
            if step is None:break
            target,sg=step;sign*=sg
        else:
            plus(accumulated,target,sign*re,sign*im)
    rows=tuple((target,a*(DEN//100),b*(DEN//100))
               for target,(a,b) in accumulated.items() if a or b)
    # Independent calibration against the old floating CAR implementation.
    approximate=old.car_step(axis,state)
    reconstructed={target:complex(float(Q(a,DEN)),float(Q(b,DEN))) for target,a,b in rows}
    assert set(reconstructed)==set(approximate)
    assert all(abs(z-approximate[target])<1e-14 for target,z in reconstructed.items())
    return rows

def H(poly):
    out={}
    for (state,p),(re,im) in poly.items():
        degree=sum(p)
        c=integer(Q(5,2)+Q(2*degree,3)-Q(degree*degree,12))+U0
        plus(out,(state,p),c*re,c*im)
        for i,d in enumerate(p):
            if d>=2:
                q=list(p);q[i]-=2;c=-DEN*d*(d-1)//2
                plus(out,(state,tuple(q)),c*re,c*im)
            q=list(p);q[i]+=2
            c=DEN*degree//6+U1[int(i==4)]
            plus(out,(state,tuple(q)),c*re,c*im)
            for j in range(5):
                qq=q.copy();qq[j]+=2
                c=-DEN//12+U2[int(i==4)][int(j==4)]
                plus(out,(state,tuple(qq)),c*re,c*im)
            q=list(p);q[i]+=1
            for target,a,b in mass_steps(i,state):
                plus(out,(target,tuple(q)),a*re-b*im,a*im+b*re)
    return {k:z for k,z in out.items() if z!=(0,0)}

def jets(state):
    ans=[{(state,ZERO):(1,0)}]
    for _ in range(6):ans.append(H(ans[-1]))
    return ans

def radial(p):
    ans={}
    for (state,x),z in p.items():
        if any(x[1:4]):continue
        ans.setdefault(state,{})[(x[0],x[4])]=z
    return ans

def density(jet,n):
    terms=list(map(radial,jet));out={}
    for k in range(n+1):
        a,b=terms[n-k],terms[k]
        # i^n is treated separately for odd n.
        factor=(-1)**k*comb(n,k)
        for state,pa in a.items():
            if state not in b:continue
            for (r,s),(ar,ai) in pa.items():
                for (rr,ss),(br,bi) in b[state].items():
                    re=ar*br+ai*bi;im=ar*bi-ai*br
                    for _ in range(n%4):re,im=-im,re
                    plus(out,(r+rr,s+ss),factor*re,factor*im)
    return {k:z for k,z in out.items() if z!=(0,0)}

def run():
    a,b=jets(0),jets((1<<30)|(1<<31));rows={}
    for n in range(7):
        da,db=density(a,n),density(b,n)
        delta={k:(db.get(k,(0,0))[0]-da.get(k,(0,0))[0],
                  db.get(k,(0,0))[1]-da.get(k,(0,0))[1]) for k in set(da)|set(db)}
        delta={k:z for k,z in delta.items() if z!=(0,0)}
        assert all(z[1]==0 for z in delta.values())
        coeff={k:Q(z[0],DEN**n) for k,z in delta.items()}
        rows[str(n)]=[dict(r=k[0],s=k[1],numerator=str(z.numerator),
                          denominator=str(z.denominator),float_value=float(z)) for k,z in sorted(coeff.items())]
    assert all(not rows[str(n)] for n in (0,1,2,3,5))
    assert len(rows['4'])==8 and len(rows['6'])==26
    # Agreement is calibration, not part of the integer certificate.
    saved=json.loads((ROOT/'round744_drafts/history_density_probe_results.json').read_text('utf8'))
    error=[]
    for n in (4,6):
        older={(v['r_power'],v['s_power']):v['real'] for v in saved[str(n)]['coefficients']}
        for v in rows[str(n)]:error.append(abs(v['float_value']-older[(v['r'],v['s'])]))
    assert max(error)<2e-9
    return dict(input_interpretation='Exact frozen decimal matching-table and Yukawa data; exact rational u=L^-1 C.',
        common_operator_denominator_digits=len(str(DEN)),terms_empty=list(map(len,a)),terms_pair=list(map(len,b)),
        exact_zero_density_orders=[0,1,2,3,5],nonzero_density_orders=[4,6],coefficients=rows,
        maximum_difference_from_floating_density=max(error),interval_integral_certified=False,
        no_finite_time_or_graph_extension_claim=True,all_checks_passed=True)

if __name__=='__main__':
    result=run()
    target=Path(__file__).with_name('exact_history_density_results.json')
    with target.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:result[k] for k in ('exact_zero_density_orders','maximum_difference_from_floating_density','common_operator_denominator_digits','all_checks_passed')}))
