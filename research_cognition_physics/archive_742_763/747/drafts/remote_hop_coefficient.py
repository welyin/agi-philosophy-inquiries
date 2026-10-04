"""Exact coefficient of two kinetic, two mass and two sterile-hop insertions.

This computes the hop-amplitude-squared part of the remote sixth jet after
the analytic word/ready-parity exclusions in the working proof. It is not the
full remote jet or a replacement of the original scalar edge Hamiltonian.
"""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'round745_drafts'))
import exact_history_density as native
DEN=300;ZERO=(0,)*5;TARGET_COUNTS=(2,2,2)
HOP_MATRIX=(((100,0),(0,0)),((0,0),(100,0)))
def add(out,key,a,b=0):
    c,d=out.get(key,(0,0));out[key]=(c+a,d+b)
def clean(out):return {k:z for k,z in out.items() if z!=(0,0)}
def kinetic(poly):
    out={}
    for (state,p),(re,im) in poly.items():
        degree=sum(p);c=int(DEN*(Q(5,2)+Q(2*degree,3)-Q(degree*degree,12)))
        add(out,(state,p),c*re,c*im)
        for i,d in enumerate(p):
            if d>=2:
                q=list(p);q[i]-=2;c=-DEN*d*(d-1)//2
                add(out,(state,tuple(q)),c*re,c*im)
            q=list(p);q[i]+=2;c=DEN*degree//6
            add(out,(state,tuple(q)),c*re,c*im)
            for j in range(5):
                qq=q.copy();qq[j]+=2;c=-DEN//12
                add(out,(state,tuple(qq)),c*re,c*im)
    return clean(out)
@lru_cache(None)
def act_terms(state,axis):
    terms=native.MASS_TERMS[axis] if axis<5 else []
    if axis==5:
        for i in range(2):
            for j in range(2):
                re,im=HOP_MATRIX[i][j]
                if not(re or im):continue
                terms.extend([(((32+j,False),(30+i,True)),re,im),
                              (((30+i,False),(32+j,True)),re,-im)])
    out={}
    for actions,re,im in terms:
        target=state;sign=1
        for index,create in actions:
            step=native.car_action(target,index,create)
            if step is None:break
            target,sg=step;sign*=sg
        else:add(out,target,3*sign*re,3*sign*im)
    return tuple((k,*v) for k,v in clean(out).items())
def mass(poly):
    out={}
    for (state,p),(re,im) in poly.items():
        for axis in range(5):
            q=list(p);q[axis]+=1
            for target,a,b in act_terms(state,axis):
                add(out,(target,tuple(q)),a*re-b*im,a*im+b*re)
    return clean(out)
def hop(poly):
    out={}
    for (state,p),(re,im) in poly.items():
        for target,a,b in act_terms(state,5):add(out,(target,p),a*re-b*im,a*im+b*re)
    return clean(out)
def jets(state):
    levels=[{(0,0,0):{(state,ZERO):(1,0)}}]
    for _ in range(6):
        new={}
        for count,poly in levels[-1].items():
            for i,op in enumerate((kinetic,mass,hop)):
                if count[i]==2:continue
                nc=list(count);nc[i]+=1;nc=tuple(nc);dest=new.setdefault(nc,{})
                for key,(a,b) in op(poly).items():add(dest,key,a,b)
        levels.append({c:clean(p) for c,p in new.items() if clean(p)})
    return levels
def radial(poly):return native.radial(poly)
def density(levels):
    out={}
    for k in range(7):
        factor=-(-1)**k*comb(6,k)
        for ca,pa in levels[6-k].items():
            cb=tuple(2-c for c in ca)
            if cb not in levels[k]:continue
            aa,bb=radial(pa),radial(levels[k][cb])
            for state,terms in aa.items():
                if state not in bb:continue
                for (r,s),(ar,ai) in terms.items():
                    for (rr,ss),(br,bi) in bb[state].items():
                        add(out,(r+rr,s+ss),factor*(ar*br+ai*bi),factor*(ar*bi-ai*br))
    return clean(out)
def run():
    a,b=jets(0),jets((1<<32)|(1<<33));pa,pb=density(a),density(b);out={}
    for key in set(pa)|set(pb):
        ar,ai=pa.get(key,(0,0));br,bi=pb.get(key,(0,0))
        assert bi-ai==0
        if br!=ar:out[key]=Q(br-ar,DEN**6)
    prior=json.loads((ROOT/'round745_drafts/exact_history_density_results.json').read_text('utf8'))
    p4={(r['r'],r['s']):Q(int(r['numerator']),int(r['denominator'])) for r in prior['coefficients']['4']}
    ratios={out[k]/v for k,v in p4.items() if k in out}
    return dict(counts=list(TARGET_COUNTS),exact_coefficients=[
        dict(r=k[0],s=k[1],numerator=str(v.numerator),denominator=str(v.denominator),value=float(v)) for k,v in sorted(out.items())],
        support_equals_local_fourth=set(out)==set(p4),ratios_to_local_fourth=sorted(map(str,ratios)),
        total_terms_vacuum=[sum(map(len,q.values())) for q in a],
        total_terms_pair=[sum(map(len,q.values())) for q in b],
        scope='Coefficient extraction only. The full remote sixth coefficient also contains scalar-edge contributions.')
if __name__=='__main__':
    r=run()
    with Path(__file__).with_name('remote_hop_coefficient_results.json').open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(r,ensure_ascii=False))
