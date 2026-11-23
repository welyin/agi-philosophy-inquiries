"""Independent proof checks from two-color loop contraction, fusion and Young paths."""
import argparse
from itertools import product
import json
import math
from pathlib import Path

import numpy as np


def matchings(points):
    if not points:
        yield ()
        return
    a,*tail=points
    for i,b in enumerate(tail):
        for rest in matchings(tail[:i]+tail[i+1:]):
            yield ((a,b),*rest)


def loops(left,right,k):
    # Each color index may independently take three values on a connected cycle.
    next_vertex=[[] for _ in range(2*k)]
    for a,b in (*left,*right):
        next_vertex[a].append(b);next_vertex[b].append(a)
    todo=set(range(2*k))
    count=0
    while todo:
        count+=1
        stack=[todo.pop()]
        while stack:
            for v in next_vertex[stack.pop()]:
                if v in todo:
                    todo.remove(v);stack.append(v)
    return count


def folded(route,k,weights):
    tops=[edge for edge in route if edge[1]<k]
    bottoms=[edge for edge in route if edge[0]>=k]
    cross=[edge for edge in route if edge[0]<k<=edge[1]]
    assert len(tops)==len(bottoms)
    for choices in product(range(3),repeat=len(tops)):
        edges=cross.copy()
        weight=1.
        for (a,b),(c,d),choice in zip(tops,bottoms,choices):
            edges+=(((a,c),(b,d)),((a,d),(b,c)),((a,b),(c,d)))[choice]
            weight*=weights[choice]
        yield tuple(sorted(edges)),weight


def fusion(k):
    spins={0:1}
    for _ in range(2*k):
        nxt={}
        for spin,amount in spins.items():
            for out in range(abs(spin-1),spin+2):
                nxt[out]=nxt.get(out,0)+amount
        spins=nxt
    shapes={(0,0,0):1}
    for _ in range(k):
        nxt={}
        for row,amount in shapes.items():
            for i in range(3):
                if i>0 and row[i]>=row[i-1]:continue
                out=list(row);out[i]+=1
                out=tuple(out)
                nxt[out]=nxt.get(out,0)+amount
        shapes=nxt
    return {'orthogonal':spins[0],'unitary':sum(c*c for c in shapes.values())}


def calculate(k):
    pairings=list(matchings(list(range(2*k))))
    index={m:i for i,m in enumerate(pairings)}
    m=len(index)
    pair_gram=np.array([[3.**(loops(a,b,k)-k) for b in pairings] for a in pairings])
    theta=np.pi/np.sqrt(3)
    antisym=np.cos(theta)
    report={'routes':m}
    for family,plus,third in [('unitary',np.cos(theta*np.sqrt(2.5)),0.),
                              ('orthogonal',-1.,2/3)]:
        coeffs=((plus+antisym)/2,(plus-antisym)/2,third)
        t=np.zeros((m,m))
        for j,route in enumerate(pairings):
            for out,w in folded(route,k,coeffs):
                t[index[out],j]+=w
        gram=t.T @ pair_gram @ t
        spec=np.linalg.eigvalsh(gram)
        positive=spec[spec>1e-10]
        rank=fusion(k)[family]
        assert len(positive)==rank
        assert positive[0]>.01 and np.max(np.abs(spec[:m-rank]),initial=0)<1e-11
        flag=.01*rank/m**2
        report[family]={'rank':rank,'nullity':m-rank,'least_positive':float(positive[0]),
                        'largest':float(positive[-1]),'mixed_flag':flag}
    return report


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    else:
        np.testing.assert_allclose(a,b,rtol=2e-10,atol=1e-13)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--write-new',action='store_true')
    g.add_argument('--check',action='store_true')
    args=p.parse_args()
    results={str(k):calculate(k) for k in (3,4)}
    path=Path(__file__).with_name('independent_results.json')
    if args.write_new:
        with path.open('x',encoding='utf8') as f:
            json.dump(results,f,indent=2,allow_nan=False)
            f.write('\n')
    else:compare(results,json.loads(path.read_bytes()))
    print(json.dumps(results))
