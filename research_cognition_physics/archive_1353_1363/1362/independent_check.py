"""Finite tensor checks for the two unconditioned readouts and a tetrahedral competitor."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def permutation_matrix(p):
    n=len(p)
    m=np.zeros((3**n,3**n))
    for word in itertools.product(range(3),repeat=n):
        pos=np.ravel_multi_index(word,(3,)*n)
        out=np.ravel_multi_index(tuple(word[t] for t in p),(3,)*n)
        m[out,pos]=1
    return m


def antisym(n):
    op=np.zeros((3**n,3**n))
    for p in itertools.permutations(range(n)):
        inversions=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        op+=(-1)**inversions*permutation_matrix(p)/math.factorial(n)
    assert np.allclose(op@op,op)
    return op


def tensor_slots(rho,n):
    out=np.eye(1)
    for _ in range(n):out=np.kron(out,rho)
    return out


def word_test(weights,basis):
    m=basis@np.diag(np.sqrt(weights))@basis.T
    w=m.reshape(-1)
    i=np.eye(9)
    s=permutation_matrix((1,0))
    p=np.outer(w,w.conj())
    t=math.pi/math.sqrt(3)
    x=(-1+math.cos(t))/2
    y=(-1-math.cos(t))/2
    k=x*i+y*s+2*p
    lcu=(k-x*i-y*s)/(1+abs(x)+abs(y))
    np.testing.assert_allclose(lcu,p,atol=5e-15)
    packet=lcu.reshape(-1)/3
    yield_one=float(np.vdot(packet,packet).real)
    assert abs(yield_one-1/9)<1e-14
    ket=packet.reshape(3,-1)/np.sqrt(yield_one)
    rho=ket@ket.conj().T
    probs=[]
    for n in (2,3):
        prob=float(np.trace(antisym(n)@tensor_slots(rho,n)).real)
        probs.append(max(prob,0))
    np.testing.assert_allclose(probs,[(1-sum(a*a for a in weights))/2,np.prod(weights)],atol=4e-15)
    vectors=[a.reshape(-1)/3 for a in (i,s,k)]
    gram=np.array([[np.vdot(a,b) for b in vectors] for a in vectors])
    return gram,{'post_corridor':probs,'all_attempts':[p2*yield_one**n for n,p2 in zip((2,3),probs)]}


def subgroups():
    cycle=np.eye(3)[[1,2,0]]
    group=[]
    for signs in itertools.product((-1,1),repeat=3):
        if np.prod(signs)==1:
            for j in range(3):group.append(np.diag(signs)@np.linalg.matrix_power(cycle,j))
    assert len({tuple(g.flatten()) for g in group})==12
    test=np.zeros(27)
    for p in itertools.permutations((0,1,2)):test[np.ravel_multi_index(p,(3,3,3))]=1/np.sqrt(6)
    assert np.linalg.norm(antisym(3)@test)<1e-14
    for g in group:
        np.testing.assert_allclose(tensor_slots(g,3)@test,test)
    rot=np.array([[1,-1,0],[1,1,0],[0,0,np.sqrt(2)]])/np.sqrt(2)
    assert np.linalg.norm(tensor_slots(rot,3)@test-test)>0.5
    tetra=[sum(round(np.trace(g))**n for g in group)//12 for n in (3,4)]
    assert tetra==[2,7]
    # Nontrivial proper rotation made real by complex Takagi coordinates.
    z=np.exp(2j*np.pi/3)
    assert np.allclose(np.linalg.det(z*rot),1) and np.allclose((z*rot)@(z*rot).T,z*z*np.eye(3))
    phase=np.exp(1j*np.pi/5)
    assert abs(np.linalg.det(phase*rot)-1)>1 and np.allclose((phase*rot)@(phase*rot).T,phase**2*np.eye(3))
    return {"tetra_ranks_at_three_and_balanced_four":tetra,
            "extra_symmetric_tensor_rotation_overlap":float(np.vdot(test,tensor_slots(rot,3)@test).real)}


def answer():
    rng=np.random.default_rng(1362)
    b=np.linalg.qr(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))[0]
    out={}
    ref=None
    for name,weights in [('equal',(1/3,1/3,1/3)),('missing',(1/2,1/2,0)),('full',(2/3,1/6,1/6))]:
        gram,detail=word_test(weights,b)
        if ref is None:ref=gram
        np.testing.assert_allclose(gram,ref,atol=1e-14)
        out[name]=detail
    return {'round':1362,'direct_readouts':out,'stabilizer_tests':subgroups(),
            'max_common_abs_error_for_two_nonuniform_triple_windows':1/(108*9**3)}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    else:np.testing.assert_allclose(a,b,rtol=1e-10,atol=3e-14)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--check',action='store_true');group.add_argument('--write-new',action='store_true')
    args=parser.parse_args()
    data=answer()
    path=Path(__file__).with_name('independent_results.json')
    if args.write_new:
        with path.open('x',encoding='utf8') as f:
            json.dump(data,f,indent=2,allow_nan=False);f.write('\n')
    else:compare(data,json.loads(path.read_bytes()))
    print(json.dumps(data))

