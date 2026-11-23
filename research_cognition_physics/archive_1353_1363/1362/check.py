"""Finite tensor sources, rank-one preparation checks and subgroup alternatives."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def sign(p):
    return (-1)**sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1, len(p)))


def permutation(k, p):
    b = list(itertools.product(range(3), repeat=k))
    index = {v:i for i,v in enumerate(b)}
    a = np.zeros((3**k, 3**k))
    for i, v in enumerate(b):
        a[index[tuple(v[j] for j in p)],i] = 1.
    return a


def power(u, n):
    a = np.ones((1,1), dtype=complex)
    for _ in range(n):
        a = np.kron(a, u)
    return a


def tensor_checks():
    cycle = np.roll(np.eye(3),1,axis=0)
    tetra = [np.diag([a,b,a*b]) @ np.linalg.matrix_power(cycle,k)
             for a,b,k in itertools.product([1.,-1.],[1.,-1.],range(3))]
    zeta = np.exp(2j*math.pi/3)
    eps, sym = np.zeros((3,3,3)), np.zeros((3,3,3))
    for p in itertools.permutations(range(3)):
        eps[p], sym[p] = sign(p)/math.sqrt(6), 1/math.sqrt(6)
    eps, sym = eps.ravel(), sym.ravel()
    phi = .21
    rot = np.array([[math.cos(phi),-math.sin(phi),0],[math.sin(phi),math.cos(phi),0],[0,0,1]])
    d = np.eye(3).ravel()/math.sqrt(3)
    s = np.outer(d,d)
    for u in [rot]+tetra:
        u2 = power(u,2)
        assert np.linalg.norm(u2 @ s @ u2.T-s) < 1.e-12
        assert np.linalg.norm(power(u,3)@eps-eps) < 1.e-12
    assert np.linalg.norm(power(rot,3)@sym-sym) > .1
    records = {}
    for name, group in [("A4", tetra),("mu3_A4",[zeta**j*u for j in range(3) for u in tetra])]:
        inv = sum(power(u,3) for u in group)/len(group)
        mix = sum(np.kron(power(u,2),power(u.conj(),2)) for u in group)/len(group)
        assert np.linalg.norm(inv@eps-eps) < 1.e-12
        assert np.linalg.norm(inv@sym-sym) < 1.e-12
        assert np.linalg.norm(inv@inv-inv) < 1.e-12
        assert np.linalg.norm(mix@mix-mix) < 1.e-12
        records[name] = [int(np.linalg.matrix_rank(inv,tol=1.e-8)),
                         int(np.linalg.matrix_rank(mix,tol=1.e-8))]
        assert records[name] == [2,7]
    even = sum(permutation(3,p) for p in itertools.permutations(range(3)))/6
    emitted = (eps+sym)/math.sqrt(2)
    scalar_rates = [float(np.linalg.norm(even@emitted)**2 * math.sin(t)**2)
                    for t in [math.pi/2,math.pi/4]]
    assert np.allclose(scalar_rates,[.5,.25])
    return records, scalar_rates


def metric_data(weights):
    identity = np.eye(9)
    swap = permutation(2,[1,0])
    w = np.diag(np.sqrt(weights)).ravel()
    pw = np.outer(w,w)
    f1 = math.cos(math.pi/math.sqrt(3))
    x,y,z = (f1-1)/2,(-1-f1)/2,2.
    k = x*identity + y*swap + z*pw
    assert np.max(np.abs(np.linalg.eigvalsh(k))) < 1+1.e-12
    old = np.column_stack([q.ravel()/3 for q in [identity,swap,k]])
    lcu = (k-x*identity-y*swap)/(1+abs(x)+abs(y))
    assert np.allclose(lcu,pw)
    prepared = lcu/3
    prep = float(np.linalg.norm(prepared)**2)
    normalized = prepared/np.sqrt(prep)
    active = normalized @ normalized.T
    rho = np.einsum("iaja->ij",active.reshape(3,3,3,3))
    assert np.allclose(rho,np.diag(weights))
    anti2 = (np.eye(9)-swap)/2
    anti3 = sum(sign(p)*permutation(3,p) for p in itertools.permutations(range(3)))/6
    p2 = float(np.trace(anti2@np.kron(rho,rho)))
    p3 = float(np.trace(anti3@np.kron(np.kron(rho,rho),rho)))
    assert np.isclose(p2,(1-np.trace(rho@rho))/2)
    assert np.isclose(p3,np.linalg.det(rho))
    assert np.isclose(np.trace(active@active),1)
    assert np.isclose(abs(np.trace(normalized)/3)**2,1/9)
    return old.T@old, {"schmidt_probabilities":weights,"single_return":prep,
                       "two_copy_minus":p2,"three_copy_antisymmetric":p3,
                       "from_all_preparations":[prep**2*p2,prep**3*p3]}


def calculate():
    weights = {"metric":[1/3]*3, "rank1":[1.,0.,0.],
               "rank2":[.5,.5,0.],"other_full_rank":[2/3,1/6,1/6]}
    grams, tests = {}, {}
    for name,v in weights.items():
        grams[name],tests[name] = metric_data(v)
    assert all(np.allclose(g,grams["metric"]) for g in grams.values())
    assert np.isclose(tests["rank2"]["two_copy_minus"],tests["other_full_rank"]["two_copy_minus"])
    assert np.isclose(tests["other_full_rank"]["three_copy_antisymmetric"],1/54)
    ranks, rates = tensor_checks()
    # Two calibrated individual metrics may fail to be a shared projective metric.
    d = np.eye(3)/math.sqrt(3)
    other = np.diag(np.exp(2j*math.pi*np.arange(3)/3))/math.sqrt(3)
    both_metric_packet_minus=(1-abs(np.vdot(d,other))**4)/2
    assert np.isclose(both_metric_packet_minus,.5)
    return {"gram_all_four_preparations":grams["metric"].tolist(),"holdouts":tests,
            "two_individual_metrics_distinct_full_packet_minus":float(both_metric_packet_minus),
            "smaller_groups_dimensions_triple_and_balanced_four":ranks,
            "example_extra_symmetric_channel_joint_probabilities":rates}


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:
            compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):
            compare(x,y)
    else:
        assert np.isclose(a,b,rtol=1.e-10,atol=1.e-12)


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--write",action="store_true")
    args=cli.parse_args()
    result=calculate()
    path=Path(__file__).with_name("results.json")
    if args.write:
        with path.open("x",encoding="utf8") as f:
            json.dump(result,f,indent=2,allow_nan=False)
            f.write("\n")
    else:
        compare(result,json.loads(path.read_text(encoding="utf8")))
    print(json.dumps(result,indent=2))
