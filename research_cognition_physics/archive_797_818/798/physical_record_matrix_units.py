"""Round 798: exact formal lifting of a parity-even two-mode record algebra.

Finite rational diagnostics of the algebraic formulas, not a continuum
Hamiltonian simulation, finite-coupling convergence, or blank preparation.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import json
import sys
sys.dont_write_bytecode = True
import numpy as np
import internal_register_probe as probe
C, matrix, scale, dag, same = probe.C, probe.matrix, probe.scale, probe.dag, probe.same
HERE = Path(__file__).resolve().parent
ORDER = 4
I = probe.eye(4)
Z = scale(I, 0)

def clean(s):
    return {k: a for k, a in s.items() if any(a.flat)}

def add(*args):
    out = {}
    for s in args:
        for k, a in s.items():
            out[k] = out.get(k, Z) + a
    return clean(out)

def smul(s, c):
    return clean({k: scale(a, c) for k, a in s.items()})

def mul(s, t):
    out = {}
    for k, a in s.items():
        for l, b in t.items():
            if k+l <= ORDER:
                out[k+l] = out.get(k+l, Z) + a @ b
    return clean(out)

def adj(s):
    return clean({k: dag(a) for k, a in s.items()})

def const(a):
    return clean({0: a})

ONE = const(I)

def sub(a, b):
    return add(a, smul(b, -1))

def equal(a, b):
    return not sub(a, b)

def invsqrt(defect, unit=ONE):
    assert 0 not in defect
    out, power, coefficient = unit, unit, Q(1)
    for n in range(1, ORDER+1):
        power = mul(power, defect)
        coefficient *= (Q(-1, 2) - (n-1))/n
        out = add(out, smul(power, coefficient))
    return out

def lift(a, b):
    assert equal(adj(a), a)
    bb = sub(smul(a, 2), ONE)
    defect = sub(mul(bb, bb), ONE)
    p = smul(add(ONE, mul(bb, invsqrt(defect))), Q(1, 2))
    q = sub(ONE, p)
    v = mul(mul(q, b), p)
    corner_defect = sub(mul(adj(v), v), p)
    u = mul(v, invsqrt(corner_defect, p))
    return {(0,0): p, (1,1): q, (1,0): u, (0,1): adj(u)}

def expectation(e, x):
    return smul(add(*(mul(mul(e[i,j], x), e[j,i])
                       for i in range(2) for j in range(2))), Q(1,2))

def first_defect(d):
    k = min(d)
    a = d[k]
    i,j = next((i,j) for i in range(4) for j in range(4) if a[i,j])
    return dict(order=k, row=i, column=j, value=str(a[i,j]))

def run():
    c1, c2 = probe.car(2)
    n1, n2 = dag(c1) @ c1, dag(c2) @ c2
    parity = (I-scale(n1,2)) @ (I-scale(n2,2))
    pair = dag(c1) @ dag(c2)
    a1 = {0:c1, 1:scale(c2,Q(1,3))+scale(dag(c1),C(0,Q(1,5))),
          2:scale(c1@n2,Q(1,7))+scale(dag(c2)@n1,Q(1,11))}
    a2 = {0:c2, 1:scale(c1,Q(1,4))+scale(dag(c2),Q(1,6)),
          2:scale(c2@n1,Q(1,9))}
    a = sub(ONE, mul(adj(a1),a1))
    b = mul(adj(a1), add(a2,adj(a2)))
    raw_defect = sub(mul(a,a),a)
    assert raw_defect
    e = lift(a,b)
    p,q,u = e[0,0],e[1,1],e[1,0]
    assert equal(add(p,q),ONE)
    assert equal(mul(adj(u),u),p) and equal(mul(u,adj(u)),q)
    relation_count = 0
    for i in range(2):
        for j in range(2):
            assert equal(adj(e[i,j]),e[j,i])
            assert equal(mul(const(parity),e[i,j]),mul(e[i,j],const(parity)))
            for k in range(2):
                for l in range(2):
                    target = e[i,l] if j==k else {}
                    assert equal(mul(e[i,j],e[k,l]),target),(i,j,k,l)
                    relation_count += 1
    # A raw interacting occupation and pairing do not satisfy matrix relations.
    raw_v = mul(mul(sub(ONE,a),b),a)
    raw_norm_defect = sub(mul(adj(raw_v),raw_v),a)
    assert raw_norm_defect
    assert equal(expectation(e,ONE),ONE)
    queries = [const(pair+dag(pair)),const(n1),
               add(const(parity),{1:n1}),const(c1+dag(c1))]
    for x in queries:
        ex = expectation(e,x)
        assert equal(expectation(e,ex),ex)
        for unit in e.values():
            assert equal(mul(ex,unit),mul(unit,ex))
    x = queries[2]
    changed = sub(expectation(e,x),x)
    assert changed and min(changed)==1
    # CP Gram formula at one nontrivial element: E(x^*x)=sum y_ij^*y_ij.
    x = add(queries[0],smul(queries[2],C(0,1)))
    gram = smul(add(*(mul(adj(mul(x,e[j,i])),mul(x,e[j,i]))
                     for i in range(2) for j in range(2))),Q(1,2))
    assert equal(expectation(e,mul(adj(x),x)),gram)
    # Naturality under one common *-automorphism, including source coefficients.
    g = pair+dag(pair)+scale(n1,Q(1,3))
    generator = {1:scale(g,C(0,1))}
    unitary, power = ONE, ONE
    for j in range(1,ORDER+1):
        power=mul(power,generator)
        unitary=add(unitary,smul(power,Q(1,factorial(j))))
    assert equal(mul(adj(unitary),unitary),ONE)
    def alpha(x):
        return mul(mul(unitary,x),adj(unitary))
    transported = lift(alpha(a),alpha(b))
    for key in e:
        assert equal(transported[key],alpha(e[key]))
    # Fixed leading state is retained; constructing matrix units does not reset it.
    rho = scale(I,Q(1,4))
    occupied = probe.expect(rho,q[0])
    assert occupied==Q(1,2)
    return dict(round=798,all_checks_passed=True,epsilon_order=ORDER,
                exact_matrix_dimension=4,
                matrix_unit_product_identities=relation_count,
                product_coefficient_identities=relation_count*(ORDER+1),
                raw_occupation_projection_defect=first_defect(raw_defect),
                raw_partial_isometry_defect=first_defect(raw_norm_defect),
                parity_preserved_on_both_CAR_sectors=True,
                unital_full_M2_not_only_even_sector_corner=True,
                conditional_expectation_unital_idempotent_commutant=True,
                conditional_expectation_Gram_identity=True,
                averaging_changes_leading_commuting_source=first_defect(changed),
                common_star_transport_preserves_entire_lift=True,
                fixed_state_leading_record_one_probability=str(occupied),
                same_state_reset_to_blank=False,
                finite_matrix_probe_proves_continuum_lift=False,
                original_Hamiltonian_realizes_instrument=False)

if __name__=='__main__':
    result=run()
    (HERE/'physical_record_matrix_units_results.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))

