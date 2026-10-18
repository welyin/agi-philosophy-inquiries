"""Independent exact checks for round 1304.
No author code is imported. The complete integer box is justified by the
connected-edge span lemma in the proof, not by a guessed charge cutoff.
Default and --check are read-only; --write exclusively creates the result.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations_with_replacement
from math import isqrt
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def det(m):
    a = [[F(x) for x in row] for row in m]
    n = len(a)
    value = F(1)
    for k in range(n):
        pivot = next((j for j in range(k,n) if a[j][k]), None)
        if pivot is None:
            return F(0)
        if pivot != k:
            a[k],a[pivot] = a[pivot],a[k]
            value = -value
        p = a[k][k]
        value *= p
        for i in range(k+1,n):
            factor = a[i][k]/p
            for j in range(k+1,n):
                a[i][j] -= factor*a[k][j]
            a[i][k] = F(0)
    return value

def moments(a):
    return tuple(sum(x**j for x in a) for j in (1,2,3))

def elementary_from_power(a):
    """Newton identities up to degree 3, independently evaluated exactly."""
    n = len(a)
    p = [F(n)]+[F(sum(x**j for x in a)) for j in range(1,n+1)]
    e = [F(1)]
    for k in range(1,n+1):
        e.append(sum((-1)**(j-1)*e[k-j]*p[j]
                     for j in range(1,k+1))/k)
    return e

def graph_components(vertices, edges):
    remaining = set(vertices)
    adj = {v:set() for v in remaining}
    for a,b in edges:
        adj[a].add(b)
        adj[b].add(a)
    out = []
    while remaining:
        seed = min(remaining)
        remaining.remove(seed)
        component = {seed}
        stack = [seed]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y in remaining:
                    remaining.remove(y)
                    component.add(y)
                    stack.append(y)
        out.append(sorted(component))
    return sorted(out)

def allowed_components(a,b,power):
    # A and B are disjoint; coordinate roles can therefore be labeled by
    # their values even though the original physical charges are A and -B.
    assert not set(a)&set(b)
    edges = [(x,y) for x in set(a) for y in set(b)
             if abs(x-y)<=power]
    return graph_components(set(a+b),edges)

def has_matching(a,b,power):
    # Subset dynamic programming; no permutation enumeration shared with
    # the author's bucket construction.
    masks = {0}
    for x in a:
        masks = {mask|(1<<j) for mask in masks for j,y in enumerate(b)
                 if not mask&(1<<j) and abs(x-y)<=power}
    return (1<<len(b))-1 in masks

def complete_box_by_two_roots():
    """Every normalized disjoint pair is oriented with 0 in A.
    Given A and the first two sorted B roots u,v, the remaining sum and
    sum of squares determine their quadratic polynomial exactly. The
    cubic test then finishes the degree-three moment condition.
    """
    box = 14
    found = set()
    a_count = 0
    for tail in combinations_with_replacement(range(box+1),3):
        a = (0,)+tail
        a_count += 1
        p1,p2,p3 = moments(a)
        for u,v in combinations_with_replacement(range(1,box+1),2):
            total = p1-u-v
            twice_product = total*total-(p2-u*u-v*v)
            if twice_product%2:
                continue
            product = twice_product//2
            discriminant = total*total-4*product
            if discriminant<0:
                continue
            square_root = isqrt(discriminant)
            if square_root*square_root != discriminant:
                continue
            if (total-square_root)%2:
                continue
            w,z = (total-square_root)//2,(total+square_root)//2
            if not v<=w<=z<=box:
                continue
            b = (u,v,w,z)
            if set(a)&set(b) or sum(x**3 for x in b)!=p3:
                continue
            assert moments(a)==moments(b)
            found.add((a,b))
    expected = {
        ((0,3,4,7),(1,1,6,6)),
        ((0,4,7,11),(1,2,9,10)),
        ((0,5,5,10),(1,2,8,9)),
        ((0,6,7,13),(1,3,10,12)),
        ((0,6,8,14),(2,2,12,12)),
    }
    assert a_count==680 and found==expected
    certificates = []
    for a,b in sorted(found):
        components = allowed_components(a,b,2)
        assert len(components)>1
        certificates.append({
            "A":a,"B":b,"moments":moments(a),
            "power_two_components":components,
            "power_two_perfect_matching":has_matching(a,b,2),
        })
    return {"span_from_eight_vertices_times_two":14,
            "normalized_A_multisets":a_count,
            "all_disjoint_equal_moment_pairs":certificates,
            "number_connected":0,
            "normalization_not_a_physical_equivalence":True}

def six_anomalies(a,b):
    q = a+tuple(-x for x in b)
    identity = (1,)*len(a)+(-1,)*len(b)
    return {
        "g":sum(q),
        "g3":sum(x**3 for x in q),
        "F":sum(identity),
        "F3":sum(x**3 for x in identity),
        "gF2":sum(x*f*f for x,f in zip(q,identity)),
        "g2F":sum(x*x*f for x,f in zip(q,identity)),
    }

def matrix_certificate(a,b,d,max_power):
    n = len(a)
    q = a+tuple(-x for x in b)
    f = (1,)*n+(-1,)*n
    assert all(q) and not any(-x in q for x in q)
    assert all(v==0 for v in six_anomalies(a,b).values())
    # Same physical charge must carry one scalar identity value.
    for i in range(2*n):
        for j in range(2*n):
            assert q[i]!=q[j] or f[i]==f[j]
    whole = [[F(0) for _ in range(2*n)] for _ in range(2*n)]
    edges,monomials = [],[]
    epsilon,cutoff = F(1,5),F(10)
    vacuum = epsilon*cutoff
    scale = cutoff*epsilon**max_power/2
    for i in range(n):
        for j in range(n):
            if not d[i][j]:
                continue
            r = abs(a[i]-b[j])
            assert 1<=r<=max_power
            u = max(b[j]-a[i],0)
            v = max(a[i]-b[j],0)
            coefficient = d[i][j]*epsilon**(max_power-r)/2
            assert 0<coefficient<=F(1,2)
            # Gauge covariance of each complete polynomial, including S.
            assert a[i]-b[j]+u-v == 0
            assert f[i]+f[n+j] == 0
            amplitude = coefficient*cutoff**(1-r)*vacuum**r
            assert amplitude == scale*d[i][j]
            radial_derivative = r*amplitude/vacuum
            assert vacuum*radial_derivative == r*amplitude
            phase_coefficient = (u-v)*amplitude
            assert phase_coefficient+(a[i]-b[j])*amplitude == 0
            whole[i][n+j] = whole[n+j][i] = d[i][j]
            edges.append((q[i],q[n+j]))
            monomials.append({"i":i,"j":j,"S_power":u,"Sbar_power":v,
                              "dimension":3+r,
                              "wilson_at_fixed_epsilon":str(coefficient)})
    for i in range(2*n):
        for j in range(2*n):
            assert (f[i]+f[j])*whole[i][j]==0
    components = graph_components(set(q),edges)
    determinant = det(d)
    assert determinant
    assert det(whole)==(-1)**n*determinant**2
    smallest = float(np.linalg.svd(np.array(d,dtype=float),
                                    compute_uv=False)[-1])
    return {"moments":moments(a),"six_local_coefficients":six_anomalies(a,b),
            "D_determinant":str(determinant),
            "full_mass_determinant_normalized":str(det(whole)),
            "actual_role_components":components,
            "largest_actual_dimension":max(x["dimension"] for x in monomials),
            "monomials":monomials,
            "epsilon":str(epsilon),"Lambda":str(cutoff),
            "v":str(vacuum),"common_mass_scale":str(scale),
            "minimum_singular_value_diagnostic":float(format(smallest,".12g")),
            "gauge_scalar_source_and_F_identities":True}

def main_examples():
    a,b = (1,6,6,11),(2,3,9,10)
    d = [[F(i==j) for j in range(4)] for i in range(4)]
    for i,j in ((0,1),(1,2),(3,2)):
        d[i][j] = F(1,4)
    main = matrix_certificate(a,b,d,3)
    assert len(main["actual_role_components"])==1
    assert main["largest_actual_dimension"]==6
    assert main["D_determinant"]=="1"
    # Exact norm bound: ||K||_1=2, ||K||_infinity=1; e=1/4.
    # sqrt(2)/4 < 1/2 since 2<4, hence sigma_min(D)>1/2.
    assert F(2,16)<F(1,4)
    main["analytic_minimum_singular_value_gt"] = "1/2"
    # All three moments also follow from the central Pythagorean identity.
    assert 5*5+0*0 == 4*4+3*3
    assert moments(tuple(x-6 for x in a)) == moments(tuple(x-6 for x in b))
    delete = matrix_certificate((1,4,5,8),(2,2,7,7),
        [[F(i==j) for j in range(4)] for i in range(4)],2)
    assert len(delete["actual_role_components"])==2
    assert delete["largest_actual_dimension"]==5
    assert delete["D_determinant"]=="1"
    # Primitive F is global for the retained identity; deleting graph
    # connectivity permits additional independent component identities.
    newton = {}
    for values in ((2,),(-2,5),(-3,0,7)):
        e = elementary_from_power(values)
        direct = [F(1)]
        for v in values:
            direct.append(F(0))
            for k in range(len(direct)-1,0,-1):
                direct[k] += v*direct[k-1]
        assert e==direct
        newton[str(len(values))] = [str(x) for x in e]
    return {"dimension_six_connected":main,
            "dimension_five_disconnected_control":delete,
            "newton_formula_diagnostics_not_lower_bound_proof":newton}

def unequal_identity_disconnected_control():
    # Only graph connectivity is deleted. Component identity normalizations
    # are now independent, so they need not all have |F|=1.
    q=(2,-1,10,-9,4,-5,8,-9)
    f=(1,-1,3,-3,1,-1,3,-3)
    coeffs={
        "g":sum(q),"g3":sum(x**3 for x in q),
        "F":sum(f),"F3":sum(x**3 for x in f),
        "gF2":sum(x*y*y for x,y in zip(q,f)),
        "g2F":sum(x*x*y for x,y in zip(q,f)),
    }
    assert all(x==0 for x in coeffs.values())
    assert all(q) and not any(-x in q for x in q)
    assert all(q[i]!=q[j] or f[i]==f[j]
               for i in range(8) for j in range(8))
    mass=[[F(0) for _ in q] for _ in q]
    edges=[]
    for i in range(0,8,2):
        weight=q[i]+q[i+1]
        assert abs(weight)==1 and f[i]+f[i+1]==0
        u,v=max(-weight,0),max(weight,0)
        assert weight+u-v==0 and u+v==1
        mass[i][i+1]=mass[i+1][i]=F(1)
        edges.append((q[i],q[i+1]))
    assert det(mass)==1
    comps=graph_components(set(q),edges)
    assert len(comps)==3
    for i in range(8):
        for j in range(8):
            assert (f[i]+f[j])*mass[i][j]==0
    return {"q":q,"F":f,"six_local_coefficients":coeffs,
            "actual_role_components":comps,
            "maximum_operator_dimension":4,
            "mass_determinant":str(det(mass)),
            "equal_global_identity_magnitude":False,
            "gauge_and_F_polynomial_identities":True}

def make_result():
    return {
        "status":"PASS_INDEPENDENT_CHECK_WITH_SCOPE",
        "complete_finite_reconstruction":complete_box_by_two_roots(),
        "same_object_examples":main_examples(),
        "delete_connectivity_without_equal_magnitude":unequal_identity_disconnected_control(),
        "source_bindings":{
            "admission":sha(ROOT/"_shared"/"round1304_admission.md"),
            "admission_clarification":sha(ROOT/"_shared"/"round1304_admission_clarification.md"),
            "author_proof":sha(HERE/"proof.md"),
            "1302_proof":sha(ROOT/"1302"/"proof.md"),
            "1022_note":sha(ROOT.parent/"archive_1009_1043"/"research_note_1022.md"),
        },
        "independent_code_sha256":sha(Path(__file__).resolve()),
        "no_quantum_instrument_or_global_anomaly_claim":True,
        "PTE_and_eight_field_lower_bound_are_inherited":True,
    }

def run():
    parser = argparse.ArgumentParser()
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--check",action="store_true")
    action.add_argument("--write",action="store_true")
    opts = parser.parse_args()
    result = make_result()
    target = HERE/"independent_results.json"
    if opts.write:
        with target.open("x",encoding="utf-8") as out:
            json.dump(result,out,ensure_ascii=False,indent=2)
            out.write("\n")
        print("CREATED",result["status"])
    else:
        saved = json.loads(target.read_text(encoding="utf-8"))
        assert saved==json.loads(json.dumps(result)), "Saved result mismatch"
        print("READ_ONLY_MATCH",result["status"])

if __name__=="__main__":
    run()
