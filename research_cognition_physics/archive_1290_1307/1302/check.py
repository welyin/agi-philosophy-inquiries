"""1302 exact identity/anomaly and all-count certificates.
All-N and all-charge conclusions are analytic, not finite sample inference.
Default/--check are read-only. --write exclusively creates results.json.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import numpy as np
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results.json'


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p or [0]


def poly_add(a, b):
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0) for i in range(max(len(a), len(b)))])


def poly_scale(a, c):
    return trim([c*v for v in a])


def poly_mul(a, b):
    p = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            p[i+j] += x*y
    return trim(p)


def divide_monic(p, divisor):
    p, divisor = trim(p), trim(divisor)
    assert divisor[-1] == 1
    q = [0] * max(1, len(p) - len(divisor) + 1)
    while p != [0] and len(p) >= len(divisor):
        j, c = len(p) - len(divisor), p[-1]
        q[j] += c
        for i, v in enumerate(divisor):
            p[j+i] -= c*v
        p = trim(p)
    return trim(q), trim(p)


def path(x, k):
    assert x[0] > 0 and x[-1] > 0 and min(x) >= 0
    assert all(x[i] + x[i+1] > 0 for i in range(len(x)-1))
    n = [([0]+x)[i] + (x+[0])[i] for i in range(len(x)+1)]
    q = [(-1)**(i+1)*(k+i) for i in range(len(n))]
    f = [(-1)**i for i in range(len(n))]
    return q, n, f


def moments(x):
    return [sum((-1)**i * i**r * v for i, v in enumerate(x)) for r in range(3)]


def anomalies(q, n, f):
    return {'g': sum(c*z for z,c in zip(q,n)),
            'g3': sum(c*z**3 for z,c in zip(q,n)),
            'F': sum(c*z for z,c in zip(f,n)),
            'F3': sum(c*z**3 for z,c in zip(f,n)),
            'gF2': sum(c*z*w*w for z,w,c in zip(q,f,n)),
            'g2F': sum(c*z*z*w for z,w,c in zip(q,f,n))}


def polynomial_certificate(x):
    # k remains a formal indeterminate: no bounded charge search.
    _, n, f = path(x, 1)
    qpoly = [poly_scale([i, 1], (-1)**(i+1)) for i in range(len(n))]
    out = {key:[0] for key in ['g','g3','F','F3','gF2','g2F']}
    for p,c,w in zip(qpoly,n,f):
        out['g'] = poly_add(out['g'], poly_scale(p,c))
        out['g3'] = poly_add(out['g3'], poly_scale(poly_mul(poly_mul(p,p),p),c))
        out['F'] = poly_add(out['F'], [c*w])
        out['F3'] = poly_add(out['F3'], [c*w**3])
        out['gF2'] = poly_add(out['gF2'], poly_scale(p,c*w*w))
        out['g2F'] = poly_add(out['g2F'], poly_scale(poly_mul(p,p),c*w))
    m0,m1,m2 = moments(x)
    expected = {'g':[m0], 'g3':[m0+3*m1+3*m2,3*m0+6*m1,3*m0],
                'F':[0], 'F3':[0], 'gF2':[m0], 'g2F':[-m0-2*m1,-2*m0]}
    assert out == {key:trim(p) for key,p in expected.items()}
    quotient, remainder = divide_monic(x, [1,3,3,1])
    assert poly_add(poly_mul(quotient,[1,3,3,1]),remainder) == trim(x)
    divisible = remainder == [0]
    assert divisible == (m0 == m1 == m2 == 0)
    if divisible:
        assert 2*sum(x) == 16*sum(quotient) and sum(quotient) > 0
    return {'moments': [m0,m1,m2], 'all_k_coefficients_low_to_high':out,
            'division_quotient':quotient,'division_remainder':remainder,
            'triple_root_iff_three_zero_moments':True,
            'total_Weyl_fields':2*sum(x), 'divisible_by_one_plus_z_cubed':divisible}


def bareiss(matrix):
    a = [[int(v) for v in row] for row in matrix]
    prev, sign, n = 1, 1, len(a)
    for k in range(n-1):
        if a[k][k] == 0:
            row = next((j for j in range(k+1,n) if a[j][k]), None)
            if row is None:
                return 0
            a[k],a[row] = a[row],a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value = a[i][j]*pivot-a[i][k]*a[k][j]
                assert value % prev == 0
                a[i][j] = value // prev
        for i in range(k+1,n):
            a[i][k] = 0
        prev = pivot
    return sign*a[-1][-1]


def mass_certificate(x,k,denominator):
    q,n,f = path(x,k)
    charges = [z for z,c in zip(q,n) for _ in range(c)]
    identities = [w for w,c in zip(f,n) for _ in range(c)]
    N = sum(n)
    groups,offset = [],0
    for count in n:
        groups.append(list(range(offset,offset+count)))
        offset += count
    unused = [v[:] for v in groups]
    P = np.zeros((N,N),dtype=int)
    pairs = []
    for i,count in enumerate(x):
        for _ in range(count):
            a,b = unused[i].pop(0),unused[i+1].pop(0)
            P[a,b] = P[b,a] = 1
            pairs.append([a,b])
    assert not any(unused) and np.array_equal(P@P,np.eye(N,dtype=int))
    A = np.array([[int(abs(a+b)==1) for b in charges] for a in charges],dtype=int)
    numerator = denominator*P + A
    ff = np.array(identities)
    identity_residual = (ff[:,None]+ff[None,:])*numerator
    assert np.count_nonzero(identity_residual) == 0
    # Every adjacent distinct-charge block is actually nonzero, even at x_i=0.
    assert all(np.count_nonzero(numerator[np.ix_(groups[i],groups[i+1])]) > 0
               for i in range(len(n)-1))
    assert np.all((A==0)|(numerator!=0))
    assert all(z != 0 and -z not in charges for z in charges)
    degree = int(max(A.sum(axis=1)))
    bound = 1-Fraction(degree,denominator)
    assert bound > 0
    det = bareiss(numerator.tolist())
    assert det != 0
    M = numerator/denominator
    sv = np.linalg.svd(M,compute_uv=False)
    assert min(sv) >= float(bound)-1e-12
    angle = .237
    U = np.diag(np.exp(1j*angle*ff))
    assert np.linalg.norm(U.T@M@U-M) < 5e-12
    qq = np.array(charges)
    sums = qq[:,None]+qq[None,:]
    S = np.exp(1j*angle)
    new_M = S*M*(sums==-1)+S.conjugate()*M*(sums==1)
    G = np.diag(np.exp(1j*angle*qq))
    assert np.linalg.norm(G.T@new_M@G-M) < 5e-12
    return {'N':N,'k':k,'role_charges':q,'multiplicities':n,'role_F':f,
            'anomalies':anomalies(q,n,f),'matching_flow':x,'matching_pairs':pairs,
            'epsilon':str(Fraction(1,denominator)), 'row_sum_bound':degree,
            'exact_mass_singular_lower_bound':str(bound),
            'integer_mass_numerator_determinant':str(det),
            'exact_F_transpose_M_plus_MF_zero':True,
            'actual_all_allowed_role_blocks_nonzero':True,
            'zero_matching_bridges_activated':[i for i,v in enumerate(x) if not v],
            'finite_F_and_gauge_transform_residual_below_5e_minus_12':True,
            'full_flavor_irreducibility_claimed':False}


# Integer polynomial differential forms on R^4; keys are (basis wedge, monomial).
# This verifies a nonconstant local example of the general differential identity.
def clean(form):
    return {key:value for key,value in form.items() if value}


def form_add(*forms):
    out = {}
    for form in forms:
        for key,value in form.items():
            out[key] = out.get(key,0)+value
    return clean(out)


def form_scale(form,c):
    return clean({key:c*value for key,value in form.items()})


def wedge(a,b):
    out = {}
    for (u,p),x in a.items():
        for (v,q),y in b.items():
            if set(u)&set(v):
                continue
            sign = (-1)**sum(i>j for i in u for j in v)
            basis = tuple(sorted(u+v))
            powers = tuple(i+j for i,j in zip(p,q))
            key = (basis,powers)
            out[key] = out.get(key,0)+sign*x*y
    return clean(out)


def exterior_d(form):
    out = {}
    for (basis,powers),value in form.items():
        for i,power in enumerate(powers):
            if power:
                lower = list(powers)
                lower[i] -= 1
                term = wedge({((i,),tuple(lower)):value*power},
                             {(basis,(0,0,0,0)):1})
                out = form_add(out,term)
    return out


def scalar(powers,coefficient=1):
    return {((),tuple(powers)):coefficient}


def exterior_certificate():
    zero = (0,0,0,0)
    A = {((1,),(1,0,0,0)):1, ((3,),(0,0,1,0)):2}
    theta = scalar((1,0,1,0))
    beta = scalar((0,1,0,1))
    alpha = scalar((2,0,0,1))
    W = form_add(A,form_scale(exterior_d(theta),-1))
    FA = exterior_d(A)
    assert exterior_d(FA) == {}
    assert exterior_d(W) == FA
    W_new = form_add(A,exterior_d(alpha),
                     form_scale(exterior_d(form_add(theta,alpha)),-1))
    assert W_new == W
    FF = wedge(FA,FA)
    assert FF == {((0,1,2,3),zero):4}
    current = wedge(W,FA)
    assert exterior_d(current) == FF
    variation = wedge(exterior_d(beta),current)
    boundary = exterior_d(wedge(beta,current))
    assert form_add(variation,form_scale(boundary,-1),wedge(beta,FF)) == {}
    assert variation and boundary and wedge(beta,FF)
    return {'nonconstant_polynomial_example':True,'F_A_wedge_F_A_volume_coefficient':4,
            'W_invariant_under_A_plus_dalpha_theta_plus_alpha':True,
            'd_W_wedge_F_A_equals_F_A_wedge_F_A':True,
            'delta_counterterm_equals_boundary_minus_beta_F_A_squared':True,
            'global_extension_or_S_zero_regularity_claimed':False}


def payload():
    templates = {}
    for k,c in [(1,1),(2,1),(7,2)]:
        x = [c*v for v in [1,3,3,1]]
        cert = mass_certificate(x,k,4*(16*c))
        assert all(v==0 for v in cert['anomalies'].values())
        assert Fraction(cert['exact_mass_singular_lower_bound']) >= Fraction(3,4)
        templates[f'k{k}_c{c}'] = {'polynomial':polynomial_certificate(x),'matrix':cert}
    deletions = {}
    # Keep the accepted 1301 mass coefficients, not a replacement apparatus.
    for name,k,x,den in [('N12',2,[2,2,0,1,1],16),('N20',1,[3,5,2],80)]:
        cert = mass_certificate(x,k,den)
        assert cert['anomalies']['g2F'] == 2
        assert all(v==0 for key,v in cert['anomalies'].items() if key!='g2F')
        deletions[name] = {'polynomial':polynomial_certificate(x),'matrix':cert,
                          'old_gauge_model_invalidated':False}
    # Same count, different role spectrum; the quotient need not be nonnegative.
    nonunique = polynomial_certificate([1,2,1,1,2,1])
    assert nonunique['division_quotient'] == [1,-1,1]
    assert nonunique['total_Weyl_fields'] == 16
    assert all(p == [0] for p in nonunique['all_k_coefficients_low_to_high'].values())
    # A non-anomaly-free flow checks all polynomial coefficients, including M0!=0.
    general = polynomial_certificate([2,1,3])
    assert general['moments'][0] != 0 and not general['divisible_by_one_plus_z_cubed']
    sources = ['proof.md','NEXT.md','../research_note_1302.md',
               '../_shared/round1302_admission.md','../1301/proof.md',
               '../1301/results.json',
               '../../archive_1009_1043/research_note_1022.md',
               '../../archive_1009_1043/research_note_1034.md',
               '../../archive_1009_1043/research_note_1036.md',
               '../../archive_1046_1062/research_note_1060.md']
    return {'status':'PASS_AUTHOR_LOCAL_IDENTITY_CLASSIFICATION',
            'all_N_theorem':'N=16c, positive integer c; analytic, not finite scan',
            'positive_templates':templates,'deletion_controls':deletions,
            'general_coefficient_diagnostic':general,
            'same_N_nonunique_spectrum_diagnostic':nonunique,
            'local_compensation_form_check':exterior_certificate(),
            'new_axioms':0, 'physical_SM_generation_interpretation':False,
            'full_global_anomalies_claimed':False,'actual_instrument_O_claimed':False,
            'reference_helicity_compensation_excluded_from_original_role':True,
            'source_bindings':{p:sha256((HERE/p).read_bytes()).hexdigest() for p in sources},
            'check_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args = parser.parse_args()
    result = payload()
    if args.write:
        with OUT.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('CREATED',result['status'])
    else:
        expected = json.loads(OUT.read_text(encoding='utf-8'))
        if expected != result:
            raise SystemExit('MISMATCH: diagnostics or source bindings changed')
        print('READ_ONLY_MATCH',result['status'])


if __name__ == '__main__':
    main()
