"""1281 independent audit; no imports from author or historical checkers.

Exact Fraction real blocks audit the algebra. Complex Gram diagnostics use
separately generated data. Default only reads and compares saved results.
--write creates independent_results.json exclusively.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]

def fm(rows):
    return [[F(x) for x in row] for row in rows]

def ident(n):
    return [[F(int(i == j)) for j in range(n)] for i in range(n)]

def zero(n, m):
    return [[F(0) for _ in range(m)] for _ in range(n)]

def trans(a):
    return [list(x) for x in zip(*a)]

def plus(a, b):
    return [[x+y for x,y in zip(u,v)] for u,v in zip(a,b)]

def minus(a, b):
    return [[x-y for x,y in zip(u,v)] for u,v in zip(a,b)]

def scale(q, a):
    return [[q*x for x in row] for row in a]

def prod(a, b):
    bt = trans(b)
    return [[sum((x*y for x,y in zip(row,col)), F(0)) for col in bt] for row in a]

def hcat(a, b):
    return [x+y for x,y in zip(a,b)]

def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))

def infinity(a):
    return max(sum(abs(x) for x in row) for row in a)

def frobenius_sq(a):
    return sum(x*x for row in a for x in row)

def inv2(a):
    det = a[0][0]*a[1][1]-a[0][1]*a[1][0]
    assert det != 0
    return scale(1/det, [[a[1][1],-a[0][1]],[-a[1][0],a[0][0]]])

def psd2(a):
    assert a == trans(a)
    return (a[0][0] >= 0 and a[1][1] >= 0
            and a[0][0]*a[1][1]-a[0][1]*a[1][0] >= 0)

def real_block(a, c, k):
    # Unitary change diag(iI,I) converts the original complex H into this block.
    return hcat(a,scale(-1,c)) + hcat(scale(-1,trans(c)),k)

def frac(q):
    return str(q)

def exact_audit():
    i = ident(2)
    a = fm([[2,F(1,3)],[F(1,3),1]])
    r = fm([[F(1,2),F(2,3)],[F(2,3),-1]])
    q = plus(i,prod(r,r)); y = r+i
    assert prod(trans(y),y) == q
    assert psd2(a)
    ar = prod(a,r)
    assert ar != trans(ar)  # Genuine noncommutation; C need not be Hermitian.
    k = prod(r,ar); h = real_block(a,ar,k)
    f = hcat(i,scale(-1,r))
    assert h == prod(trans(f),prod(a,f))
    assert prod(h,y) == zero(4,2)

    sa = fm([[1,0],[0,0]])
    sc = prod(sa,r); sk = prod(r,sc)
    sh = real_block(sa,sc,sk)
    assert sh == prod(trans(f),prod(sa,f))
    assert prod(sh,y) == zero(4,2)

    e0 = fm([[F(1,7),-F(2,9)],[F(3,11),F(1,13)]])
    m0 = fm([[F(1,17),F(1,29)],[F(1,29),F(1,19)]])
    assert psd2(m0)
    certified = []
    for t in (F(1,10),F(1,3),F(1),F(2)):
        e = scale(t,e0); m = scale(t*t,m0)
        c = plus(ar,e)
        k = plus(prod(trans(c),prod(inv2(a),c)),m)
        h = real_block(a,c,k)
        d = prod(trans(y),prod(h,y))
        expected = plus(m,prod(trans(e),prod(inv2(a),e)))
        assert d == expected
        assert d == minus(minus(minus(k,prod(r,ar)),prod(r,e)),prod(trans(e),r))
        assert psd2(d) and psd2(minus(d,m))
        # rho=Q/TrQ, mu=1/TrQ; epsilon=TrD/TrQ exactly (no square roots).
        eps = trace(d)/trace(q)
        mu = 1/trace(q)
        q_bound = infinity(q); a_bound = infinity(a); r_bound = infinity(r)
        tau = eps*q_bound/mu; delta_sq = a_bound*tau
        assert psd2(minus(scale(tau,i),d))
        # Conservative rational norm certificates, independent of eigensolvers.
        assert frobenius_sq(e) <= delta_sq
        kgap = minus(k,prod(r,ar)); k_bound = infinity(kgap)
        assert max(F(0),k_bound-tau)**2 <= 4*r_bound*r_bound*delta_sq
        assert infinity(m) <= tau
        certified.append(dict(t=frac(t),epsilon=frac(eps),mu=frac(mu),
            tau_upper=frac(tau),delta_squared_upper=frac(delta_sq),
            C_frobenius_squared=frac(frobenius_sq(e)),K_infinity=frac(k_bound)))

    # The sqrt order is exact in the fixed role, for rational nonzero t.
    x = fm([[0,1],[1,0]]); z = fm([[1,0],[0,-1]])
    yx = x+i
    sqrt_rows=[]
    for t in (F(1,3),F(1,9),F(1,27)):
        c=plus(x,scale(t,z)); k=prod(c,c)
        h=real_block(i,c,k)
        d=prod(trans(yx),prod(h,yx))
        assert d == scale(t*t,i)
        assert prod(minus(c,x),minus(c,x)) == scale(t*t,i)
        sqrt_rows.append(dict(t=frac(t),epsilon=frac(t*t/2),
                              squared_C_gap=frac(t*t)))

    # Exact deletion witnesses, with data different from the author's rank test.
    hbad=real_block(i,zero(2,2),scale(-1,i))
    assert prod(trans(yx),prod(hbad,yx)) == zero(2,2)
    assert prod(hbad,yx) != zero(4,2) and hbad[2][2] == -1

    rr=fm([[1,0],[0,3]]); cc=fm([[1,0],[0,4]])
    yy=rr+i; hh=real_block(i,cc,prod(cc,cc))
    dd=prod(trans(yy),prod(hh,yy))
    assert dd == fm([[0,0],[0,1]])
    # Q=diag(2,10): rank-one rho on first state sees zero, second costs 1/10.
    missed=dd[1][1]/(1+rr[1][1]**2)
    assert missed == F(1,10)

    # Fixed nonzero baseline is not the stipulated zero positive comparison cost.
    offset=F(1,7)
    ko=plus(i,scale(2*offset,i))
    do=prod(trans(yx),prod(real_block(i,x,ko),yx))
    assert do == scale(2*offset,i)  # N=offset I, all input energies equal.
    # Re-fitting R to C lets incompatible references both pass.
    hz=real_block(i,z,i); yz=z+i
    assert prod(hz,yz)==zero(4,2)
    cross=prod(trans(yx),prod(hz,yx))
    assert cross == scale(2,i)  # N=I on fixed X code.

    # Reconstruct the author's chosen positive mismatch without using its code.
    c=plus(x,scale(F(1,4),z)); k=scale(F(9,8),i)
    m=minus(k,prod(trans(c),c))
    d=prod(trans(yx),prod(real_block(i,c,k),yx))
    assert m==scale(F(1,16),i) and d==scale(F(1,8),i)
    return dict(ar_nonhermitian=True,singular_A_kernel_exact=True,
        rational_finite_cases=certified,square_root_order=sqrt_rows,
        deletion=dict(negative_quadratic_witness="-1",
            rank_one_missed_cost=frac(missed),refitted_reference_fixed_code_cost="1",
            constant_nonzero_baseline=frac(offset)),
        author_positive_mismatch=dict(schur_remainder="I/16",energy="1/16"))

def spectral_norm(a):
    return float(np.linalg.norm(a,2))

def encode(r):
    q=np.eye(len(r))+r@r
    lam,v=np.linalg.eigh(q)
    return np.vstack((1j*r,np.eye(len(r))))@((v*(lam**-.5))@v.conj().T)

def numerical_audit():
    rng=np.random.default_rng(202610171281)
    residual=0.0; c_ratios=[]; k_ratios=[]; min_eig=1.0
    cases=0; singular_cases=0; zero_A_cases=0
    def random_mat(n,m):
        return rng.normal(size=(n,m))+1j*rng.normal(size=(n,m))
    for d in (2,3,5):
        for rank in (0,d-1,d):
            for sample in range(3):
                raw=random_mat(d,d); r=(raw+raw.conj().T)/5
                # H derives from a rectangular Gram matrix, including rank-deficient A.
                l=random_mat(rank,d)/2
                f=random_mat(rank,d)/11
                mraw=random_mat(d,d)/13
                m=mraw.conj().T@mraw
                a=l.conj().T@l
                w=l@r+f
                c=l.conj().T@w
                k=w.conj().T@w+m
                h=np.block([[a,-1j*c],[1j*c.conj().T,k]])
                gram=np.hstack((l,-1j*w))
                correction=np.zeros_like(h); correction[d:,d:]=m
                assert spectral_norm(h-(gram.conj().T@gram+correction)) < 1e-11
                b=encode(r); y=np.vstack((1j*r,np.eye(d)))
                dd=y.conj().T@h@y
                err=spectral_norm(dd-(f.conj().T@f+m))
                residual=max(residual,err,spectral_norm(b.conj().T@b-np.eye(d)))
                pp=random_mat(d,d); rho=pp.conj().T@pp+np.eye(d)
                rho/=np.trace(rho)
                mu=float(np.linalg.eigvalsh(rho)[0])
                n=b.conj().T@h@b
                eps=float(np.trace(rho@n).real)
                tau=eps*spectral_norm(np.eye(d)+r@r)/mu
                delta=np.sqrt(spectral_norm(a)*tau)
                e=c-a@r; kg=k-r@a@r
                cn=spectral_norm(e); kn=spectral_norm(kg)
                assert cn <= delta+1e-10
                assert kn <= tau+2*spectral_norm(r)*delta+1e-10
                assert np.linalg.eigvalsh(tau*np.eye(d)-dd)[0] >= -1e-10
                he=float(np.linalg.eigvalsh(h)[0]); min_eig=min(min_eig,he)
                assert he >= -1e-10
                if delta>1e-14:c_ratios.append(cn/delta)
                k_ratios.append(kn/(tau+2*spectral_norm(r)*delta))
                cases+=1; singular_cases+=rank<d; zero_A_cases+=rank==0

    # A noncommuting exact model, with an explicit partner and actual evolution.
    a=np.array([[2.,1/3],[1/3,1.]])
    r=np.array([[.5,2/3],[2/3,-1.]])
    c=a@r; k=r@a@r
    h=np.block([[a,-1j*c],[1j*c.conj().T,k]])
    b=encode(r)
    bell=np.array([1,0,0,1],complex)/np.sqrt(2)
    initial=np.kron(b,np.eye(2))@bell
    ev,v=np.linalg.eigh(h)
    u=(v*np.exp(-1j*.413*ev))@v.conj().T
    final=np.kron(u,np.eye(2))@initial
    reshaped=initial.reshape(4,2)
    partner=reshaped.conj().T@reshaped
    partner_residual=spectral_norm(partner-np.eye(2)/2)
    evolution_residual=float(np.linalg.norm(final-initial))
    assert partner_residual<1e-12 and evolution_residual<1e-12
    return dict(cases=cases,singular_A_cases=singular_cases,zero_A_cases=zero_A_cases,
        largest_identity_residual=float(format(residual,'.12g')),
        largest_C_bound_ratio=round(max(c_ratios),12),
        largest_K_bound_ratio=round(max(k_ratios),12),
        smallest_H_eigenvalue=float(format(min_eig,'.12g')),
        partner_marginal_residual=float(format(partner_residual,'.12g')),
        exact_code_evolution_residual=float(format(evolution_residual,'.12g')),
        diagnostics_are_not_universal_proof=True)

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def calculate():
    exact=exact_audit()
    numerical=numerical_audit()
    # Bind and validate saved author provenance; never import or execute author calculation.
    saved=json.loads((HERE/'results.json').read_text('utf-8'))
    checked=0
    for rel,want in saved['source_sha256'].items():
        assert sha(RESEARCH/rel)==want, 'Author source hash mismatch: '+rel
        checked+=1
    assert saved['positive_mismatch']['exact_calibration_energy']==exact['author_positive_mismatch']['energy']
    assert saved['new_axioms']==0
    bindings=[HERE/'proof.md',HERE/'check.py',HERE/'results.json',Path(__file__).resolve()]
    return dict(round=1281,status='PASS_WITH_SCOPE',independent_groups=1,new_axioms=0,
        implementation='Independent Fraction realification and rectangular complex Gram matrices; no author imports',
        exact_fraction=exact,complex_diagnostics=numerical,author_source_bindings_verified=checked,
        scope=dict(fixed_finite_comparison_class=True,arbitrary_partner_encoding=True,
            singular_A_allowed=True,positive_zero_fixed_independently=True,
            A_R_or_gauge_group_selected=False,spacetime_field_or_Gauss_dynamics_constructed=False,
            full_cognitive_axiom_model_certified=False,zero_cost_is_zero_total_resources=False),
        source_sha256={str(p.relative_to(RESEARCH)).replace('\\','/'):sha(p) for p in bindings})

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true')
    args=parser.parse_args()
    payload=calculate(); target=HERE/'independent_results.json'
    if args.write:
        with target.open('x',encoding='utf-8',newline='\n') as f:
            json.dump(payload,f,ensure_ascii=False,indent=2);f.write('\n')
        mode='CREATED_ONCE'
    else:
        assert json.loads(target.read_text('utf-8'))==payload,'Independent replay/source mismatch'
        mode='READ_ONLY_MATCH'
    print(json.dumps(dict(status=payload['status'],mode=mode,round=1281,new_axioms=0,
        rational_cases=len(payload['exact_fraction']['rational_finite_cases']),
        complex_cases=payload['complex_diagnostics']['cases'],
        author_bindings_verified=payload['author_source_bindings_verified']),ensure_ascii=False))

if __name__=='__main__':
    main()
