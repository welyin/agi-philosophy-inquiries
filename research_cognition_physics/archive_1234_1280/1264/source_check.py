"""Round 1264: one source-author group, read-only reproduction by default.
Finite analytic majorants only. No PDE/Fourier quadrature/QFT truncation.
"""
from fractions import Fraction as F
from math import comb, factorial, isqrt
import json
import argparse
import hashlib
from pathlib import Path

HERE=Path(__file__).resolve().parent

# eta^(n)(t)=exp(-1/t) P_n(1/t), P_(n+1)=v^2(P_n-P_n').
P = {0: 1}
eta = []
for n in range(7):
    eta.append(sum((F(abs(c)) * (F(k, 2)**k if k else 1)
                    for k, c in P.items()), F(0)))
    nxt = {}
    for k, c in P.items():
        nxt[k + 2] = nxt.get(k + 2, 0) + c
        if k:
            nxt[k + 1] = nxt.get(k + 1, 0) - k*c
    P = {k:c for k,c in nxt.items() if c}

# D rho=u, D>=1/2, |D^(k)|<=2 eta[k].
rho = [F(1)]
for n in range(1, 7):
    rho.append(2*eta[n] + 4*sum(
        (comb(n,k)*eta[k]*rho[n-k] for k in range(1,n+1)), F(0)))

def multi(r):
    return [(a,b,r-a-b) for a in range(r+1) for b in range(r-a+1)]

# Exact chain-rule polynomials:
# derivative^alpha F(x)=sum c*x^power*rho^(k)(|x|^2).
polys = {(0,0,0): {(0,0,0,0): 1}}
for r in range(1,7):
    for alpha in multi(r):
        i = next(i for i,x in enumerate(alpha) if x)
        prev=list(alpha); prev[i]-=1
        out={}
        for (k,*powers),c in polys[tuple(prev)].items():
            if powers[i]:
                ex=powers.copy(); ex[i]-=1
                key=(k,*ex)
                out[key]=out.get(key,0)+c*powers[i]
            ex=powers.copy(); ex[i]+=1
            key=(k+1,*ex)
            out[key]=out.get(key,0)+2*c
        polys[alpha]={k:c for k,c in out.items() if c}
bounds = {a:sum((abs(c)*2**sum(key[1:])*rho[key[0]]
                 for key,c in p.items()), F(0)) for a,p in polys.items()}

def sobolev_C(j):
    sj=F(0)
    for r in range(j+3):
        for a in multi(r):
            multin=factorial(r)//(factorial(a[0])*factorial(a[1])*factorial(a[2]))
            sj += comb(j+2,r)*multin*64*bounds[a]**2
    # Cj^2 <= Sj/(8*pi) < Sj/24; x+1>=sqrt(x).
    return sj/F(24)+1

C2,C4=sobolev_C(2),sobolev_C(4)
def ceil_root(q,n):
    k=isqrt(q.numerator//q.denominator)
    if n==4:
        k=isqrt(k)
    assert n in (2,4)
    while k**n*q.denominator < q.numerator:
        k+=1
    assert k**n >= q
    return k

L=max(10000,ceil_root(F(625,52)*C2,2),
      ceil_root(F(125,281216)*C4,4))
R=2*L
D=F(5,208)*C2/R**2+F(1,281216)*C4/R**4
assert D<=F(1,1000)

# The same C2/L budget supplies Hessian and gradient controls.
mmin,mmax=F(26),F(32)
exp2=F(101,100)**2
h=F(13,625)
grad_squared=2*h
assert C2/R**2<=h
Hprime_upper=mmax*exp2+(h+2*F(1,100)*grad_squared)/(2*mmin)
Hsecond_upper=2*mmax*exp2+grad_squared/mmin
Hsecond_negative_magnitude=F(33**2,26)
assert Hprime_upper<33
assert Hsecond_upper<66
assert Hsecond_negative_magnitude<42
r=F(1,10000)
packet=F(36,13)/L**2
source=4*(r*D+33*r*r)
contrast=F(9,16)*r-packet-source
assert contrast>r/2
def q(x): return str(x.numerator)+"/"+str(x.denominator)
payload={
    "status":"PASS_AUTHOR_WITH_SCOPE",
    "round":1264,"source_author_groups":1,
    "scope":"finite same-source static background comparison; not generated geometry",
    "profile":"rho(|x|^2/R^2), rho=eta(4-s)/(eta(4-s)+eta(s-1))",
    "C2_upper":q(C2),"C4_upper":q(C4),
    "finite_L":str(L),"finite_R":str(R),
    "source_difference_upper":q(D),
    "Hprime_upper":q(Hprime_upper),
    "Hsecond_upper":q(Hsecond_upper),
    "Hsecond_negative_magnitude":q(Hsecond_negative_magnitude),
    "declared_taylor_coefficient":33,
    "source_probability_upper":q(source),
    "packet_probability_upper":q(packet),
    "contrast_lower":q(contrast),
    "declared_contrast":q(r/2),
    "display_only":{"D":float(D),"contrast":float(contrast),
                    "Hsecond_upper":float(Hsecond_upper)}
}
bindings={
    "source_code":HERE/"source_check.py",
    "source_proof":HERE/"source_proof.md",
    "round1262_proof":HERE.parent/"1262/proof.md",
    "round1253_material":HERE.parent/"1253/proof.md",
    "round1186_contact":HERE.parent.parent/"archive_1161_1200/1186/proof.md",
}
payload["source_bindings"]={name:hashlib.sha256(path.read_bytes()).hexdigest()
                            for name,path in bindings.items()}
ap=argparse.ArgumentParser()
ap.add_argument("--write-results",action="store_true")
args=ap.parse_args()
out=HERE/"source_results.json"
if args.write_results:
    with out.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(payload,f,ensure_ascii=False,indent=2)
        f.write("\n")
else:
    saved=json.loads(out.read_text(encoding="utf-8"))
    if saved!=payload:
        raise SystemExit("FAIL: readonly payload mismatch")
print(json.dumps({"status":payload["status"],
                  "source_author_groups":payload["source_author_groups"],
                  "mode":"exclusive_first_write" if args.write_results else "readonly_match",
                  "L":str(L),"R":str(R),"display_only":payload["display_only"]},
                 ensure_ascii=False,indent=2))

