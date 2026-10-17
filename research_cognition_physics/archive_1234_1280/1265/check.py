"""1265: exact static conformal tensor and one finite constraint budget.
One author group. Defaults to read-only result matching; no PDE or QFT simulation.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent

# Sparse polynomials: monomial is a sorted tuple of jet variable names.
def c(a):
    a=F(a)
    return {():a} if a else {}
def var(name): return {(name,):F(1)}
def add(*ps):
    out={}
    for p in ps:
        for k,v in p.items(): out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def scale(a,p): return {k:F(a)*v for k,v in p.items() if a*v}
def mul(p,q):
    out={}
    for a,x in p.items():
        for b,y in q.items():
            k=tuple(sorted(a+b))
            out[k]=out.get(k,F(0))+x*y
    return {k:v for k,v in out.items() if v}
eta=[-1,1,1,1]
grad=[{}]+[var("g"+str(i)) for i in range(1,4)]
hess=[[{} if (i==0 or j==0) else var("h"+str(min(i,j))+str(max(i,j)))
       for j in range(4)] for i in range(4)]
def Gamma(a,u,v):
    return add(grad[v] if a==u else {},
               grad[u] if a==v else {},
               scale(-eta[u]*eta[a],grad[a]) if u==v else {})
def dGamma(k,a,u,v):
    return add(hess[v][k] if a==u else {},
               hess[u][k] if a==v else {},
               scale(-eta[u]*eta[a],hess[a][k]) if u==v else {})
Ricci=[]
for u in range(4):
    row=[]
    for v in range(4):
        p={}
        for a in range(4):
            p=add(p,dGamma(a,a,u,v),scale(-1,dGamma(v,a,u,a)))
            for b in range(4):
                p=add(p,mul(Gamma(a,a,b),Gamma(b,u,v)),
                      scale(-1,mul(Gamma(a,v,b),Gamma(b,u,a))))
        row.append(p)
    Ricci.append(row)
Rflat=add(*(scale(eta[i],Ricci[i][i]) for i in range(4)))
Einstein=[[add(Ricci[i][j],scale(F(-eta[i],2),Rflat) if i==j else {})
           for j in range(4)] for i in range(4)]
lap=add(*(hess[i][i] for i in range(1,4)))
grad2=add(*(mul(grad[i],grad[i]) for i in range(1,4)))
target=[[{} for j in range(4)] for i in range(4)]
target[0][0]=add(scale(-2,lap),scale(-1,grad2))
for i in range(1,4):
    for j in range(1,4):
        target[i][j]=add(scale(-2,hess[i][j]),scale(2,mul(grad[i],grad[j])),
                        add(scale(2,lap),grad2) if i==j else {})
assert Einstein==target
assert Rflat==add(scale(-6,lap),scale(-6,grad2))
spatial_trace=add(*(Einstein[i][i] for i in range(1,4)))
assert add(spatial_trace,scale(2,Einstein[0][0]))==scale(3,grad2)

old=HERE.parent/"1264"
record=json.loads((old/"results.json").read_text(encoding="utf-8"))
source=json.loads((old/"source_results.json").read_text(encoding="utf-8"))
R=int(record["finite_R"])
assert R==int(source["finite_R"])==38201213950
r=F(record["rational_authoritative"]["r"])
assert r==F(1,10000)

# Radial CS: int_R^(2R) dr/r^2=1/(2R);
# I >= 8*pi*R, Z <= (256*pi/3)*R^3.
cs_integral=F(1,2*R)
I_over_pi=F(4)/cs_integral
Z_over_pi_upper=F(256,3)*R**3
linear_coeff=2*I_over_pi/Z_over_pi_upper
assert linear_coeff==F(3,16*R**2)
lambda_lower=r/(1+r)  # log(1+r) strictly greater for r>0.
d=linear_coeff*lambda_lower
assert d==F(3,233516588884543995193640000)
# An explicit tolerance on either side of the necessary gate.
rejected_tolerance=d/2
not_rejected_by_this_bound=2*d
assert rejected_tolerance<d<not_rejected_by_this_bound
def q(x): return f"{x.numerator}/{x.denominator}"
def serial(p): return {"*".join(k) if k else "1":q(v) for k,v in sorted(p.items())}
paths={
    "code":HERE/"check.py","proof":HERE/"proof.md",
    "round1264_source_proof":old/"source_proof.md",
    "round1264_record_proof":old/"record_proof.md",
    "round1264_source_results":old/"source_results.json",
    "round1264_record_results":old/"results.json",
    "round1194_source_scope":HERE.parent.parent/"archive_1161_1200/1194/proof.md",
    "round1245_constraint_identity":HERE.parent/"1245/proof.md",
}
payload={
    "status":"PASS_AUTHOR_WITH_SCOPE","round":1265,"author_groups":1,
    "group":"same_profile_full_tensor_and_finite_positive_source_constraint",
    "tensor_identity":"all 16 components, exact polynomial coefficients",
    "G00_polynomial":serial(Einstein[0][0]),
    "spatial_trace_plus_2G00":serial(add(spatial_trace,scale(2,Einstein[0][0]))),
    "finite_R":str(R),
    "rational_authoritative":{
        "r":q(r),"lambda_strict_lower":q(lambda_lower),
        "radial_cs_integral":q(cs_integral),"I_over_pi_lower":q(I_over_pi),
        "Z_over_pi_upper":q(Z_over_pi_upper),"linear_coefficient":q(linear_coeff),
        "strict_abs_constraint_residual_lower":q(d),
        "example_rejected_tolerance":q(rejected_tolerance),
        "example_tolerance_not_rejected_by_this_bound":q(not_rejected_by_this_bound)
    },
    "scope":[
        "prescribed static compact conformal metric, no inner boundary",
        "kappa positive and effective cosmological term zero",
        "minimal positive-mass one-particle relative-vacuum source with T00 nonnegative",
        "extra vacuum or negative-source sectors not silently included",
        "exact extension excluded; finite tolerance excluded only below the stated bound",
        "not a probability bound, full quantum constraint or effective GR no-go"
    ],
    "display_only":{"constraint_lower":float(d)},
    "bindings_sha256":{name:hashlib.sha256(path.read_bytes()).hexdigest()
                       for name,path in paths.items()}
}
ap=argparse.ArgumentParser()
ap.add_argument("--write-results",action="store_true")
args=ap.parse_args()
out=HERE/"results.json"
if args.write_results:
    with out.open("x",encoding="utf-8",newline="\n") as f:
        json.dump(payload,f,ensure_ascii=False,indent=2)
        f.write("\n")
else:
    if json.loads(out.read_text(encoding="utf-8"))!=payload:
        raise SystemExit("FAIL: readonly result mismatch")
print(json.dumps({"status":payload["status"],"author_groups":1,
    "mode":"exclusive_first_write" if args.write_results else "readonly_match",
    "exact_tensor_components":16,"strict_residual_lower":q(d),
    "display_only":payload["display_only"]},ensure_ascii=False,indent=2))
