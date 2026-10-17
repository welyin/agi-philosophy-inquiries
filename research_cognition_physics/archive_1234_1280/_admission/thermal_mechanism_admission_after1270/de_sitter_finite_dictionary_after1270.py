"""Finite de Sitter null dictionary: exact rational Laurent-polynomial audit.
No field discretization, detector model, or nonlinear Einstein certificate.
Default is read-only; --write creates results exclusively.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, sys

HERE=Path(__file__).resolve()
ROOT=next(p for p in HERE.parents if (p/"research_cognition_physics").is_dir() and (p/"scripts").is_dir())
MD=HERE.with_suffix(".md")
OUT=HERE.with_suffix(".json")

def clean(p): return {j:F(c) for j,c in p.items() if c}
def add(*ps):
    z={}
    for p in ps:
        for j,c in p.items(): z[j]=z.get(j,F(0))+c
    return clean(z)
def scale(p,c): return clean({j:c*a for j,a in p.items()})
def mul(p,q):
    z={}
    for i,a in p.items():
        for j,b in q.items(): z[i+j]=z.get(i+j,F(0))+a*b
    return clean(z)
def power(p,n):
    z={0:F(1)}
    for _ in range(n): z=mul(z,p)
    return z
def diff(p): return clean({j-1:j*c for j,c in p.items()})
def ev(p,x): return sum((c*x**j for j,c in p.items()),F(0))
def primitive(p):
    assert p.get(-1,F(0))==0
    return clean({j+1:c/F(j+1) for j,c in p.items()})
def integral(p,a,b):
    # Exact rational + logarithm coefficient * log(abs(b/a)).
    rational=F(0); log=F(0)
    for j,c in p.items():
        if j==-1: log+=c
        else: rational+=c/F(j+1)*(b**(j+1)-a**(j+1))
    return rational,log
def pair_add(x,y): return x[0]+y[0],x[1]+y[1]
def pair_scale(x,c): return x[0]*c,x[1]*c
def pair_num(x): return [str(x[0]),str(x[1])]
def operator_L(p): return add(diff(diff(p)),scale(mul({-1:F(1)},diff(p)),2))
def edge(p,a,b): return ev(p,b)-ev(p,a)
def boundary(p,w,ell,a,b): return edge(add(p,scale(mul(w,diff(p)),1/(2*ell))),a,b)
def log_two_bounds(n=16):
    low=2*sum((F(1,3)**(2*j+1)/F(2*j+1) for j in range(n)),F(0))
    tail=2*F(1,3)**(2*n+1)/(F(2*n+1)*(1-F(1,9)))
    return low,low+tail
def pair_bounds_log_half(x):
    lo,hi=log_two_bounds()
    vals=(x[0]-x[1]*lo,x[0]-x[1]*hi)
    return min(vals),max(vals)

def compute():
    a,b,ell,xi=F(-2),F(-1),F(2),F(1,6)
    u={1:F(1)}; w={2:F(1),1:2*ell}
    g=mul(power({1:F(1),0:F(2)},3),power({1:F(1),0:F(1)},3))
    gf=mul(g,g); can=mul(diff(g),diff(g))
    imp=add(can,scale(operator_L(gf),-xi))
    # qhat = delta A/(8*pi*G); L qhat = -J_imp, area per unit solid angle H=1.
    raw=primitive(mul({-2:F(1)},primitive(scale(mul({2:F(1)},imp),-1))))
    c2=(ev(raw,b)-ev(raw,a))/(1/a-1/b)
    c1=-ev(raw,a)-c2/a
    q=add(raw,{0:c1,-1:c2})
    assert ev(q,a)==ev(q,b)==0
    assert operator_L(q)==scale(imp,-1)
    assert all(ev(diff(g) if j==1 else diff(diff(g)) if j==2 else g,x)==0 for j in range(3) for x in (a,b))
    # Physical Phi=g and flat phi=Omega*Phi=-2g/u on v=0.
    flat=scale(mul({-1:F(1)},g),-2)
    tf=add(mul(diff(flat),diff(flat)),scale(diff(diff(mul(flat,flat))),-xi))
    assert scale(mul({2:F(1)},tf),F(1,4))==imp
    kin=pair_scale(integral(mul(w,imp),a,b),-1/ell)
    kc=pair_scale(integral(mul(w,can),a,b),-1/ell)
    geometry=(2*boundary(q,w,ell,a,b),F(0))
    # S_xi/pi=-2 xi F_A; S_area/pi=2 qhat.
    sx=scale(gf,-2*xi)
    sw=add(scale(q,2),sx)
    assert kin==geometry
    assert kc==(boundary(sw,w,ell,a,b),F(0))
    assert kin==kc and kin[1]==0 and kin[0]>0
    # Nonzero-contact interior window of the SAME solution/wave.
    ai,bi=F(-7,4),F(-5,4)
    ki=pair_scale(integral(mul(w,imp),ai,bi),-1/ell)
    kci=pair_scale(integral(mul(w,can),ai,bi),-1/ell)
    contact=boundary(sx,w,ell,ai,bi)
    assert ki==pair_add(kci,(-contact,F(0)))
    assert contact!=0
    assert ki==(2*boundary(q,w,ell,ai,bi),F(0))
    # Literal printed (18) if T_uu is literally the flat covariant component.
    lit=pair_scale(integral(mul(w,tf),a,b),-2/ell)
    v=mul(w,{-2:F(1)})
    bnd=add(mul(v,diff(q)),mul(add({-1:F(2),-2:6*ell}),q))
    density={-2:F(2),-3:12*ell}
    lit_from_area=pair_scale(pair_add((edge(bnd,a,b),F(0)),integral(mul(density,q),a,b)),F(8)/ell)
    assert lit==lit_from_area
    # Printed (20)'s kernel naively applied to physical A on this finite segment.
    naive=pair_scale(integral(mul(w,diff(diff(q))),a,b),F(1)/(2*ell))
    naive_difference=pair_add(naive,pair_scale(kin,-1))
    ndlo,ndhi=pair_bounds_log_half(naive_difference)
    assert -F(3,2500)<ndlo<=ndhi<-F(1,1000)
    # Affine conversion lambda=-1/u: xi^lambda/pi=-(1+2ell/u)/ell.
    affine_xi=scale(mul(w,{-2:F(1)}),-1/ell)
    assert affine_xi=={0:-1/ell,-1:F(-2)}
    # Characteristic polynomial data and exact physical finite gap.
    return {
        "status":"PASS_FINITE_DICTIONARY_ONLY",
        "scope":"linear Einstein input; physical null area and coherent relative source; no actual detector/QFT implementation",
        "algebra_groups":1,
        "parameters":{"H":"1","G":"1","ell":str(ell),"a":str(a),"b":str(b),"xi":str(xi),"solid_angle":"1","field_amplitude":"eta; displayed fluxes divided by pi*eta^2"},
        "wave":"Phi(u)=(u+2)^3(u+1)^3 on [-2,-1], zero outside; C2",
        "qhat_coefficients":{str(j):str(c) for j,c in sorted(q.items())},
        "area_endpoint_qhat":[str(ev(q,a)),str(ev(q,b))],
        "area_endpoint_qhat_prime":[str(ev(diff(q),a)),str(ev(diff(q),b))],
        "K_improved_over_pi":pair_num(kin),
        "K_canonical_over_pi":pair_num(kc),
        "physical_area_boundary_over_pi":pair_num(geometry),
        "interior_window":[str(ai),str(bi)],
        "interior_K_improved_over_pi":pair_num(ki),
        "interior_K_canonical_over_pi":pair_num(kci),
        "interior_Wald_boundary_over_pi":str(contact),
        "literal_eq18_over_pi":pair_num(lit),
        "literal_eq18_correct_weighted_area_over_pi":pair_num(lit_from_area),
        "naive_eq20_physical_area_over_pi":pair_num(naive),
        "naive_minus_physical_lower":str(ndlo),
        "naive_minus_physical_upper":str(ndhi),
        "log_convention":"for main window: log(abs(b/a))=log(1/2); all entries exact",
        "identities":{"conformal_complete_stress":True,"nonaffine_raychaudhuri":True,"physical_current_area":True,"finite_Wald_contact":True,"literal18_weighted_area":True,"naive20_dictionary_rejected":True},
    }
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def output():
    r=compute()
    sources={
        "code":HERE,"memo":MD,
        "round318":ROOT/"research_cognition_physics/archive_301_341/research_note_318.md",
        "round1250":ROOT/"research_cognition_physics/archive_1234_/1250/proof.md",
        "previous_screen":HERE.parent/"geometric_thermal_identity_candidate_after1270.md",
    }
    r["source_bindings"]={k:{"path":str(p.relative_to(ROOT)).replace("\\","/"),"sha256":sha(p)} for k,p in sources.items()}
    return r
if __name__=="__main__":
    r=output()
    if sys.argv[1:]==["--write"]:
        with OUT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(r,f,ensure_ascii=False,indent=2,sort_keys=True);f.write("\n")
        print("FIRST_WRITE_PASS")
    elif not sys.argv[1:]:
        old=json.loads(OUT.read_text(encoding="utf-8"))
        assert old==r,"saved result mismatch"
        print("READ_ONLY_MATCH")
    else: raise SystemExit("use no argument or --write")
    print(r["K_improved_over_pi"], r["interior_Wald_boundary_over_pi"])

