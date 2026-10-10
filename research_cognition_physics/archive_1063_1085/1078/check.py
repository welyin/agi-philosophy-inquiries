"""1078: exact interval audit of the unchanged 1077 source family.
Prints JSON only. Uses frozen integer source/reader, no new model or controls.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,importlib.util,json
import numpy as np

BASE=Path(__file__).resolve().parent
def load_source():
    p=BASE.parent/"1077/check.py"
    spec=importlib.util.spec_from_file_location("frozen1077_source",p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    return m

def original_data():
    source=load_source();reader=source.load_reader()
    _,g,dg,gs,tail,e=source.rational_source(reader)
    *_,en,raw=reader.exact_reader()
    nr,ni,ns=raw
    br=[nr[:6,:6]+nr[6:,6:],nr[:6,6:]+nr[6:,:6],
        ni[6:,:6]-ni[:6,6:],nr[:6,:6]-nr[6:,6:]]
    bi=[ni[:6,:6]+ni[6:,6:],ni[:6,6:]+ni[6:,:6],
        nr[:6,6:]-nr[6:,:6],ni[:6,:6]-ni[6:,6:]]
    def contract(a,b,c,d):return F(int(np.sum(a*c.T-b*d.T)),gs*2*ns)
    v=[contract(a,b,*g) for a,b in zip(br,bi)]
    j=[[contract(a,b,*der) for der in dg] for a,b in zip(br,bi)]
    return v,j,e,en

def interval(center,radius):return (center-radius,center+radius)
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def neg(a):return (-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    p=[x*y for x in a for y in b]
    return min(p),max(p)
def scale(a,s):return mul(a,(F(s),F(s)))
def total(values):
    s=(F(0),F(0))
    for v in values:s=add(s,v)
    return s
def outward(a,den=10**12):
    l=a[0]*den;u=a[1]*den
    lo=F(l.numerator//l.denominator,den)
    hi=F(-((-u.numerator)//u.denominator),den)
    assert lo<=a[0]<=a[1]<=hi
    return [str(lo),str(hi)]

def run():
    values,j,e,en=original_data()
    tau=F(1,10**10);r0=F(1,246240000);m=(9,9,1)
    ev=en+(1+en)*e*(2+e)+18*tau
    assert ev<F(1,500000000)
    ej=[2*a*en+(1+en)*2*a*e*(2+e)+36*a*tau for a in m]
    bv=[interval(v,ev+38*r0) for v in values]
    jj=[[interval(x,ej[a]+76*m[a]*r0) for a,x in enumerate(row)] for row in j]
    ds=[scale(total(mul(bv[k],jj[k][a]) for k in (1,2,3)),2) for a in range(3)]
    # Columns (u,theta): determinant of D(b0, |b|²).
    minor=sub(mul(jj[0][1],ds[2]),mul(jj[0][2],ds[1]))
    assert F(-1,50000)<minor[0]<=minor[1]<F(-1,100000)
    coarse=outward(minor)
    assert F(coarse[1])<F(-1,100000)
    # Reuse 1077 global parameter embedding, including all cube boundary points.
    embedding=F(1,2000)-2052*r0
    assert embedding==F(59,120000)>0
    detcenter=j[0][1]*2*sum(values[k]*j[k][2] for k in (1,2,3))-j[0][2]*2*sum(values[k]*j[k][1] for k in (1,2,3))
    # The small control cube's entire absolute Bloch image is in one hemisphere.
    # This is an additional scope diagnostic, not a new dimension theorem.
    assert bv[3][0]>F(4,25) # 0.16
    hashes={}
    for rel in ("1077/check.py","1077/proof.md","1077/results.json","1076/check.py"):
        hashes[rel]=hashlib.sha256((BASE.parent/rel).read_bytes()).hexdigest()
    return {
      "round":1078,"status":"strict_certificate_passed",
      "source_family":"unchanged 1077 finite H-on W_tau(p), fixed preparation and fixed reader",
      "parameter_order":["t","u","theta"],"pulse_tau":str(tau),"reader_time":"7/10",
      "domain":"entire closed 480 control cube about (1/5,2/5,pi/2)",
      "control_radius":str(r0),
      "spectral_invariants":["c=b0","s=bx^2+by^2+bz^2"],
      "center_effect_value_error_upper":str(F(1,500000000)),
      "center_effect_value_error_float":float(ev),
      "center_jacobian_error_by_column_float":[float(x) for x in ej],
      "uniform_effect_boxes": [outward(v) for v in bv],
      "uniform_effect_derivative_boxes":[[outward(v) for v in row] for row in jj],
      "uniform_spectral_derivative_boxes":[[outward(v) for v in jj[0]],[outward(v) for v in ds]],
      "minor_columns":["u","theta"],
      "uniform_minor_enclosure":coarse,
      "uniform_minor_float":[float(v) for v in minor],
      "uniform_minor_strict_upper":"-1/100000",
      "center_minor_diagnostic":float(detcenter),
      "uniform_spectral_jacobian_rank":2,
      "ambient_local_isospectral_dimension":1,
      "within_closed_cube_dimension_upper":1,
      "smooth_extension_minimum_extra_parameter_if_old_source_retained":1,
      "reused_bloch_embedding_lower":str(embedding),
      "bloch_z_uniform_strict_lower":"4/25",
      "source_hashes":hashes,
      "scope":{"local_rank_two_isospectral_angular_surface_excluded":True,
        "noncommuting_continuous_covariant_action_inside_same_cube_excluded":True,
        "whole_orbit_truncation_is_not_the_new_result":True,
        "arbitrary_context_transport_excluded":False,
        "extension_count_requires_C1_manifold_and_embedded_old_source":True,
        "all_control_parameters_excluded":False,
        "all_cognition_or_spatial_dimension_refuted":False,
        "new_controller_or_dephasing":False,
        "empirical_results":0,"new_accepted_cognitive_axioms":0}
    }
if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,indent=2))
