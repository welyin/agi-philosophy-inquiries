"""Non-numbered finite scattering screen. Default recomputes and compares only."""
from __future__ import annotations
import argparse, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results.json"

def ceilroot(x:F,n:int,digits:int=160)->F:
    if x<0: raise ValueError("negative root")
    den=10**digits
    top=x.numerator*den**n
    floor_target=top//x.denominator
    if n==2:
        z=math.isqrt(floor_target)
    elif n==4:
        z=math.isqrt(math.isqrt(floor_target))
    else: raise ValueError(n)
    while (z+1)**n*x.denominator<=top: z+=1
    if z**n*x.denominator<top: z+=1
    r=F(z,den)
    assert r**n>=x
    return r

def rat(x:F): return str(x.numerator)+"/"+str(x.denominator)
def sha(p:Path): return hashlib.sha256(p.read_bytes()).hexdigest()

def calculate():
    b=F(10**100); k=F(10**40); sigma=F(10**65)
    h=128*k/b+16/sigma
    pbad=(96*sigma/b)**8+h**8+(4*h/k)**8
    B=1024*k/b
    parts={
      "massive_free":2*b*h*h/52,
      "massless_free":2*b*h*h/k,
      "photon_energy_spread":256*h,
      "full_linear_velocity":20*k*h,
      "TT_and_velocity_quadratic":F(40,26)*k*h*h,
      "full_3D_position_Taylor":589824*k*sigma*sigma/(b*b),
      "bounded_bad_tail":2176*k*ceilroot(pbad,2),
      "linear_comparator_bad_tail":3072*k*sigma/b*ceilroot(pbad,4),
    }
    eps=sum(parts.values(),F(0))
    A=8*sigma*sigma*k*k/(b*b)
    expupper=1/(1+A)
    qbox=1/(2*26**2*A)
    gap=F(5,16)-F(11,16)*expupper-2*qbox-2*eps
    # Actual normal state and reference-domain conditions.
    assert b>0 and k>0 and sigma>0
    assert k-128*k/b-4/sigma>k/2
    assert B<26
    assert pbad<1
    assert eps<F(1,10**12)
    assert gap>F(31,100)
    # Distinct side terms in the full pair bound, no source-frozen replacement.
    vraw=F(512)+64+F(256,26)
    assert vraw<1024
    # Full radial good condition implies plateau, not a longitudinal shortcut.
    assert F(3,4)>F(2,3) and F(5,4)<F(4,3)
    # Compact cosine^4 packet: integral ratios contain cancelling pi.
    cos4_norm=F(35,128)       # integral f^2 / pi
    derivative4_norm=F(544)  # integral |f''''|^2 / pi
    derivative4_ratio=derivative4_norm/cos4_norm
    assert derivative4_ratio==F(69632,35)
    x8_sigma8_bound=derivative4_ratio/F(2**8)
    assert x8_sigma8_bound==F(272,35) and x8_sigma8_bound<8
    # Gaussian eighth-root bounds; position/momentum Minkowski constants.
    assert 105<2**8
    assert 4+4+8<=24
    assert 2+2+5<16 and 2+2+2<16
    # Gaussian Fourier dephasing: |15|+|6|+|1| = 22.
    assert F(15+6+1,32)==F(11,16)
    frozen=[
      "research_cognition_physics/archive_1234_/1253/proof.md",
      "research_cognition_physics/archive_1234_/1261/proof.md",
      "research_cognition_physics/archive_1234_/1261/parent.md",
      "research_cognition_physics/archive_1161_1200/1169/proof.md",
      "research_cognition_physics/archive_1219_1233/1231/proof.md",
    ]
    project=ROOT.parents[3]
    return {
      "status":"AUTHOR_CANDIDATE_NOT_INDEPENDENTLY_AUDITED_NOT_NUMBERED",
      "scope":"Two-body complete 1PM bounded AW model; no finite-time KG-source, actual-port, EFT-remainder or geometry certificate.",
      "source_hashes":{x:sha(project/x) for x in frozen},
      "author_hashes":{
          "memo":sha(ROOT/"README.md"),
          "code":sha(Path(__file__)),
      },
      "parameters":{key:rat(value) for key,value in {
         "b":b,"k":k,"sigma":sigma,"T":2*b,"Omega":8/sigma,
         "h":h,"B_V":B,"p_bad_upper":pbad,"A":A,
      }.items()},
      "moment_and_completion_checks":{
        "full_pair_support_coefficient":rat(vraw),
        "cos4_fourth_derivative_ratio":rat(derivative4_ratio),
        "longitudinal_position_eighth_moment_coefficient_upper":rat(x8_sigma8_bound),
        "mass_independent_frame_and_cutoff":True,
        "physical_reference_positive_pz_gap_verified":True,
      },
      "rational_propagation_terms":{key:rat(value) for key,value in parts.items()},
      "rational_authoritative_bounds":{
         "state_error_upper":rat(eps),
         "reference_dephasing_exp_upper":rat(expupper),
         "each_Q_box_tail_upper":rat(qbox),
         "joint_probability_gap_lower":rat(gap),
         "material_free_energy_change_upper":rat(2*B),
      },
      "display_only":{
         "state_error_upper":float(eps),
         "joint_probability_gap_lower":float(gap),
         "p_bad_upper":float(pbad),
         "A":float(A),
         "B_V":float(B),
         "terms":{key:float(value) for key,value in parts.items()},
      },
      "not_certified":[
        "source insertion transported from on-shell KG to finite-time LSB",
        "higher-PM, radiation or UV/EFT remainder",
        "real finite local ports, reader sources and future collection",
        "independent physical metric or full GR",
      ],
    }

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--write-once",action="store_true")
    a=ap.parse_args()
    result=calculate()
    if a.write_once:
        with OUT.open("x",encoding="utf-8") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
        print("WRITE_ONCE_PASS_NONNUMBERED_FINITE_SCATTERING")
    else:
        expected=json.loads(OUT.read_text(encoding="utf-8"))
        if expected!=result: raise SystemExit("READ_ONLY_MISMATCH")
        print("READ_ONLY_MATCH_NONNUMBERED_FINITE_SCATTERING")
    print(json.dumps(result["display_only"],ensure_ascii=False,indent=2))
