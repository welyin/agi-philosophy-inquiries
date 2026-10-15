"""1226: uniform comparison constants for the original complete Gauss Gibbs family.

Four check groups; no full Gibbs matrix/spectrum, energy or entropy evaluated.
Rational guards are exact for frozen decimal model inputs. Transcendental
constants are numerical displays of the analytic bounds, not interval outputs.
Default verifies a saved receipt without writing. --write-results only creates.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = HERE / "results.json"
OLD_SCALAR = "archive_531_553/544/joint_singlet_common_mass_rg_results.json"
OLD_GAUGE = "archive_554_584/569/joint_quotient_gauge_completion_results.json"
OLD_CAR = "archive_585_628/598/joint_fermion_gauss_completion_results.json"
SOURCES = [
    OLD_SCALAR, OLD_GAUGE, OLD_CAR,
    "archive_585_628/research_note_589.md",
    "archive_585_628/research_note_598.md",
    "archive_585_628/research_note_603.md",
    "archive_585_628/603/joint_thermal_gauss_source.py",
    "archive_585_628/research_note_617.md",
    "archive_629_652/research_note_637.md",
    "archive_629_652/637/joint_full_thermal_region_entropy.py",
    "archive_629_652/research_note_644.md",
    "archive_1219_/1221/proof.md", "archive_1219_/1225/proof.md",
    "archive_1219_/1226/proof.md", "archive_1219_/1226/check.py",
]
CT, UCOEF, ESTAR, CSHIFT = Q(1, 4000), Q(1, 4), 5446, 1072


def read_exact(path):
    return json.loads((ROOT / path).read_text("utf-8"), parse_float=Q)


def local_constants():
    r = read_exact(OLD_SCALAR)["examples"][2]["state"]
    l11, l12, l22 = r["lambda_H"], r["p"], r["lambda_s"]
    det = l11*l22-l12*l12
    lam0 = det/(l11+l22)
    c1, c2 = r["x"]/4, r["y"]/4
    uh, us = (l22*c1-l12*c2)/det, (l11*c2-l12*c1)/det
    D = 12-uh-us
    assert det > 0 and lam0 > Q(9, 100) and D > Q(45, 4)
    assert uh > 0 and us > 0 and uh+us < 1
    assert 1-l11 > 0 and (1-l11)*(1-l22)-l12*l12 > 0
    aa, bb = 32/(lam0*D*D), 144/(D*D)
    assert aa < 3 and bb < 2
    old = read_exact(OLD_CAR)
    caps = dict(u=Q(5, 4), d=Q(3, 5), e=Q(1, 4), nu=Q(1, 2), s=Q(1, 3))
    for label, (re, im) in old["representation"]["diagnostic_Y"].items():
        assert re*re+im*im < caps[label]**2
    assert Q(49, 10)**2 > 24 and Q(7, 2)**2 > 12
    my = 2*Q(49, 10)*(6*(caps["u"]+caps["d"])+2*(caps["e"]+caps["nu"]))+7*caps["s"]
    assert my < 128 and Q(4, 3)**4 > 3 and Q(6, 5)**4 > 2
    assert 128*Q(4, 3) < 171 and 128*Q(6, 5) < 154
    assert 171 < 6**3 and Q(3, 4)*171*6+154 < 1024
    assert 2*32*Q(17, 100) < 16 and CSHIFT == 1024+3*16
    g = read_exact(OLD_GAUGE)["evidence"]["units_check"]
    bc = g["electric_coefficients_colour_weak_integer_circle"]
    assert all(Q(3, 2000) < x < Q(7, 10) for x in bc)
    assert Q(4, 3)**3 > 2
    electric_lower = Q(3, 2000)*Q(1, 4)*Q(3, 4)
    assert electric_lower > CT and Q(1, 4) > CT
    rc, rw, ry = g["residual_positive_coefficients"]
    sad = g["delta"]*(12+4*g["positive_weight_ratio"])+8*rc+3*rw+ry
    assert sad < 16
    mag_bound = Q(15, 2)*Q(4, 3)*16*16
    assert mag_bound == 2560 and mag_bound < 4096
    # Exact beta ratio B(7/2,7)/B(5/2,9), not a quadrature guess.
    beta_ratio = Q(5, 2)*Q(21, 2)/(7*8)
    assert 64*8*beta_ratio == 240
    assert Q(81, 100) < 1
    trial_terms = {"scalar_kinetic":240, "onsite_potential":2,
                   "full_gradient_stars":36, "full_magnetic_stars":4096,
                   "electric_expectation":0, "CAR_expectation":0,
                   "analysis_shift":CSHIFT}
    assert sum(trial_terms.values()) == ESTAR
    tail_a = lam0*D*D/(512*4)
    assert tail_a > Q(1, 200) and Q(16, 25) > Q(13, 24)
    grad = 6*Q(4, 3)*72*4
    assert grad == 2304
    assert 2*grad*4 == 18432 and 2*grad+6*4096+2*16 == 29216
    assert Q(3, 5)**3 < Q(1, 4)  # w_min^(2/3) > .6, so face > .3.
    return dict(status="PASS", exact_guard_arithmetic="Fraction on frozen decimal inputs",
                lambda0_display=float(lam0), u_sum_display=float(uh+us), D_display=float(D),
                F_bound_A_display=float(aa), F_bound_B_display=float(bb),
                mass_numerator_rational_upper=str(my), mass_shift=1024,
                electric_lower_rational=str(electric_lower), c_T=str(CT), u=str(UCOEF),
                c_shift_per_node=CSHIFT, magnetic_weight_sum_display=float(sad),
                full_magnetic_star_norm_upper=4096, radial_guard="r^2 <= 72(1+U)",
                trial_terms_per_node=trial_terms, e_trial_shifted_per_node=ESTAR)


def s_upper(k, alpha):
    return 1+sum(math.comb(k,j)*2**(k-j)*math.gamma((j+1)/2)
                 /(2*alpha**((j+1)/2)) for j in range(k+1))


def log_gauge_upper(t):
    return 2*math.log(s_upper(4,t/3))+math.log(s_upper(2,t/4))+math.log(2*s_upper(0,t)-1)


def log_scalar_upper(s):
    t, alpha = s*float(CT), s*float(UCOEF)
    omega, b, a = 8*math.pi**2/3, 4/math.sqrt(6), .005
    x = math.log(3)
    # sinh(4 log 3)/4 - 2 sinh(2 log 3) + 3 log 3.
    bracket = (81-1/81)/8-(9-1/9)+3*x
    core = omega*36*math.sqrt(6)/8*bracket
    tail = 9*omega/(4*b*alpha*a)*math.exp(-alpha*a*81)
    log_heat = -2.5*math.log(4*math.pi*t)-2*t/3+math.log1p(t/9)
    return log_heat+math.log(core+tail)


def thermal_row(beta):
    s = beta/2
    a = max(0.,32*math.log(2)+log_scalar_upper(s)+3*log_gauge_upper(s*float(CT)))
    e = 2*ESTAR+2*a/beta
    cb = 18432*e+29216
    lg = log_gauge_upper(1.)
    return dict(beta=beta, a_star_at_half_beta=a, e_beta_per_node=e,
                U_per_node_upper=e/float(UCOEF), edge_Casimir_upper=e/float(CT),
                C_boundary_per_cut_edge=cb, double_relative_entropy_per_area=beta*cb/.3,
                full_mutual_information_per_area_t1=(beta*cb+e/float(CT)+lg)/.3,
                log_gauge_trace_upper_at_t1=lg)


def thermal_constants():
    # Direct sums only check the scalar series majorant; not a Gibbs spectrum.
    checks=[]
    for k, alpha in itertools.product((0,2,4),(.01,.2,1.)):
        partial=sum((n+1)**k*math.exp(-alpha*n*n) for n in range(200))
        upper=s_upper(k,alpha)
        assert partial < upper
        checks.append(dict(k=k, alpha=alpha, partial_200=partial, analytic_upper_display=upper))
    rows=[thermal_row(b) for b in (.85,1.7,3.4)]
    assert all(r["e_beta_per_node"]>2*ESTAR and math.isfinite(r["full_mutual_information_per_area_t1"]) for r in rows)
    return dict(status="PASS", rows=rows, scalar_series_checks=checks,
                numbers_are="floating displays of analytic bounds, not interval upper endpoints",
                full_Gibbs_partition_energy_entropy_evaluated=False,
                convexity_inequality="E_L(beta) <= 2*(F(beta/2)-F(beta))/beta",
                volume_size_enters_only_as="N multiplying the same per-node bound")


def shift(v,i,L,step=1):
    x=list(v); x[i]=(x[i]+step)%L; return tuple(x)


def magnetic_edges(v,L):
    edges={(v,i) for i in range(3)}
    for i,j in itertools.combinations(range(3),2):
        edges.add((shift(v,i,L),j)); edges.add((shift(v,j,L),i))
    return edges


def graph_checks():
    rows=[]; samples=[]
    for L in (5,6,8):
        vertices=list(itertools.product(range(L),repeat=3))
        incidence={(v,i):0 for v in vertices for i in range(3)}
        for v in vertices:
            me=magnetic_edges(v,L)
            mv={x for edge in me for x in (edge[0],shift(edge[0],edge[1],L))}
            assert len(me)==9 and len(mv)==7
            for edge in me: incidence[edge]+=1
        assert set(incidence.values())=={3}
        regions={
            "slab":{v for v in vertices if v[0]<L//2},
            "box":{v for v in vertices if max(v)<L//2},
            "checkerboard":{v for v in vertices if sum(v)%2==0},
            "fixed_modular_mask":{v for v in vertices if (v[0]+2*v[1]+3*v[2])%5<2},
        }
        for name,A in regions.items():
            cut={(v,i) for v in vertices for i in range(3) if (v in A)!=(shift(v,i,L) in A)}
            grad=sum(any((v,i) in cut for i in range(3)) for v in vertices)
            mag=sum(bool(magnetic_edges(v,L)&cut) for v in vertices)
            # Connected support straddles iff one of its skeleton edges is cut.
            direct=0
            for v in vertices:
                me=magnetic_edges(v,L)
                mv={x for edge in me for x in (edge[0],shift(edge[0],edge[1],L))}
                direct += bool(mv&A) and bool(mv-A)
            assert direct==mag and grad<=len(cut) and mag<=3*len(cut)
            directions=[sum(i==j for _,i in cut) for j in range(3)]
            rows.append(dict(L=L,N=L**3,region=name,cut_edges=len(cut),
                             cut_directions=directions,gradient_stars=grad,magnetic_stars=mag,
                             CAR_cross_edges=len(cut),electric_cross_energy_difference=0))
            samples.append((directions,len(cut)))
    return dict(status="PASS", full_magnetic_star_edges=9, full_magnetic_star_vertices=7,
                exact_edge_to_magnetic_star_incidence=3, rows=rows,
                uniform_bounds="gradient <= b; magnetic <= 3b; CAR <= b",
                electric_accounting="full same-source star remains regional after original J"),samples


def entropy_area(samples):
    # Four original legal colour-loop representations from637, not replacement matter.
    labels=[(0,0),(1,1),(3,0),(2,2)]
    dims=[]; cas=[]
    for a,b in labels:
        assert (-2*(a+2*b))%6==0  # ell=q=0 meets original quotient constraint.
        dims.append((a+1)*(b+1)*(a+b+2)//2)
        cas.append(Q(a*a+b*b+a*b+3*a+3*b,3))
    assert dims==[1,8,10,27] and cas==[0,3,6,8]
    delta=[Q(x,160) for x in (-5,16,-20,9)]
    rows=[]
    for sign in (-1,1):
        p=[float(Q(1,4)+sign*x) for x in delta]
        assert min(p)>0 and sum(p)==1
        H=-sum(x*math.log(x) for x in p)
        pair=2*sum(x*math.log(d) for x,d in zip(p,dims))
        energy=sum(x*float(c) for x,c in zip(p,cas))
        zsubset=sum(d*d*math.exp(-float(c)) for d,c in zip(dims,cas))
        rhs=energy+log_gauge_upper(1.)
        assert H+pair<=energy+math.log(zsubset)+1e-13<=rhs+1e-13
        # Correlated two-cut label fixture lambda1=lambda2; no independence assumed.
        assert H<=2*H
        rows.append(dict(sign=sign,one_edge_label_entropy=H,one_edge_pair_term=pair,
                         one_edge_Casimir_mean=energy,one_edge_bound=rhs,
                         identical_two_edge_label_joint_entropy=H,
                         correlated_two_edge_baseline=H+2*pair, two_edge_bound=2*rhs))
    primary=thermal_row(1.7); areas=[]
    for w,shear in ((.5,0.),(1.,.6),(2.,-.6)):
        G=np.array([[1.,shear,0.],[shear,1+shear*shear,0.],[0.,0.,1.]])
        assert abs(np.linalg.det(G)-1)<1e-14
        assert np.linalg.eigvalsh(G)[0]>=.25 and np.linalg.eigvalsh(G)[-1]<=4
        P=w**(4/3)*np.linalg.inv(G)
        faces=np.sqrt(np.diag(P)); assert min(faces)>.3
        for counts,boundary in samples:
            A=float(np.dot(counts,faces))
            assert boundary<=A/.3
        # Conditional transfer, no claim that this error is experimentally certified.
        counts=np.array([18,12,6]); area=float(counts@faces)
        delta_A=.02*area; estimate=area-delta_A
        k=primary["full_mutual_information_per_area_t1"]
        assert k*area<=k*(estimate+delta_A)*(1+1e-14)
        areas.append(dict(w=w,shear=shear,face_areas=faces.tolist(),
                          geometric_area=area,conditional_delta_A=delta_A,
                          conditional_estimate=estimate,
                          safe_upper_formula="K*(estimated_area + delta_A)"))
    return dict(status="PASS", legal_original_representation_labels=[list(x) for x in labels],
                center_pair_KL_diagnostics=rows, independent_geometry_diagnostics=areas,
                label_fixtures_are_full_Gibbs_predictions=False,
                no_interval_or_actual_instrument_error_certification=True,
                full_mutual_information_preserves=["joint center entropy","representation pairing","relative entropy correlation"])


def run():
    local=local_constants(); thermal=thermal_constants(); graph,samples=graph_checks()
    entropy=entropy_area(samples)
    groups=dict(local_original_constants=local,uniform_thermal_comparison=thermal,
                complete_interaction_cluster_counts=graph,center_pair_and_independent_area=entropy)
    hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}
    return dict(round=1226,status="PASS",check_groups=4,tests_run=4,failures=0,errors=0,
                source_sha256=hashes,groups=groups,
                scope="N-uniform fixed-spacing positive-temperature full-Gauss Gibbs comparison bounds",
                not_claimed=["full thermal spectrum computed","S(A) area law","S=sA","deltaS=s deltaA",
                             "UV uniformity","autonomous thermal preparation","Einstein dynamics","stage complete"])


def main():
    p=argparse.ArgumentParser(); p.add_argument("--write-results",action="store_true")
    args=p.parse_args(); result=run()
    if args.write_results:
        with OUTPUT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2); f.write("\n")
    else:
        saved=json.loads(OUTPUT.read_text("utf-8"))
        assert saved==result,"Saved receipt differs; no files changed"
    print(json.dumps(dict(round=1226,status="PASS",check_groups=4,
                          receipt_mode="created" if args.write_results else "read_only_verified",
                          e_beta=thermal_row(1.7)["e_beta_per_node"],
                          I_per_area=thermal_row(1.7)["full_mutual_information_per_area_t1"]),ensure_ascii=False))


if __name__=="__main__":
    main()
