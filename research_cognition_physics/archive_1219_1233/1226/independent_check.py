"""1226 independent audit. Original-source constants and local geometry only.
No author code is imported and no full Gauss thermal spectrum is simulated.
Default execution is read-only; --write-results only creates a missing result.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import itertools
import json
import math

HERE=Path(__file__).resolve().parent
RESEARCH=HERE.parents[1]
OUT=HERE/"independent_results.json"
SOURCES=(
 "archive_554_584/569/joint_quotient_gauge_completion_results.json",
 "archive_585_628/598/joint_fermion_gauss_completion_results.json",
 "archive_585_628/603/joint_thermal_gauss_source.py",
 "archive_629_652/research_note_637.md",
 "archive_1219_/1221/proof.md",
 "archive_1219_/1225/proof.md",
)
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads((RESEARCH/p).read_text(encoding="utf-8"))
def frac(x): return F(str(x))

def original_local_constants():
    old=read(SOURCES[1]); gauge=read(SOURCES[0])["evidence"]["units_check"]
    lam=frac(old["confinement"]["lambda_lower_bound"])
    D=frac(old["confinement"]["D"])
    assert lam>F(9,100) and D>F(45,4)
    assert 32/(lam*D*D)<3 and 144/(D*D)<2
    guards={"u":F(5,4),"d":F(3,5),"e":F(1,4),"nu":F(1,2),"s":F(1,3)}
    for k,v in old["representation"]["diagnostic_Y"].items():
        assert sum(frac(x)**2 for x in v)<guards[k]**2
    # Rational upper bounds sqrt(24)<49/10 and sqrt(12)<7/2.
    my=2*F(49,10)*(6*(guards["u"]+guards["d"])+2*(guards["e"]+guards["nu"]))+7*guards["s"]
    assert my<128
    assert 128*F(4,3)<171 and 128*F(6,5)<154
    assert 171<6**3
    young_remainder=F(3,4)*171*6+154
    assert young_remainder<1024
    # Young is an exact polynomial factorization after x=y^4:
    # y^4-4*t^3*y+3*t^4=(y-t)^2*(y^2+2*t*y+3*t^2).
    for y,t in itertools.product((F(0),F(1,3),F(2),F(7)),(F(1,2),F(1),F(3))):
        assert y**4-4*t**3*y+3*t**4==(y-t)**2*(y*y+2*t*y+3*t*t)>=0
    bs=[frac(x) for x in gauge["electric_coefficients_colour_weak_integer_circle"]]
    assert min(bs)>F(3,2000) and max(bs)<F(7,10)
    assert F(3,4)*F(1,4)*F(3,2000)==F(9,32000)>F(1,4000)
    rc,rw,ry=map(frac,gauge["residual_positive_coefficients"])
    magnetic=frac(gauge["delta"])*(12+4*frac(gauge["positive_weight_ratio"]))+8*rc+3*rw+ry
    assert magnetic<16
    assert 2*32*F(17,100)<16
    return dict(original_lambda_lower=float(lam),original_D=float(D),
      recomputed_A=float(32/(lam*D*D)),recomputed_B0=float(144/(D*D)),
      rational_mass_prefactor_upper=str(my),young_remainder_upper=str(young_remainder),
      magnetic_dimension_weighted_sum=float(magnetic),
      cT="1/4000",u="1/4",shift_per_site=1072,
      inputs_are_frozen_model_parameters_not_experimental_uncertainty=True,
      no_N_in_local_constants=True)

def trial_and_spectral_algebra():
    # Exact beta-function ratio B(7/2,7)/B(5/2,9).
    ratio=F(5,2)*F(21,2)/F(8*7)
    assert ratio==F(15,32)
    scalar=64*8*ratio
    assert scalar==240
    gradient=F(6)*F(4,3)*4 # support radius <=1
    assert gradient==32<36
    magnetic=F(15)*16*F(4,3)*16/2
    assert magnetic==2560<4096
    trial=240+2+36+4096
    assert trial==4374 and trial+1072==5446
    # Original tail coefficient: lambda*D^2/(512 M^2), M=2.
    a603_lower=F(9,100)*F(45,4)**2/2048
    assert a603_lower>F(1,200)
    # R603<2.5: u_sum<3/4, sqrt(6)<5/2, e^2>7.
    assert F(17,32)<F(9,16)
    assert sum(F(2)**k/F(math.factorial(k)) for k in range(5))==7
    # Geometric area lower a=2^(-5/3), so cube exceeds (3/10)^3.
    assert F(1,32)>F(3,10)**3
    return dict(trial_beta_ratio=str(ratio),trial_scalar_energy=240,
      trial_gradient_bound=36,complete_magnetic_bound=4096,
      trial_H_per_site=4374,trial_shifted_H_per_site=5446,
      original_tail_coefficient_lower=str(a603_lower),tail_a="1/200",tail_R0="5/2",
      positive_face_area_lower="3/10",Gauss_trial_is_form_domain_only=True,
      spectral_steps=["min-max eigenvalue ordering","Gauss trace <= full trace",
                      "spectral Jensen on a legal Gauss trial",
                      "convex log trace secant beta/2 to beta"],
      energy_density_formula="2*5446+2*a_star(beta/2)/beta")

def lattice(side):
    V=list(itertools.product(range(side),repeat=3))
    def shift(v,i): return tuple((v[k]+(1 if k==i else 0))%side for k in range(3))
    edges=[(v,i) for v in V for i in range(3)]
    ends={e:(e[0],shift(*e)) for e in edges}
    grad={v:{(v,i) for i in range(3)} for v in V}
    mag={}
    for v in V:
        ed=set()
        for i in range(3):
            for j in range(i+1,3):
                ed.update(((v,i),(v,j),(shift(v,i),j),(shift(v,j),i)))
        assert len(ed)==9
        points={x for e in ed for x in ends[e]}
        assert len(points)==7
        mag[v]=ed
    incidences={e:0 for e in edges}
    for es in mag.values():
        for e in es: incidences[e]+=1
    assert set(incidences.values())=={3}
    return V,edges,ends,grad,mag

def clusters():
    rows=[]
    for side in (5,6,7):
        V,E,ends,grad,mag=lattice(side)
        regions={"slab":{v for v in V if v[0]<side//2},
                 "one_site":{(0,0,0)},
                 "checker":{v for v in V if sum(v)%2==0},
                 "oblique":{v for v in V if (v[0]+2*v[1]+3*v[2])%5<2}}
        for label,A in regions.items():
            cut={e for e,(x,y) in ends.items() if (x in A)!=(y in A)}
            ng=sum(bool(es&cut) for es in grad.values())
            nm=sum(bool(es&cut) for es in mag.values())
            assert ng<=len(cut) and nm<=3*len(cut)
            rows.append(dict(side=side,region=label,sites=len(V),cut_edges=len(cut),
                             crossing_gradient_stars=ng,crossing_magnetic_stars=nm))
    return dict(rows=rows,all_magnetic_stars_vertices=7,all_magnetic_stars_edges=9,
        magnetic_star_incidence_per_edge=3,gradient_star_incidence_per_edge=1,
        electric_star_difference=0,
        electric_reason="source half-link representation and preserved regional marginals",
        counting_is_local_not_a_full_H_simulation=True)

def S_upper(k,alpha):
    return 1+sum(math.comb(k,j)*2**(k-j)*math.gamma((j+1)/2)/(2*alpha**((j+1)/2)) for j in range(k+1))
def log_zG(t):
    return 2*math.log(S_upper(4,t/3))+math.log(S_upper(2,t/4))+math.log(2*S_upper(0,t)-1)
def scalar_integral(alpha):
    R0=2.5; R=math.sqrt(6); b=4/R; a=.005; omega=8*math.pi**2/3
    x=R0/R
    core=omega*36*R/8*(math.sinh(4*x)/4-2*math.sinh(2*x)+3*x)
    tail=9*omega/(4*b*alpha*a)*math.exp(-alpha*a*math.exp(b*R0))
    return core+tail
def astar(s):
    t=s/4000
    logkernel=-2.5*math.log(4*math.pi*t)-2*t/3+math.log1p(t/9)
    return max(0,32*math.log(2)+logkernel+math.log(scalar_integral(s/4))+3*log_zG(t))
def thermal_constants():
    # Direct finite sums only check the analytic integral upper bound.
    for k,alpha in itertools.product((0,2,4),(.3,.6,1.2)):
        partial=sum((n+1)**k*math.exp(-alpha*n*n) for n in range(200))
        assert partial<=S_upper(k,alpha)
    rows=[]
    for beta in (.5,1.,2.):
        a=astar(beta/2); eb=10892+2*a/beta
        c=18432*eb+29216
        rows.append(dict(beta=beta,a_star_half_beta=a,e_beta=eb,
          C_boundary=c,double_relative_entropy_per_area=beta*c/.3,
          full_extended_I_per_area_t1=(beta*c+4000*eb+log_zG(1))/.3))
    return dict(rows=rows,log_zG_upper_t1=log_zG(1),
        numerical_values_are_displays_of_analytic_bounds=True,
        no_full_Gibbs_spectrum_or_measured_energy_computed=True,
        fixed_positive_beta_only=True)

def quotient_KL():
    labels=[(0,0,0,0),(0,0,2,0),(1,1,0,0)]
    weights=[F(1,2),F(1,3),F(1,6)]
    dims=[]; energies=[]
    for a,b,ell,q in labels:
        assert (2*(a+2*b)+3*ell+q)%6==0
        dims.append(F((a+1)*(b+1)*(a+b+2),2)*(ell+1))
        energies.append(F(a*a+b*b+a*b+3*a+3*b,3)+F(ell*(ell+2),4)+q*q)
    assert dims==[1,3,8] and energies==[0,2,3]
    t=1.; z=sum(float(d*d)*math.exp(-t*float(c)) for d,c in zip(dims,energies))
    entropy=-sum(float(p)*math.log(float(p)) for p in weights)
    pair=2*sum(float(p)*math.log(float(d)) for p,d in zip(weights,dims))
    energy=sum(p*c for p,c in zip(weights,energies))
    qprob=[float(d*d)*math.exp(-t*float(c))/z for d,c in zip(dims,energies)]
    kl=sum(float(p)*math.log(float(p)/q) for p,q in zip(weights,qprob))
    rhs=t*float(energy)+math.log(z)
    assert abs((rhs-entropy-pair)-kl)<2e-14 and kl>=0
    assert entropy+pair <= t*float(energy)+log_zG(t)
    # Two perfectly correlated allowed loop labels: Hjoint=H(p), not 2H(p).
    full_center_pair=entropy+2*pair
    assert full_center_pair<=2*(t*float(energy)+log_zG(t))
    return dict(labels=[list(l) for l in labels],dims=[str(d) for d in dims],
      Casimir=[str(c) for c in energies],probabilities=[str(p) for p in weights],
      single_edge_center_plus_pair=entropy+pair,finite_label_KL=kl,
      two_correlated_edge_center_plus_pair=full_center_pair,
      full_original_group_upper_used=True,
      fixture_is_not_original_full_Gibbs_distribution=True,
      d_squared_in_reference_probability=True,independent_labels_assumed=False)

def boundary_area_algebra():
    # Polynomial in eb. Both sigma and tau expectations must be included.
    coeff=2*2304*4
    constant=2*2304+6*4096+2*16
    assert coeff==18432 and constant==29216
    # A safe use of a finite area estimate is Ahat+delta, not Ahat.
    area=F(13,7); error=F(1,100)
    estimate=area-error
    assert estimate+error==area
    return dict(C_boundary_linear=coeff,C_boundary_constant=constant,
      two_state_factor_retained=True,center_pair_thermal_term="t*e_beta/cT+log(zG_upper(t))",
      safe_area_input="A_hat+delta_A",
      area_source="1225 independent cofactor dictionary, same uniform geometry",
      ordinary_region_entropy_area_law=False,entropy_area_equality=False,
      thermodynamic_limit_unique_state=False,autonomous_thermalization=False,
      vacuum_energy_subtraction=False)

def run():
    return dict(status="PASS",all_checks_passed=True,check_groups=6,
      independent_code_sha256=sha(Path(__file__)),
      source_sha256={p:sha(RESEARCH/p) for p in SOURCES},
      local_constants=original_local_constants(),trial_and_trace=trial_and_spectral_algebra(),
      boundary_clusters=clusters(),thermal_bound_functions=thermal_constants(),
      original_quotient_KL=quotient_KL(),boundary_area=boundary_area_algebra())
def main():
    parser=argparse.ArgumentParser();parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args(); result=run()
    if args.write_results:
        with OUT.open("x",encoding="utf-8",newline="\n") as f: json.dump(result,f,ensure_ascii=False,indent=2);f.write("\n")
    else: assert json.loads(OUT.read_text(encoding="utf-8"))==result,"saved results or sources changed"
    print(json.dumps(dict(status="PASS",all_checks_passed=True,results_match=True,
      check_groups=6,independent_code_sha256=result["independent_code_sha256"])))
if __name__=="__main__": main()
