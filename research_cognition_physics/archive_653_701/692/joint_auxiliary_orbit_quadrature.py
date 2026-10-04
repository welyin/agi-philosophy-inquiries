"""692: full-functional auxiliary orbit reduction and positive radial cubature.

The exact theorem keeps all bosonic/Gauss averages. Fixed-background radial
checks are not a sphere-integral replacement or a positivity certificate.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_gauss_boundary_functional as base

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_auxiliary_orbit_quadrature_results.json'
spec=importlib.util.spec_from_file_location('entry692',HERE/'round692_drafts/odd_gauss_source_probe.py')
entry=importlib.util.module_from_spec(spec);spec.loader.exec_module(entry)


def moment(k):return F(12,(k+3)*(k+4))


def exact_solve(a,b):
    a=[list(row)+[v] for row,v in zip(a,b)]
    n=len(b)
    for j in range(n):
        pivot=next(i for i in range(j,n) if a[i][j])
        a[j],a[pivot]=a[pivot],a[j]
        d=a[j][j];a[j]=[v/d for v in a[j]]
        for i in range(n):
            if i!=j:
                d=a[i][j];a[i]=[v-d*w for v,w in zip(a[i],a[j])]
    return [row[-1] for row in a]


def rule(n):
    # Jacobi(alpha=1,beta=2), x=2r-1; normalized Beta(3,2) measure.
    a,b=1,2
    diagonal=[];off=[]
    for k in range(n):
        s=2*k+a+b
        diagonal.append((1+(b*b-a*a)/(s*(s+2)))/2)
        if k:
            off.append(np.sqrt(k*(k+a)*(k+b)*(k+a+b)/(s*s*(s*s-1))))
    jacobi=np.diag(diagonal)+np.diag(off,1)+np.diag(off,-1)
    points,vec=np.linalg.eigh(jacobi)
    weights=vec[0]**2
    assert np.min(points)>0 and np.max(points)<1 and np.min(weights)>0
    return points,weights


def rule_check():
    coefficients=exact_solve([[moment(i+j) for j in range(9)] for i in range(9)],
                             [-moment(i+9) for i in range(9)])+[F(1)]
    for j in range(9):assert sum(c*moment(i+j) for i,c in enumerate(coefficients))==0
    r,w=rule(9)
    residual=max(abs(np.sum(w*r**k)-float(moment(k))) for k in range(18))
    assert residual<2e-14
    p=np.array([float(c) for c in coefficients])
    root_residual=max(abs(np.polynomial.polynomial.polyval(r,p)))
    assert root_residual<1e-12
    # Independent integration of invariant spherical monomials. Exact beta
    # moments are also checked against the 6+4 Cartesian sphere moments.
    mixed=0.
    import math
    for a in range(17):
        for b in range(17-a):
            exact=F(math.factorial(a+2)*math.factorial(b+1)*12,math.factorial(a+b+4))
            mixed=max(mixed,abs(np.sum(w*r**a*(1-r)**b)-float(exact)))
    assert mixed<2e-14
    return dict(nodes=r.tolist(),positive_normalized_weights=w.tolist(),
        exact_monic_degree9_polynomial=[str(c) for c in coefficients],
        exact_orthogonality_tests=9,degree0_through17_moment_error=float(residual),
        mixed_invariant_moment_error=float(mixed),polynomial_root_error=float(root_residual),
        local_E_degree_bound=32,local_invariant_radial_degree_bound=16,
        nodes_per_spacetime_site=9,four_site_auxiliary_assignments=9**4,
        Gaussian_quadrature_rule_exactness_is_analytic=True)


def original_orbits():
    center=base.prior.rep(np.eye(3),np.eye(2),-1.+0j)
    o=base.prior.vector_rotation(center)
    pc=(np.eye(10)+o)/2;pw=np.eye(10)-pc
    assert np.max(abs(pc@pc-pc))<1e-13 and abs(np.trace(pc)-6)<1e-13
    assert abs(np.trace(pw)-4)<1e-13
    ci=int(np.argmax(np.diag(pc)));wi=int(np.argmax(np.diag(pw)))
    c=pc[:,ci]/np.linalg.norm(pc[:,ci]);v=pw[:,wi]/np.linalg.norm(pw[:,wi])
    r=np.array([.16,.37,.61,.82])
    canonical=np.sqrt(r)[:,None]*c+np.sqrt(1-r)[:,None]*v
    links,_,phis=entry.sources.fixture()
    groups=[base.prior.group(69250+i,.41) for i in range(4)]
    lt,et,pt,reps=base.transform(links,canonical,phis,groups)
    orbit_errors=[]
    for gi,e,ri in zip(reps,et,r):
        og=base.prior.vector_rotation(gi)
        orbit_errors.append(float(max(np.max(abs(og@pc-pc@og)),abs(e@pc@e-ri),abs(e@e-1))))
    assert max(orbit_errors)<3e-13
    first=entry.evaluate(links,canonical,phis)
    changed=entry.evaluate(lt,et,pt)
    source_error=float(abs(changed['source']/first['source']-1))
    weight_error=float(abs(changed['weight']/first['weight']-1))
    assert source_error<2e-9 and weight_error<2e-9
    # Vary one radial coordinate in the ACTUAL128-mode auxiliary Pfaffian.
    # Sign-average makes the conditional expression polynomial in r. This
    # check alone does not replace the S9 integral at fixed gauge background.
    u,_,_,_,gap=base.kernel(links)
    coeff=[]
    for x in range(4):
        ux=u[64*x:64*(x+1)]
        coeff.append(np.array([ux.T@np.kron(base.internal.B,t)@ux for t in base.internal.T]))
    coeff=np.array(coeff)
    fixed=np.einsum('xa,xaij->ij',canonical[1:],coeff[1:])
    ac=np.einsum('a,aij->ij',c,coeff[0]);aw=np.einsum('a,aij->ij',v,coeff[0])
    def radial(x):
        return sum(base.pf(fixed+sc*np.sqrt(x)*ac+sw*np.sqrt(1-x)*aw)
                   for sc in (-1,1) for sw in (-1,1))/4
    points,weights=rule(9)
    low=sum(w*radial(x) for x,w in zip(points,weights))
    # Independent Legendre rule integrates polynomial * explicit cubic weight.
    xx,ww=np.polynomial.legendre.leggauss(12)
    high=sum(w*6*((x+1)/2)**2*(1-(x+1)/2)*radial((x+1)/2) for x,w in zip(xx,ww))
    radial_error=float(abs(low-high)/max(abs(low),abs(high),1e-280))
    assert max(abs(low),abs(high))>1e-12 and radial_error<2e-9
    return dict(original_vector_color_projector=pc.real.tolist(),color_rank=6,weak_rank=4,
        canonical_color_axis=ci,canonical_weak_axis=wi,local_orbit_errors=orbit_errors,
        original_scalar_weight_covariance_error=weight_error,
        original_six_field_B_source_covariance_error=source_error,
        original_nonflat_Wilson_gap=gap,auxiliary_spectral_dimension=u.shape[1],
        sign_symmetrized_conditional_radial_integral=base.old.cpair(low),
        independent_Legendre_integral=base.old.cpair(high),radial_relative_error=radial_error,
        full_Hb_Gauss_integral_not_numerically_evaluated=True,
        radial_rule_not_claimed_as_fixed_background_sphere_rule=True)


def run():
    deps=('research_note_614.md','research_note_653.md','research_note_657.md','research_note_673.md',
          'research_note_675.md','research_note_676.md','research_note_679.md','research_note_680.md',
          'research_note_691.md','round692_drafts/odd_gauss_source_probe.py',
          'round692_drafts/odd_gauss_source_probe_results.json','round692_drafts/entry_checks.json')
    return dict(date='2026-10-02',round=692,tests_run=2,failures=0,errors=0,
        exact_radial_rule=rule_check(),original_orbit_and_full_source=original_orbits(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(full_functional_auxiliary_orbit_reduction=True,
            all_E_independent_Gauss_invariant_physical_sources_retained=True,
            original_auxiliary16_measure_Hb_double_Haar_and_normalization_retained=True,
            nine_positive_radial_nodes_exact_after_complete_original_average=True,
            general_E_dependent_extended_test_algebra_covered_by_nine_nodes=False,
            fixed_background_sphere_integral_replaced=False,
            original_dynamic_RP_or_HF_identification=False,
            continuum_or_quantum_GR_completed=False,old_space_and_goal_unchanged=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    out=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    else:assert out==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=692,tests=2,all_checks_passed=True,
        radial_moment_error=out['exact_radial_rule']['degree0_through17_moment_error'],
        original_source_error=out['original_orbit_and_full_source']['original_six_field_B_source_covariance_error'],
        original_radial_error=out['original_orbit_and_full_source']['radial_relative_error'])))
