"""1038 limited verifier: independent 4-vector boosts and exact certificates.

Different mathematical constructions, not a separate author's scientific review.
Writes only its own initial receipt with --write; otherwise checks existing data.
"""
from pathlib import Path
from fractions import Fraction as F
from urllib.parse import unquote
import argparse
import ast
import hashlib
import json
import math
import re
import numpy as np
import moving_direction_transport as science

HERE=Path(__file__).resolve().parent
BASE=HERE.parents[1]
NOTE=HERE.parent/'research_note_1038.md'
OUT=HERE/'research_round_1038_checks.json'
OWN=['moving_direction_transport.py','moving_direction_transport_results.json',
     'drafts/moving_direction_transport_derivation.md','selection_audit.md',
     'input_dependency_update.md','NEXT.md','review.md','verify_round1038.py']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def boost4(p):
    p=np.asarray(p,float);e=np.sqrt(1+p@p)
    b=np.eye(4);b[0,0]=e;b[0,1:]=p;b[1:,0]=p
    b[1:,1:]+=np.outer(p,p)/(e+1)
    return b


def axis_boost4(n,xi):
    n=np.asarray(n,float);c,s=np.cosh(xi),np.sinh(xi)
    a=np.eye(4);a[0,0]=c;a[0,1:]=s*n;a[1:,0]=s*n
    a[1:,1:]+=(c-1)*np.outer(n,n)
    return a


def four_vector_checks():
    rng=np.random.default_rng(41038)
    pauli=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
    eta=np.diag([-1.,1.,1.,1.])
    max_lorentz=max_rotation=max_time=max_dv=0.
    for _ in range(28):
        p=rng.normal(size=3)*1.2
        n=rng.normal(size=3);n/=np.linalg.norm(n)
        xi=rng.uniform(-1.2,1.2)
        a=axis_boost4(n,xi)
        ep=a@np.r_[np.sqrt(1+p@p),p]
        w4=np.linalg.solve(boost4(ep[1:]),a@boost4(p))
        w2=science.wigner_closed(p,n,xi)
        r2=np.array([[np.trace(si@w2@sj@w2.conj().T).real/2 for sj in pauli] for si in pauli])
        max_lorentz=max(max_lorentz,float(np.max(abs(a.T@eta@a-eta))))
        max_rotation=max(max_rotation,float(np.max(abs(w4[1:,1:]-r2))))
        max_time=max(max_time,float(np.max(abs(w4[0]-[1,0,0,0]))),float(np.max(abs(w4[:,0]-[1,0,0,0]))))
        e=np.sqrt(1+p@p)
        dv=np.eye(3)/(e+1)-np.outer(p,p)/(e*(e+1)**2)
        max_dv=max(max_dv,float(np.linalg.eigvalsh(dv).max()))
        assert np.linalg.eigvalsh(dv).min()>0 and np.linalg.eigvalsh(dv).max()<=.5+1e-13
    assert max(max_lorentz,max_rotation,max_time)<3e-13
    return dict(cases=28,lorentz_metric_residual=max_lorentz,
                induced_rotation_residual=max_rotation,little_group_time_residual=max_time,
                sampled_velocity_derivative_max=max_dv,analytic_derivative_bound=.5)


def exact_witness():
    # Rational Bloch rotation, obtained independently from a=4/sqrt17, b=-y/sqrt17.
    co=F(15,17);si=F(8,17)
    assert co*co+si*si==1
    pa=(1+si)/2;pb=(1-si)/2
    lower=pa-pb-F(1,100)
    assert (pa,pb,lower)==(F(25,34),F(9,34),F(783,1700))
    rho_a=np.diag([.5,0,0,.5]);rho_b=np.diag([0,.5,.5,0])
    def trace_spin(r):return np.trace(r.reshape(2,2,2,2),axis1=1,axis2=3)
    def trace_p(r):return np.trace(r.reshape(2,2,2,2),axis1=0,axis2=2)
    assert np.array_equal(trace_spin(rho_a),trace_spin(rho_b))
    assert np.array_equal(trace_p(rho_a),trace_p(rho_b))
    return dict(central_probabilities=[str(pa),str(pb)],normal_packet_gap_lower=str(lower),
                any_marginals_only_error_lower=str(lower/2),both_marginals_equal=True)


def independent_cube_quadrature():
    # Cartesian tensor quadrature, no spinor matrices and no spherical-rule reuse.
    grid,weight=np.polynomial.legendre.leggauss(28)
    x,y,z=np.meshgrid(grid,grid,grid,indexing='ij')
    rr=x*x+y*y+z*z;sel=rr<1
    pp=np.stack([x[sel],y[sel],z[sel]],axis=-1)*.01+np.array([0,0,4/3])
    ee=np.sqrt(1+(pp*pp).sum(axis=1))
    v=pp/(ee[:,None]+1)
    wt=np.einsum('i,j,k->ijk',weight,weight,weight)[sel]*np.exp(-2/(1-rr[sel]))/(2*ee)
    wt/=wt.sum()
    def x_component(vv):
        a=1+.5*vv[:,0]
        return a*vv[:,2]/(1+vv[:,0]+.25*np.sum(vv*vv,axis=1))
    gap=float(wt@((x_component(v)-x_component(-v))/2))
    saved=json.loads(science.OUT.read_text('utf8'))['normal_packet']['fine']['probability_gap']
    assert abs(gap-saved)<2e-8 and gap>float(F(783,1700))
    return dict(points=len(wt),probability_gap=gap,difference_from_science=abs(gap-saved),
                quadrature_is_calibration_not_error_proof=True)


def checks():
    saved=json.loads(science.OUT.read_text('utf8'))
    science.compare(science.run(),saved)
    for p,h in saved['historical_source_sha256'].items():assert sha(BASE/p)==h,p
    assert len(saved['historical_source_sha256'])==15
    assert saved['new_calibration_groups']==1 and saved['new_adopted_cognitive_axioms']==0
    assert not saved['goal_completed']
    assert not saved['algebra_closure']['single_shot_linear_span_classified']
    assert saved['algebra_closure']['algebra_closure_is_an_extra_task_contract']
    for path in HERE.glob('*.py'):ast.parse(path.read_text('utf8'))
    files=[HERE/p for p in OWN]+[NOTE]
    missing=[];links=0
    for p in files:
        assert p.is_file(),p
        if p.suffix!='.md':continue
        text=p.read_text('utf8')
        assert '2026-10-09' not in text,p
        text=re.sub(r'\$\$.*?\$\$|^```.*?^```\s*$',lambda m:' '*len(m[0]),text,flags=re.S|re.M)
        for raw in re.findall(r'(?<!!)\[[^\]\n]*\]\(([^)\n]+)\)',text):
            target=unquote(raw.strip().split('#',1)[0])
            if not target or re.match(r'^[a-zA-Z]+:',target):continue
            links+=1
            target_path=(p.parent/target).resolve()
            if target_path==OUT:continue
            if not target_path.is_file():missing.append([str(p),target])
    assert not missing,missing
    independent=dict(four_vectors=four_vector_checks(),rational_witness=exact_witness(),
                     cartesian_packet=independent_cube_quadrature())
    return dict(round=1038,date='2026-10-08',status='passed',new_calibration_groups=1,
                global_cumulative_count='assigned_by_root',new_adopted_cognitive_axioms=0,
                independent_algorithms=independent,
                historical_source_sha256=saved['historical_source_sha256'],
                frozen_file_sha256={str(p.relative_to(BASE)).replace('\\','/'):sha(p) for p in files},
                frozen_files=len(files),local_links_checked=links,
                same_author_verifier=True,independent_agent_review='see review.md; root signs separately',
                no_physical_detector_or_full_parent_certification=True,goal_completed=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=checks()
    if args.write:
        with OUT.open('x',encoding='utf8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2,allow_nan=False);stream.write('\n')
    elif OUT.exists():science.compare(result,json.loads(OUT.read_text('utf8')))
    print(json.dumps(dict(round=1038,status=result['status'],frozen_files=result['frozen_files'],
                          local_links=result['local_links_checked'],
                          independent_algorithms=result['independent_algorithms']),ensure_ascii=False))
