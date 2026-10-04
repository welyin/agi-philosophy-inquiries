"""Executed740 entry: original total neutral response, not the pure spectral I.

One-loop resummation is tested as a mathematical all-frequency equation.
Positive-Laplace poles do not prove the microscopic theory is unstable.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_covariant_response_closure as full
import joint_mixed_neutral_response as mixed
TARGET=HERE/'total_pole_entry_results.json'


def constants():
    C,B,_=mixed.decomposition(1.2+.3j,256)
    return C,B


def exact_remainder(logs,n=192):
    """Exact integral up to u=32; omitted matrix norm <= sum||Cj|| exp(-64).

    u=-log(x)=t^2 regularizes the threshold square roots. No exp(logs)
    is evaluated, so the second pole can be represented by its logarithm.
    """
    z,w=np.polynomial.legendre.leggauss(n)
    t=(z+1)*np.sqrt(32)/2;q=w*np.sqrt(32)/2*2*t;u=t*t
    ans=np.zeros((2,2))
    for _,a,r,kind,Cj in mixed.terms():
        gap=1-mixed.profile(np.exp(-u),r,kind)
        suppression=np.exp(-np.logaddexp(0.,2*(logs-np.log(a)-u)))
        ans+=Cj*(.5*np.log1p(np.exp(2*(np.log(a)-logs)))+np.dot(q,gap*suppression))
    return ans


def remainder_log_bound(logs):
    rows=[]
    for _,a,r,kind,Cj in mixed.terms():
        norm=float(np.linalg.norm(Cj,2))
        if norm==0:continue
        rows.append(np.log(norm)+2*(np.log(a)-logs)+
                    np.log(.5+np.logaddexp(0.,2*(logs-np.log(a)))))
    return float(np.logaddexp.reduce(rows))


def reduced_schur(logs,n=192):
    x,G,V=full.local_data();C,B=constants();I=C*logs+B+exact_remainder(logs,n)
    e=np.exp(-2*logs);a=x@G@x-6+e*(x@V@x)
    b=-G@x-e*V@x
    assert a<0  # this probe only uses logs >= 5, beyond its auxiliary zero
    return G-I+e*V-np.outer(b,b)/a


def bisect_eigen(which,lo,hi):
    left=float(np.linalg.eigvalsh(reduced_schur(lo))[which])
    right=float(np.linalg.eigvalsh(reduced_schur(hi))[which])
    assert left>0 and right<0
    initial=[lo,hi]
    for _ in range(52):
        middle=(lo+hi)/2
        value=float(np.linalg.eigvalsh(reduced_schur(middle))[which])
        if value>0:lo=middle
        else:hi=middle
    root=(lo+hi)/2
    eigen,vectors=np.linalg.eigh(reduced_schur(root))
    direction=vectors[:,which]
    eps=1e-4
    slope=float(direction@(reduced_schur(root+eps)-reduced_schur(root-eps))@direction/(2*eps))
    return dict(eigen_index=which,numerical_initial_log_bracket=initial,
                log_positive_Laplace_pole=root,log10_positive_Laplace_pole=root/np.log(10),
                pole_value_if_representable=float(np.exp(root)) if root<700 else None,
                scalar_Schur_eigenvalues=eigen.tolist(),crossing_slope=slope,
                exact_remainder_norm_log_upper=remainder_log_bound(root),
                quadrature_refinement_error=float(np.max(abs(reduced_schur(root)-reduced_schur(root,384)))),
                normalized_active_direction=direction.tolist())


def run():
    C,B=constants();x,G,V=full.local_data()
    # Current q0=0.27 fixture is not used: all data are632's matched vacuum.
    assert min(np.linalg.eigvalsh(V))>0 and min(np.linalg.eigvalsh(C))>0
    A=6-x@G@x
    assert np.max(abs(G+np.outer(G@x,G@x)/A-np.eye(2)))<1e-14
    roots=[bisect_eigen(0,5.,12.),bisect_eigen(1,500.,1000.)]
    for row in roots:
        assert row['crossing_slope']<0
        assert row['quadrature_refinement_error']<2e-13
    L=roots[0]['log_positive_Laplace_pole']
    # At representable first pole check the ORIGINAL delta-x/sigma total block.
    # This avoids mistaking a Schur denominator zero for a true pole.
    I=C*L+B+exact_remainder(L);e=np.exp(-2*L)
    M=np.block([[G+e*V-I,(-I@x)[:,None]],[(-x@I)[None,:],np.array([[-6-x@I@x]])]])
    ev,vec=np.linalg.eigh(M);idx=int(np.argmin(abs(ev)))
    residual=float(np.linalg.norm(M@vec[:,idx]))
    assert residual<2e-14
    pure_min=float(min(np.linalg.eigvalsh(I)))
    assert pure_min>0
    names=('research_note_601.md','research_note_630.md','research_note_632.md',
           'research_note_738.md','research_note_739.md','joint_mixed_neutral_response.py',
           'joint_covariant_response_closure.py')
    return dict(entry_round=740,latest_completed_round=739,formal_count_unchanged=3417,
                roots=roots,full_original_scaled_block_eigenvalues=ev.tolist(),
                full_block_null_residual=residual,pure_I_minimum_at_total_pole=pure_min,
                baseline_local_scheme=dict(f=1,ell=[0,0],Z='original G',
                                           alpha=0,beta=0,potential='original632 matched one-loop Hessian'),
                omitted_u_tail_norm_upper=float(sum(np.linalg.norm(t[-1],2) for t in mixed.terms())*np.exp(-64)),
                numerical_brackets_not_interval_arithmetic_certificates=True,
                analytic_existence_is_by_inertia_and_positive_log_coefficient=True,
                no_physical_UV_validity_or_microscopic_instability_claim=True,
                dependencies={name:hashlib.sha256((ARCHIVE/name).read_bytes()).hexdigest() for name in names},
                all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
