"""678 entry: sparse original-Wilson recurrence, boundary map and bulk factor.
Linear-system identity only. Not a complete Grassmann measure or physical time.
"""
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_rational_physical_limit as current
base=current.base
TARGET=HERE/'local_block_probe_results.json'


def block_system(h,g5,a,layers):
    n=len(h);x=g5@h;eye=np.eye(n)
    b=2*eye+a*x
    ap=b+a*g5@x;am=b-a*g5@x
    k=np.zeros(((layers+1)*n,(layers+1)*n),complex)
    k[:n,:n]=b;k[:n,-n:]=b
    for s in range(1,layers+1):
        k[s*n:(s+1)*n,(s-1)*n:s*n]=-am
        k[s*n:(s+1)*n,s*n:(s+1)*n]=ap
    rhs=np.zeros(((layers+1)*n,n),complex);rhs[:n]=eye
    solution=np.linalg.solve(k,rhs)
    out=b@(solution[:n]-solution[-n:])
    sk,lk=np.linalg.slogdet(k)
    sa,la=np.linalg.slogdet(ap);sb,lb=np.linalg.slogdet(b)
    r=np.linalg.solve(ap,am)
    sr,lr=np.linalg.slogdet(eye+np.linalg.matrix_power(r,layers))
    return dict(eps=out,logdet=float(lk),logdet_factorized=float(layers*la+lb+lr),
        determinant_phase=base.old.cpair(sk),
        determinant_factorization_error=float(abs(lk-layers*la-lb-lr)),
        phase_factorization_error=float(abs(sk-sa**layers*sb*sr)),
        recurrence_residual=current.norm(k@solution-rhs),
        dimension=len(k),layers=layers)


def evaluate(links,e,phi,a=.23,layers=3):
    mat=base.fixed_matrices(e,phi);_,_,_,h,gap=base.kernel(links)
    g5=mat[1]@mat[1].conj().T-mat[0]@mat[0].conj().T
    q=block_system(h,g5,a,layers)
    spectral,_,_=current.regulate(h,g5,a,layers)
    direct=current.soft(spectral,mat);boundary=current.soft(q['eps'],mat)
    q.update(Wilson_gap=gap,
        rational_boundary_error=current.norm(q['eps']-spectral),
        actual_Phi_error=current.norm(boundary['Phi']-direct['Phi']),
        actual_full_N_error=current.norm(boundary['N']-direct['N']))
    q.pop('eps')
    assert max(q[key] for key in ('determinant_factorization_error','phase_factorization_error',
        'recurrence_residual','rational_boundary_error','actual_Phi_error','actual_full_N_error'))<2e-10
    return q


def run():
    links,e,phi=current.fixture()
    balanced=evaluate(links,e,phi)
    flux=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    for i,(t,x) in enumerate(base.prior.SITES):
        for mu,z in enumerate(((-1.)**t if x==1 else 1.,1j if x==1 else 1.)):
            flux[mu,i]=base.prior.rep(np.eye(3),np.eye(2),z)
    failed_seed=evaluate(flux,e,phi)
    rows=[];step=1e-4
    for angle in (-step,step):
        moved=links.copy()
        moved[1,0]=base.prior.rep(np.eye(3),np.eye(2),np.exp(1j*angle))@moved[1,0]
        data=evaluate(moved,e,phi)
        rows.append(dict(angle=angle,logdet=data['logdet'],
            factorized=data['logdet_factorized']))
    derivative=(rows[1]['logdet']-rows[0]['logdet'])/(2*step)
    return dict(date='2026-10-02',entry_round=678,not_formal_round=True,
        original_balanced=balanced,original_676_seed_obstruction=failed_seed,
        original_temporal_hypercharge_source=rows,
        bulk_logdet_source_response=float(derivative),
        source_derivative_finite_difference_diagnostic_only=True,
        original_fields_unchanged_except_declared_link_source=True,
        same_16_internal_channels=True,complete_Grassmann_identity_not_yet_proved=True,
        no_physical_time_or_RP_identity=True,
        dependency_hashes={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
            for p in ('research_note_676.md','research_note_677.md','joint_rational_physical_limit.py')})


if __name__=='__main__':
    result=run()
    if '--check' in sys.argv:assert result==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result))

