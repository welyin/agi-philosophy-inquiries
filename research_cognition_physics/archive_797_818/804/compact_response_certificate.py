"""804: analytic infinite smoothing operator calibrates response certificates.

This is an error-budget diagnostic, NOT a discretization of the original PDE.
The exact infinite tail bound, rather than a larger numerical matrix, certifies
the truncation. The larger matrix is only an independent numerical check.
"""
from pathlib import Path
import argparse, json
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'compact_response_certificate_results.json'

def block(n):
    j=np.arange(n)
    v=np.sqrt(1-.4**2)*.4**j
    return np.diag(.12*(-.5)**j)+.3*np.outer(v,v)

def tail_bound(n):
    # Exact upper bound for ||L - Q_n L Q_n|| on the infinite Hilbert space.
    return .12*.5**n+.3*(2*.4**n+.4**(2*n))

def run():
    reference=block(80)
    large_norm=float(np.max(abs(np.linalg.eigvalsh(reference))))
    rows=[]
    for n in (2,4,8,16):
        b=block(n)
        # A known Hermitian numerical error; its Frobenius norm bounds op norm.
        j=np.arange(1,n+1,dtype=float)
        e=np.outer(j,j)*1e-5/(j@j)
        computed=b+e
        ev,u=np.linalg.eigh(computed)
        k=int(np.argmax(abs(ev)));f=u[:,k]
        error=float(np.linalg.norm(e,'fro'))
        observed=float((f.conj()@b@f).real)
        lower=float(abs(ev[k])-error)
        assert lower>0 and abs(observed)+1e-14>=lower
        norm_radius=error+tail_bound(n)
        assert abs(large_norm-abs(ev[k]))<=norm_radius+1e-14
        padded=np.zeros_like(reference);padded[:n,:n]=b
        finite_tail=float(np.linalg.norm(reference-padded,2))
        assert finite_tail<=tail_bound(n)+1e-14
        rows.append(dict(modes=n,computed_eigenvalue=float(ev[k]),
            matrix_error_upper_bound=error,witness_lower_bound=lower,
            direct_rayleigh_value=observed,
            infinite_tail_upper_bound=tail_bound(n),
            global_norm_interval=[float(max(0,abs(ev[k])-norm_radius)),
                                  float(abs(ev[k])+norm_radius)],
            large_matrix_tail_check=finite_tail))
    diagonals=[dict(mode=j+1,response=float(block(j+1)[j,j])) for j in (0,3,7,15,31)]
    assert abs(diagonals[-1]['response'])<1e-10
    hidden=np.zeros((12,12));hidden[8,8]=.7
    assert np.count_nonzero(hidden[:8,:8])==0 and np.linalg.norm(hidden,2)==.7
    # Recover the Hermitian response from a complex ordered pairing.
    b=block(8);a=.5j*b+np.diag(np.arange(8)/13.)
    response=(a-a.conj().T)/1j
    extraction=float(np.max(abs(response-b)))
    assert extraction<1e-14
    return dict(round=804,all_checks_passed=True,
        scope='Analytic infinite smoothing operator; not the original joint PDE or its state.',
        finite_block_certificates=rows,weakly_null_mode_responses=diagonals,
        hidden_mode_counterexample=dict(first_eight_modes_zero=True,ninth_response=.7),
        hermitian_response_extraction_residual=extraction,
        analytic_tail_proof='D tail <= .12*.5**n; rank-one tail <= .3*(2*.4**n+.4**(2*n)).',
        original_PDE_evaluated=False,all_checks_are_numerical_calibration=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
