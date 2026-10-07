"""875: local product truncations of the same original fixed-graph process.
The theorem covers the full original graph. Numerics reuse706's declared two-node
radial/neutral-CAR, frozen-link diagnostic; they do not simulate full Gauss fields.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
TARGET=HERE/'local_graph_process_bridge_results.json'

def matrix_function(H,G,z):
    e,v=np.linalg.eigh(H);f=np.exp(z*e)
    delta=e[:,None]-e[None,:]
    div=np.empty_like(delta,dtype=complex)
    mask=np.abs(delta)>1e-10
    div[mask]=(f[:,None]-f[None,:])[mask]/delta[mask]
    div[~mask]=(z*np.exp(z*(e[:,None]+e[None,:])/2))[~mask]
    M=(v*f)@v.conj().T
    dM=v@(div*(v.conj().T@G@v))@v.conj().T
    return M,dM
def thermal(H,G,beta):
    # A fixed shift at this evaluation improves numerical range.
    shift=float(np.linalg.eigvalsh(H).min())
    M,dM=matrix_function(H-shift*np.eye(len(H)),G,-beta)
    Z=float(np.trace(M).real)
    rho=M/Z
    drho=dM/Z-rho*np.trace(dM)/Z
    return rho,drho
def solve(d,bands,s=0.):
    r=d['r'];nloc=d['nloc'];rank=4*bands
    keep=np.array([a*nloc+b for a in range(rank) for b in range(rank)])
    parts=[a[np.ix_(keep,keep)] for a in d['parts']]
    H=sum(np.exp(k*s)*a for k,a in zip(d['powers'],parts))
    G=sum(k*np.exp(k*s)*a for k,a in zip(d['powers'],parts))
    beta=1.2;hbar=.7
    rho,drho=thermal(H,G,beta)
    U1,dU1=matrix_function(H,G,-1j*.12/hbar)
    U2,dU2=matrix_function(H,G,-1j*.19/hbar)
    occ=np.diag(np.tile([0.,1.,0.,1.],bands))
    P=np.kron(occ,np.eye(rank));Q=np.kron(np.eye(rank),occ)
    probs=[];prep=[];dyn=[]
    for p in (P,np.eye(len(H))-P):
        for q in (Q,np.eye(len(H))-Q):
            K=q@U2@p@U1
            dK=q@(dU2@p@U1+U2@p@dU1)
            probs.append(float(np.trace(K@rho@K.conj().T).real))
            prep.append(float(np.trace(K@drho@K.conj().T).real))
            dyn.append(float(2*np.trace(dK@rho@K.conj().T).real))
    return dict(keep=keep,H=H,G=G,rho=rho,drho=drho,prob=np.array(probs),
        prep=np.array(prep),dyn=np.array(dyn),total=np.array(prep)+np.array(dyn),
        source_mean=float(np.trace(rho@G).real),
        record_completeness_error=float(abs(sum(probs)-1)),
        unitarity_error=float(max(np.linalg.norm(U1.conj().T@U1-np.eye(len(H))),
                                  np.linalg.norm(U2.conj().T@U2-np.eye(len(H))))))
def trace_norm(A):
    return float(np.abs(np.linalg.eigvalsh((A+A.conj().T)/2)).sum())
def run():
    with ResearchRuntime(Layout()).installed():
        import joint_local_history_source_limit as inherited
        d=inherited.data()
    full=solve(d,d['r']);rows=[];source_fd=[]
    for k in range(1,d['r']+1):
        x=full if k==d['r'] else solve(d,k)
        keep=x['keep'];embed=np.zeros_like(full['rho']);embed[np.ix_(keep,keep)]=x['rho']
        embed_d=np.zeros_like(full['drho']);embed_d[np.ix_(keep,keep)]=x['drho']
        tail=1-float(np.trace(full['rho'][np.ix_(keep,keep)]).real)
        eps=1e-5
        plus=solve(d,k,eps);minus=solve(d,k,-eps)
        fd=(plus['prob']-minus['prob'])/(2*eps)
        err=float(np.max(np.abs(fd-x['total'])))
        assert err<2e-7
        assert x['record_completeness_error']<1e-11 and x['unitarity_error']<1e-10
        assert np.min(x['prob'])>=0 and abs(x['total'].sum())<1e-10
        row=dict(local_boson_bands=k,dimension=len(keep),
          own_Gibbs_trace_distance=trace_norm(embed-full['rho']),
          Gibbs_source_derivative_trace_distance=trace_norm(embed_d-full['drho']),
          original_state_omitted_weight=tail,
          joint_record_probabilities=x['prob'].tolist(),
          preparation_response=x['prep'].tolist(),dynamical_response=x['dyn'].tolist(),
          total_geometry_response=x['total'].tolist(),
          record_error=float(np.max(np.abs(x['prob']-full['prob']))),
          total_geometry_response_error=float(np.max(np.abs(x['total']-full['total']))),
          geometry_source_mean=x['source_mean'],
          source_mean_error=abs(x['source_mean']-full['source_mean']),
          own_Gibbs_and_real_history_source_finite_difference_error=err,
          record_completeness_error=x['record_completeness_error'],
          unitarity_error=x['unitarity_error'])
        rows.append(row)
    # Actual Hamiltonian is compressed, not replaced by the comparison operator.
    two=solve(d,2);keep=two['keep'];outside=np.setdiff1d(np.arange(len(d['H'])),keep)
    leakage=float(np.linalg.norm(d['H'][np.ix_(outside,keep)]))
    assert leakage>1e-3 and rows[-1]['record_error']==0
    # Product truncation commutes with the full CAR matrix algebra, whereas H does not.
    comparison_diagonal=d['full_energy']
    p=np.zeros(len(d['H']));p[keep]=1
    nloc=d['nloc'];N=np.kron(np.diag(np.tile([0.,1.,0.,1.],d['r'])),np.eye(nloc))
    comm=float(np.linalg.norm(p[:,None]*N-N*p[None,:]))
    assert comm==0
    # Source deletion negative control: same dynamics, frozen preparation response.
    missing_preparation=float(np.max(np.abs(two['prep'])))
    assert missing_preparation>1e-4
    return dict(round=875,date='2026-10-06',formal_reports=875,
        cumulative_numbered_groups=3660,fresh_numbered_groups=1,all_checks_passed=True,
        original_fixed_graph_theorem_is_analytic=True,
        numerical_scope='Original706 frozen-link aligned radial and neutral-CAR 256-dimensional diagnostic; not a full Gauss or spatial-continuum simulation.',
        inherited_dimension=len(d['H']),inherited_CAR_error=d['CAR_error'],
        inherited_minimum_F=d['minimum_F'],source_geometry_powers=list(d['powers']),
        original_hamiltonian_cutoff_leakage=leakage,
        product_projection_CAR_commutator=comm,
        omitted_preparation_response_defect=missing_preparation,
        product_cutoff_checks=rows,
        product_preserves_original_support_and_Gauss_analytically=True,
        fixed_graph_own_thermal_and_first_source_real_history_limit_proved=True,
        second_real_time_source_limit_for_this_projection_proved=False,
        strict_one_sided_finite_compression_claimed=False,
        same_original_Q_to_continuous_E_mapping_completed=False,
        spatial_refinement_uniform_bounds_proved=False,full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    data=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert data==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(data,ensure_ascii=False,indent=2))
