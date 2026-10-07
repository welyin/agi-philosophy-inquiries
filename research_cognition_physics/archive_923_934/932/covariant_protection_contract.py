"""932: finite covariant encoding/known erasure. No continuum or physical geometry.
Generalized W encoding is a known construction, Faist et al. arXiv:1902.07714 VII.2.
--write creates one result; default recomputes and compares without overwriting.
"""
from pathlib import Path
import numpy as np
import json, hashlib, argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'covariant_protection_contract_results.json'
I2=np.eye(2,dtype=complex)
PAULI=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.diag([1,-1]).astype(complex)]
EMBED=np.array([[0,0],[1,0],[0,1]],complex)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a))
def trnorm(a):return float(np.linalg.svd(a,compute_uv=False).sum())
def positive_sqrt(a):
 w,u=np.linalg.eigh((a+a.conj().T)/2)
 assert w.min()>-2e-13
 return (u*np.sqrt(np.maximum(w,0)))@u.conj().T

def encoding(weights):
 n=len(weights);v=np.zeros((3**n,2),complex)
 for i,p in enumerate(weights):
  for a in range(2):v[(a+1)*3**(n-1-i),a]=np.sqrt(p)
 return v

def local_action(v,u,i,n):
 a=np.moveaxis(v.reshape((3,)*n+(2,)),i,0)
 b=np.tensordot(u,a,axes=(1,0))
 return np.moveaxis(b,0,i).reshape(-1,2)

def choi(kraus):
 vectors=[k.T.reshape(-1)/np.sqrt(2) for k in kraus]
 return sum(np.outer(v,v.conj()) for v in vectors)

def case(weights):
 weights=np.asarray(weights,float);n=len(weights);assert abs(weights.sum()-1)<1e-14
 v=encoding(weights);iso=norm(v.conj().T@v-I2)
 cov=0.
 for axis,theta,phase in [(np.array([1,0,0]),.37,.13),(np.array([0,1,0]),.73,-.19),(np.array([1,2,3])/np.sqrt(14),1.17,.29)]:
  gen=sum(axis[k]*PAULI[k] for k in range(3))
  u=np.exp(1j*phase)*(np.cos(theta)*I2-1j*np.sin(theta)*gen)
  g=np.eye(3,dtype=complex);g[1:,1:]=u
  acted=v.copy()
  for i in range(n):acted=local_action(acted,g,i,n)
  cov=max(cov,norm(acted-v@u))
 digits=np.indices((3,)*n).reshape(n,-1)
 total_charge=(digits==2).sum(axis=0);total_occupation=(digits!=0).sum(axis=0)
 charge_error=norm(total_charge[:,None]*v-v@np.diag([0,1]))
 occupation_error=norm(total_occupation[:,None]*v-v)
 pv=np.swapaxes(v.reshape((3,)*n+(2,)),0,n-1).reshape(-1,2)
 perm_error=norm(pv-v) if np.ptp(weights)<1e-14 else None
 ideal_vector=EMBED.T.reshape(-1)/np.sqrt(2);ideal=np.outer(ideal_vector,ideal_vector.conj())
 lost=np.kron(I2/2,np.diag([1,0,0]));rows=[]
 for i,p in enumerate(weights):
  e=np.moveaxis(v.reshape((3,)*n+(2,)),i,0).reshape(3,3**(n-1),2)
  remaining=np.delete(weights,i)/(1-p);vbar=encoding(remaining)
  decoder=np.zeros((3,3**(n-1)),complex);decoder[0,0]=1;decoder[1:]=vbar.conj().T
  ks=[decoder@ea for ea in e]
  support=max(norm(ea-decoder.conj().T@decoder@ea) for ea in e)
  j=choi(ks);expected=(1-p)*ideal+p*lost
  tp=norm(sum(k.conj().T@k for k in ks)-I2)
  rho=[e[:,:,a]@e[:,:,a].conj().T for a in range(2)]
  delta=trnorm(rho[1]-rho[0])/2
  rootf=trnorm(positive_sqrt(rho[0])@positive_sqrt(rho[1]))
  a=sum(np.outer(ea[:,0],ea[:,1].conj()) for ea in e)
  cross_norm=trnorm(a)
  fe=float(np.vdot(ideal_vector,j@ideal_vector).real)
  half_trace=trnorm(j-ideal)/2
  error_bound=(1-np.sqrt(max(0.,1-delta**2)))/2
  record=dict(erased_site=i,probability=float(p),local_trace_distance=delta,
   root_fidelity=rootf,complement_cross_trace_norm=cross_norm,
   recovery_support_error=support,trace_preserving_error=tp,
   choi_formula_error=norm(j-expected),minimum_choi_eigenvalue=float(np.linalg.eigvalsh(j).min()),
   entanglement_fidelity=fe,entanglement_infidelity=1-fe,
   entangled_witness_half_trace_distance=half_trace,
   all_task_half_diamond_error_analytical=float(p),
   general_recovery_infidelity_lower_bound=float(error_bound))
  assert abs(delta-p)<2e-13 and abs(rootf-(1-p))<2e-13
  assert abs(cross_norm-rootf)<2e-13 and abs(fe-(1-p))<2e-13
  assert abs(half_trace-p)<2e-13 and record['choi_formula_error']<2e-13
  assert support<2e-13 and tp<2e-13
  assert record['minimum_choi_eigenvalue']>-2e-13 and 1-fe+2e-13>=error_bound
  rows.append(record)
 assert iso<2e-13 and cov<2e-13 and charge_error<2e-13 and occupation_error<2e-13
 if perm_error is not None:assert perm_error<2e-13
 assert abs(sum(x['local_trace_distance'] for x in rows)-1)<2e-13
 return dict(n=n,weights=weights.tolist(),physical_tensor_dimension=3**n,
  isometry_error=iso,full_U2_covariance_sample_error=cov,
  exact_covariance_general_proof_in_note=True,total_generator_intertwining_error=charge_error,
  total_occupation_error=occupation_error,uniform_permutation_error=perm_error,
  logical_generator_span=1.,sum_local_generator_spans=float(n),
  weighted_leakage_sum=sum(x['local_trace_distance'] for x in rows),erasures=rows)

def internal_swap_check():
 s=np.zeros((9,9),complex)
 for a in range(3):
  for b in range(3):s[b*3+a,a*3+b]=1
 num=np.diag([0,1,1]);nt=np.kron(num,np.eye(3))+np.kron(np.eye(3),num)
 charge=np.diag([0,0,1]);qt=np.kron(charge,np.eye(3))+np.kron(np.eye(3),charge)
 blank=np.array([1,0,0],complex)
 v=s@np.kron(np.eye(3),blank[:,None])
 expected=np.kron(blank[:,None],np.eye(3))
 result=dict(swap_unitarity_error=norm(s.conj().T@s-np.eye(9)),
  total_occupation_commutator=norm(s@nt-nt@s),total_charge_commutator=norm(s@qt-qt@s),
  unknown_site_transferred_into_internal_inaccessible_register_error=norm(v-expected))
 assert max(result.values())<1e-14
 return result

def run():
 cases=[case(np.ones(n)/n) for n in (2,3,4,6)]
 cases.append(case([.1,.2,.3,.4]))
 src=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/operation_selection_review.md',HERE/'drafts/covariant_protection_working.md',STAGE/'research_note_924.md',STAGE/'research_note_829.md',STAGE/'research_note_830.md',STAGE/'research_note_929.md']
 return dict(round=932,date='2026-10-06',all_scientific_checks_passed=True,
  one_bounded_contract_test_group=True,cases=cases,internal_erasure=internal_swap_check(),
  analytic_scaling_examples=[dict(n=n,half_diamond_error=1/n,logical_span=1,sum_local_spans=n) for n in (8,32,128,1000)],
  exact_erasure_plus_nontrivial_transversal_continuity_excluded_conditionally=True,
  approximate_accepted_contract_has_finite_dimension_for_each_nonzero_accuracy=True,
  generator_span_is_not_consumed_energy=True,
  mature_theorem_and_generalized_W_construction_not_claimed_original=True,
  active_permissions_and_single_share_erasure_are_extra_inputs=True,
  prepared_blanks_access_rules_and_encoding_decoding_permissions_are_inputs=True,
  erased_information_kept_inside_overall_system=True,
  natural_autonomous_controller_derived=False,active_permissions_forced_by_six_protocols=False,
  all_six_protocols_physically_implemented=False,spacetime_dimension_gauge_group_or_GR_derived=False,
  exact_strong_contract_failure_is_whole_program_failure=False,full_goal_completed=False,
  branch_stops_without_code_optimization=True,
  source_hashes={str(p.relative_to(ROOT)):sha(p) for p in src})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();data=run()
 if args.write:
  with TARGET.open('x',encoding='utf-8') as f:json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
 else:
  before=json.loads(TARGET.read_text('utf-8'));assert before==data,'Recomputed result differs'
 print(json.dumps(dict(round=data['round'],passed=data['all_scientific_checks_passed'],
  cases=[dict(n=c['n'],max_channel_error=max(x['choi_formula_error'] for x in c['erasures']),max_recovery_error=max(x['entanglement_infidelity'] for x in c['erasures'])) for c in data['cases']],
  swap=data['internal_erasure']),ensure_ascii=False,indent=2))
