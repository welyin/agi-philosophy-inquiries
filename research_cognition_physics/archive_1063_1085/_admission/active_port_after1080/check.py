"""One active two-port diagnostic; count-zero only, no control/readout optimization."""
from pathlib import Path
import argparse,sys,math,json,hashlib
import numpy as np
ROOT=Path(r'D:\workspace\AGI的哲学思考')
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime

def tensor(items):
 o=np.ones((1,1),complex)
 for a in items:o=np.kron(o,a)
 return o

def realvec(a):return np.r_[a.real.ravel(),a.imag.ravel()]

def run():
 with ResearchRuntime(Layout(),'late').installed():
  import correlated_reference_local_access as source
  import relational_direction_distance_response as reader
  trees,H,*_=reader.system();K=source.prior.s12();Xi=source.pure_columns()
  kd=K[::6,::6]
  I=np.eye(2,dtype=complex);pauli=(I,*reader.P)
  raw=np.array([tensor([I,s,I,I+pauli[1],I+pauli[2],I+pauli[3]])/64 for s in pauli])
  raw2=np.array([kd@r@kd for r in raw])
  M1=np.kron(tensor([I,(I+pauli[3])/2,I,I,I,I]),np.eye(6))
  M2=K@M1@K
  vals,vec=np.linalg.eigh(H)
  def uh(t):return (vec*np.exp(-1j*t*vals))@vec.conj().T
  U=uh(.7)
  B=[]
  for rawi,M in ((raw,M1),(raw2,M2)):
   O=U.conj().T@M@U
   B.append(np.einsum('diej,ked->kij',O.reshape(64,6,64,6),rawi))
  B1,B2=B
  canonical=[frozenset(tuple(sorted(e)) for e in tree) for tree in trees]
  T=np.zeros((6,6),dtype=int);perm=[]
  def swaplabel(a):return 3-a if a in (1,2) else a
  for g,tree in enumerate(trees):
   renamed=frozenset(tuple(sorted((swaplabel(a),swaplabel(b)))) for a,b in tree)
   h=canonical.index(renamed);T[h,g]=1;perm.append(h)
  joint=np.kron(kd,T)
  exact_joint_covariance=np.array_equal(joint@H@joint.T,H)
  exact_T_involution=np.array_equal(T@T,np.eye(6,dtype=int))
  graph_covariance_residual=max(float(np.max(abs(bb-T@ba@T.T))) for ba,bb in zip(B1,B2))
  span=np.column_stack([realvec(np.eye(6))]+[realvec(bb) for bb in B1[1:]])
  target=np.column_stack([realvec(bb) for bb in B2[1:]])
  coeff,_,rank,sv=np.linalg.lstsq(span,target,rcond=None)
  residual=target-span@coeff
  allsv=np.linalg.svd(np.column_stack([span,target]),compute_uv=False)
  tau=1e-10;p0=(.2,.4,math.pi/2)
  vv,vc=np.linalg.eigh(tau*H+p0[2]*K)
  P=(vc*np.exp(-1j*vv))@vc.conj().T;L=uh(p0[1]);R=uh(p0[0]);W=L@P@R
  delta=vv[:,None]-vv[None,:];mid=(vv[:,None]+vv[None,:])/2
  kernel=-1j*np.exp(-1j*mid)*np.sinc(delta/(2*math.pi))
  dP=vc@(kernel*(vc.conj().T@K@vc))@vc.conj().T
  dW=[W@(-1j*H),(-1j*H)@W,L@dP@R]
  Y=(W@Xi).reshape(64,6,28)
  gamma=np.einsum('dga,dha->gh',Y,Y.conj())/384
  dg=[]
  for dw in dW:
   dY=(dw@Xi).reshape(64,6,28)
   a=np.einsum('dga,dha->gh',dY,Y.conj())/384
   dg.append(a+a.conj().T)
  dg=np.array(dg)
  e=np.array([np.einsum('mgh,hg->m',b,gamma).real for b in B])
  je=np.array([np.einsum('mgh,khg->mk',b,dg).real for b in B])
  transition=je[1,1:]@np.linalg.inv(je[0,1:])
  ts=np.linalg.svd(transition,compute_uv=False)
  scale=float(np.sqrt(np.trace(transition.T@transition)/3))
  conformality_residual=float(np.linalg.norm(transition.T@transition-scale**2*np.eye(3)))
  assert exact_joint_covariance and exact_T_involution and graph_covariance_residual<3e-13
  return {
   'status':'diagnostic_passed','reader_time':.7,'source_parameters':list(p0),'finite_pulse_tau':tau,
   'graph_internal_swap_permutation':perm,'graph_T':T.tolist(),
   'integer_joint_H_covariance':exact_joint_covariance,'integer_T_squared_identity':exact_T_involution,
   'B2_equals_T_B1_T_max_residual':graph_covariance_residual,
   'all_gamma_affine_fit_coefficients_rows_I_B1xyz':coeff.tolist(),
   'all_gamma_affine_fit_residual_frobenius_by_component':np.linalg.norm(residual,axis=0).tolist(),
   'all_gamma_affine_fit_max_residual':float(np.max(abs(residual))),
   'I_B1_span_singular_values':sv.tolist(),'I_B1_B2_span_singular_values':allsv.tolist(),
   'B1_real':B1.real.tolist(),'B1_imag':B1.imag.tolist(),'B2_real':B2.real.tolist(),'B2_imag':B2.imag.tolist(),
   'effect_centers_port1_port2':e.tolist(),'effect_Jacobians_port1_port2':je.tolist(),
   'fixed_probe_diagnostics':{
      'plus_X_probabilities':[float(x[0]+x[1]) for x in e],
      'plus_X_probability_port2_minus_port1':float((e[1,0]+e[1,1])-(e[0,0]+e[0,1])),
      'mixed_probe_probability_port2_minus_port1':float(e[1,0]-e[0,0]),
      'plus_Y_minus_mixed_contrast_port2_minus_port1':float(e[1,2]-e[0,2])},
   'bloch_Jacobian_singular_values':[np.linalg.svd(x[1:],compute_uv=False).tolist() for x in je],
   'local_transition_J2_J1_inverse':transition.tolist(),'local_transition_singular_values':ts.tolist(),
   'local_transition_determinant':float(np.linalg.det(transition)),
   'common_scale_least_squares':scale,'local_conformality_residual':conformality_residual,
   'source_graph_swap_distance_frobenius':float(np.linalg.norm(gamma-T@gamma@T.T)),
   'physical_graph_changed_by_protocol':False,'leaf_reference_states_changed':False,
   'unknown_probe_copied':False,'mathematical_graph_T_is_physical_protocol':False,
   'scopes':{'all_gamma_affine_test_is_universal_operator_span_question':True,
      'all_gamma_fit_failure_implies_tiny_cube_affine_failure':False,
      'local_J_transition_tests_constant_rotation_gain_necessary_condition_only':True,
      'floating_diagnostics_are_rigorous_certificate':False,
      'two_port_failure_excludes_other_contact_tasks':False,
      'spatial_dimension_derived_or_refuted':False},
   'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
   'fixed_dependency_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in (
    'research_cognition_physics/archive_467_530/480/correlated_reference_local_access.py',
    'research_cognition_physics/archive_467_530/474/relational_direction_distance_response.py',
    'research_cognition_physics/archive_1063_/1077/proof.md')}
  }
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--check',action='store_true',help='Read-only recomputation and comparison against results.json')
 args=parser.parse_args();actual=run()
 if args.check:
  saved=json.loads(Path(__file__).with_name('results.json').read_text(encoding='utf-8'))
  if actual!=saved:
   keys=sorted(k for k in set(actual)|set(saved) if actual.get(k)!=saved.get(k))
   raise SystemExit('Saved diagnostic mismatch in: '+', '.join(keys))
  print(json.dumps({'status':'saved_diagnostic_reproduced','scientific_round_increment':0,
    'integer_joint_H_covariance':actual['integer_joint_H_covariance'],
    'local_transition_singular_values':actual['local_transition_singular_values']},indent=2))
 else:print(json.dumps(actual,ensure_ascii=False,indent=2))
