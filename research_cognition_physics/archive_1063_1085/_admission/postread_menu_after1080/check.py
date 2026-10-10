"""One fixed post-read source menu diagnostic. Count-zero STOP; no formal round.
All four Pauli coefficients are linear operator evaluations, not preparations
of four copies of an unknown probe. Full DG poststates and graph coherences
remain. Floating ranks are not proofs.
"""
from pathlib import Path
import argparse,json,math,sys,hashlib
import numpy as np
ROOT=Path(r'D:\workspace\AGI的哲学思考')
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime

def tensor(items):
 x=np.ones((1,1),complex)
 for a in items:x=np.kron(x,a)
 return x

def cross(a):
 x,y,z=a
 return np.array([[0,-z,y],[z,0,-x],[-y,x,0.]])

def run():
 with ResearchRuntime(Layout(),'late').installed():
  import correlated_reference_local_access as source
  import relational_direction_distance_response as reader
  H=source.frame.system()[1];K=source.prior.s12();Xi=source.pure_columns()
  h,V=np.linalg.eigh(H)
  def uh(t):return (V*np.exp(-1j*t*h))@V.conj().T
  I=np.eye(2,dtype=complex);pauli=(I,*reader.P)
  raw=np.array([tensor([I,s,I,I+pauli[1],I+pauli[2],I+pauli[3]])/64 for s in pauli])
  M=np.kron(tensor([I,(I+pauli[3])/2,I,I,I,I]),np.eye(6))
  tau=1e-10;p0=(.2,.4,math.pi/2)
  def word(p,derivative=False):
   t,u,theta=p
   ev,vec=np.linalg.eigh(tau*H+theta*K)
   P=(vec*np.exp(-1j*ev))@vec.conj().T
   L,R=uh(u),uh(t);W=L@P@R
   if not derivative:return W
   midpoint=(ev[:,None]+ev[None,:])/2;difference=ev[:,None]-ev[None,:]
   divided=-1j*np.exp(-1j*midpoint)*np.sinc(difference/(2*math.pi))
   dP=vec@(divided*(vec.conj().T@K@vec))@vec.conj().T
   return W,np.array([W@(-1j*H),(-1j*H)@W,L@dP@R])
  W,dW=word(p0,True)
  col=(W@Xi).reshape(64,6,28)
  gamma=np.einsum('dga,dha->gh',col,col.conj())/384
  U=uh(.7);O=U.conj().T@M@U
  B=np.einsum('diej,ked->kij',O.reshape(64,6,64,6),raw)
  bo=np.einsum('kij,ji->k',B[1:],gamma).real
  # Independent original preparation, then retain the actual z=0 poststate.
  rho=[]
  for r in raw:
   x=U@np.kron(r,gamma)@U.conj().T
   rho.append(M@x@M)
  rho=np.array(rho)
  probabilities=np.array([np.trace(x).real for x in rho])
  expected_p=np.einsum('kij,ji->k',B,gamma).real
  y=np.zeros((3,4));J=np.zeros((3,4,3))
  for j in range(3):
   A=np.kron(np.eye(64),B[j+1]);pulled=W.conj().T@A@W
   y[j]=np.einsum('ij,mji->m',pulled,rho).real
   for k in range(3):
    dO=dW[k].conj().T@A@W+W.conj().T@A@dW[k]
    J[j,:,k]=np.einsum('ij,mji->m',dO,rho).real
  a=y-bo[:,None]*probabilities[None,:]
  matrix=np.vstack([np.column_stack([J[:,mu,:],cross(a[:,mu]),-a[:,mu]]) for mu in range(4)])
  left,sv,right=np.linalg.svd(matrix,full_matrices=True)
  threshold=1e-10*sv[0];rank=int(np.sum(sv>threshold));null=right[rank:].T
  no_gain=matrix[:,:6]
  sv_no_gain=np.linalg.svd(no_gain,compute_uv=False)
  asp=np.linalg.svd(a,compute_uv=False)
  normalized_J=(a[:,1:]*probabilities[0]-a[:,0,None]*probabilities[None,1:])/probabilities[0]**2
  sv_input=np.linalg.svd(normalized_J,compute_uv=False)
  postgraph=np.trace(rho[0].reshape(64,6,64,6),axis1=0,axis2=2)/probabilities[0]
  # Finite directional derivative diagnoses only the formula at the one point.
  direction=np.array([.6,-.3,.2]);pred=np.einsum('jma,a->jm',J,direction)
  fd=[]
  for step in (1e-4,5e-5):
   out=[]
   for sign in (-1,1):
    ww=word(np.array(p0)+sign*step*direction)
    yy=np.array([np.einsum('ij,mji->m',ww.conj().T@np.kron(np.eye(64),bb)@ww,rho).real for bb in B[1:]])
    out.append(yy)
   estimate=(out[1]-out[0])/(2*step)
   fd.append({'step':step,'max_error_from_Frechet':float(np.max(abs(estimate-pred)))})
  saved=json.loads((ROOT/'research_cognition_physics/archive_1063_/1080/independent_results.json').read_text(encoding='utf-8'))
  assert max(abs(probabilities-expected_p))<2e-13
  assert max(abs(expected_p-np.array(saved['sources'][0]['effect_coefficients'])))<2e-13
  assert np.max(abs(rho-rho.conj().transpose(0,2,1)))<2e-13
  assert fd[-1]['max_error_from_Frechet']<2e-10
  source_trace=float(np.trace(gamma).real)
  # Physical inputs include the six cardinal pure probes and maximally mixed.
  physical_post_min=[]
  for r in [np.zeros(3)]+[s*np.eye(3)[k] for k in range(3) for s in (-1,1)]:
   post=rho[0]+np.einsum('i,ijk->jk',r,rho[1:])
   physical_post_min.append(float(np.linalg.eigvalsh(post)[0]))
  return {
   'status':'diagnostic_passed','source_parameters':list(p0),'append_word_base':list(p0),
   'finite_pulse_tau':tau,'reader_time':.7,'raw_outcome':0,'source_trace':source_trace,
   'source_graph_coherence':float(np.linalg.norm(gamma-np.diag(np.diag(gamma)))),
   'post_graph_coherence_at_mixed_probe':float(np.linalg.norm(postgraph-np.diag(np.diag(postgraph)))),
   'anchor_b':bo.tolist(),'probability_Pauli_coefficients':probabilities.tolist(),
   'unnormalized_y_j_mu':y.tolist(),'relative_a_j_mu':a.tolist(),
   'Jacobian_j_mu_parameter':J.tolist(),
   'common_linear_system_12x7':matrix.tolist(),
   'variable_order':['dt','du','dtheta','omega_x','omega_y','omega_z','kappa'],
   'equation':'J_mu dv + cross(a_mu) omega - a_mu kappa = 0',
   'rank_diagnostic':rank,'singular_values':sv.tolist(),'rank_threshold':float(threshold),
   'nullspace_columns':null.tolist(),'smallest_right_singular_vector':right[-1].tolist(),
   'no_gain_singular_values':sv_no_gain.tolist(),
   'relative_a_span_singular_values':asp.tolist(),
   'normalized_relative_input_derivative_singular_values':sv_input.tolist(),
   'finite_difference_diagnostics':fd,
   'physical_post_minimum_eigenvalues':physical_post_min,
   'fixed_source_probability_max_residual':float(np.max(abs(probabilities-expected_p))),
   'scope':{'one_source_one_branch_one_append_context':True,'graph_reset':False,
            'unknown_probe_estimated_or_used_for_control':False,
            'full_graph_coherences_retained':True,'gain_constancy_proved_here':False,
            'floating_rank_is_rigorous_certificate':False,
            'all_history_or_all_feedback_strategies_excluded':False,
            'spatial_dimension_or_cognitive_program_refuted':False},
   'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
   'fixed_dependency_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in (
    'research_cognition_physics/archive_467_530/480/correlated_reference_local_access.py',
    'research_cognition_physics/archive_467_530/474/relational_direction_distance_response.py',
    'research_cognition_physics/archive_1063_/1077/proof.md',
    'research_cognition_physics/archive_1063_/1080/independent_results.json')},
  }
if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--check',action='store_true',help='Read-only recomputation and comparison with saved results.json')
 args=parser.parse_args()
 actual=run()
 if args.check:
  expected=json.loads(Path(__file__).with_name('results.json').read_text(encoding='utf-8'))
  if actual!=expected:
   keys=sorted(k for k in set(actual)|set(expected) if actual.get(k)!=expected.get(k))
   raise SystemExit('Saved result mismatch in: '+', '.join(keys))
  print(json.dumps({'status':'saved_diagnostic_reproduced','scientific_round_increment':0,
                    'rank_diagnostic':actual['rank_diagnostic'],
                    'smallest_singular_value':actual['singular_values'][-1]},indent=2))
 else:
  print(json.dumps(actual,ensure_ascii=False,indent=2))
