"""933: content-blind isolated links versus coherent endpoint autonomy.
Recomputable finite matrix evidence; the all-path classification is in the note.
Graph connections are mature mathematics (Kenyon, arXiv:1001.4028, section 3).
"""
from pathlib import Path
import numpy as np
import json,hashlib,argparse
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'transport_holonomy_contract_results.json'
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex);Z=np.diag([1,-1]).astype(complex)
EDGE_PAIRS=[(1,0),(3,1),(2,0),(3,2)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a))
def opnorm(a):return float(np.linalg.norm(a,2))
def unitary(h,t):
 w,v=np.linalg.eigh(h);return (v*np.exp(-1j*t*w))@v.conj().T

def rot(axis,theta):
 axis=np.asarray(axis,float);axis=axis/np.linalg.norm(axis)
 return np.cos(theta)*I-1j*np.sin(theta)*(axis[0]*X+axis[1]*Y+axis[2]*Z)

def block_diag(blocks):
 out=np.zeros((2*len(blocks),2*len(blocks)),complex)
 for i,b in enumerate(blocks):out[2*i:2*i+2,2*i:2*i+2]=b
 return out

def graph_h(edges):
 h=np.zeros((8,8),complex)
 for (v,u),k in zip(EDGE_PAIRS,edges):h[2*v:2*v+2,2*u:2*u+2]=k;h[2*u:2*u+2,2*v:2*v+2]=k.conj().T
 return h+3*np.eye(8)

def loop(edges):
 a,b,c,d=edges;return (b@a).conj().T@(d@c)

def source(a):
 v=np.zeros(8,complex);v[a]=1;return v

def task_spread(w,phase):
 a=(np.exp(1j*phase)*w+np.exp(-1j*phase)*w.conj().T)/2
 eig=np.linalg.eigvalsh(a);return float((eig[-1]-eig[0])/2)

def single_edge(k):
 h=np.block([[np.zeros((2,2)),k.conj().T],[k,np.zeros((2,2))]])
 rows=[]
 for t in [.17,.43,.91]:
  f=unitary(h,t)[2:4,:2];effect=f.conj().T@f
  rows.append(dict(t=t,probability_eigenvalues=np.linalg.eigvalsh(effect).real.tolist(),
    equal_rate_formula_error=norm(effect-np.sin(t)**2*I)))
 return rows

def run():
 edges=[I,I,I,-Z];h=graph_h(edges);hs=h-3*np.eye(8);w=loop(edges)
 assert norm(w+Z)<1e-14 and norm(h-h.conj().T)<1e-14
 link=[]
 for k in edges:
  r=single_edge(k);assert max(x['equal_rate_formula_error'] for x in r)<2e-14;link.append(r)
 # Negative control: unequal singular values filter unknown contents.
 filter_rows=single_edge(np.diag([1,2]));assert filter_rows[1]['probability_eigenvalues'][1]-filter_rows[1]['probability_eigenvalues'][0]>.3
 tstar=np.pi/2;ut=unitary(h,tstar);expected=[0.,1.];samples=[]
 for a in [0,1]:
  initial=source(a);final=ut@initial
  p=float(np.vdot(final[6:8],final[6:8]).real)
  moment=[float(np.vdot(initial,np.linalg.matrix_power(hs,k)@initial).real) for k in range(1,5)]
  full_energy=float(np.vdot(initial,h@initial).real)
  energy_error=abs(float(np.vdot(final,h@final).real)-full_energy)
  assert abs(p-expected[a])<2e-14 and energy_error<2e-13
  samples.append(dict(input_Z_eigenvalue=1-2*a,target_probability=p,
    initial_full_energy_mean=full_energy,initial_full_energy_variance=moment[1],
    initial_shifted_H_moments_1_to_4=moment,total_energy_conservation_error=energy_error))
 # Exact block certificates prove the zero-versus-one result for all real t.
 plus=hs[np.ix_([0,2,4,6],[0,2,4,6])]
 minus=hs[np.ix_([1,3,5,7],[1,3,5,7])]
 blocked_identity_error=norm(plus@plus-2*np.eye(4));assert blocked_identity_error==0
 open_formula_errors=[]
 for t in [.0,.17,.43,.91,np.pi/2]:
  u=unitary(hs,t)
  open_formula_errors.append(abs(abs(u[7,1])**2-np.sin(t)**4))
  assert abs(u[6,0])<2e-14
 assert max(open_formula_errors)<2e-14
 # Same H gives the local current; no additional independently chosen rates.
 t=.43;v=unitary(h,t)@((source(0)+source(1))/np.sqrt(2));vdot=-1j*h@v
 prob=np.array([np.vdot(v[2*j:2*j+2],v[2*j:2*j+2]).real for j in range(4)])
 derivative=np.array([2*np.vdot(v[2*j:2*j+2],vdot[2*j:2*j+2]).real for j in range(4)])
 divergence=np.zeros(4);current_rows=[]
 for (b,a),k in zip(EDGE_PAIRS,edges):
  cur=float(2*np.imag(np.vdot(v[2*b:2*b+2],k@v[2*a:2*a+2])))
  bound=float(2*opnorm(k)*np.sqrt(prob[a]*prob[b]));assert abs(cur)<=bound+2e-14
  divergence[b]+=cur;divergence[a]-=cur
  current_rows.append(dict(source=a,target=b,current=cur,content_independent_norm_bound=bound))
 continuity_error=norm(derivative-divergence);assert continuity_error<2e-14
 # Simultaneous frame changes include source, Hamiltonian, and output dictionary.
 frames=[rot([1,2,3],.19),rot([2,1,0],.31),rot([1,0,2],-.27),rot([0,1,3],.47)]
 q=block_diag(frames);ge=[frames[b]@k@frames[a].conj().T for (b,a),k in zip(EDGE_PAIRS,edges)]
 hg=graph_h(ge);gauge_error=norm(hg-q@h@q.conj().T)
 loop_gauge_error=norm(loop(ge)-frames[0]@w@frames[0].conj().T)
 gauge_prob_error=0.
 for a in [0,1]:
  v=unitary(hg,tstar)@(q@source(a));gauge_prob_error=max(gauge_prob_error,abs(float(np.vdot(v[6:],v[6:]).real)-expected[a]))
 assert gauge_error<2e-14 and loop_gauge_error<2e-14 and gauge_prob_error<2e-14
 # Central-holonomy positive case: nontrivial node frames, same scalar edge rates.
 phases=[.0,.0,.0,.37];central_edges=[frames[b]@(np.exp(1j*phase)*I)@frames[a].conj().T for (b,a),phase in zip(EDGE_PAIRS,phases)]
 hc=graph_h(central_edges);scalar_h=graph_h([np.exp(1j*p)*I for p in phases])
 central_gauge_error=norm(q.conj().T@hc@q-scalar_h)
 central_loop_error=norm(loop(central_edges)-np.exp(1j*.37)*I)
 central_arrival_spread=0.
 for t in [.17,.43,.91,tstar]:
  u=unitary(hc,t)
  for b in range(4):
   f=u[2*b:2*b+2,:2];e=f.conj().T@f
   central_arrival_spread=max(central_arrival_spread,float(np.ptp(np.linalg.eigvalsh(e))))
 assert max(central_gauge_error,central_loop_error,central_arrival_spread)<2e-14
 # Two interference quadratures matter; one can hide noncentral SU(2) transport.
 iZ=1j*Z;quadrature_spreads=[task_spread(iZ,p) for p in [0,np.pi/2]]
 assert abs(quadrature_spreads[0])<1e-14 and abs(quadrature_spreads[1]-1)<1e-14
 # Finite precision inequality: two quadrature spreads <= eps imply near-scalar.
 approximate=[]
 for theta in [.001,.02,.11]:
  wtest=np.exp(.23j)*rot([1,2,3],theta)
  a=(wtest+wtest.conj().T)/2;b=(wtest-wtest.conj().T)/(2j)
  av=np.linalg.eigvalsh(a);bv=np.linalg.eigvalsh(b)
  eps=max(np.ptp(av),np.ptp(bv))/2;c=(av[0]+av[-1])/2+1j*(bv[0]+bv[-1])/2
  distance=opnorm(wtest-c/abs(c)*I);assert 2*eps<1 and distance<=4*eps+1e-14
  approximate.append(dict(theta=theta,two_quadrature_probability_spread=float(eps),
    operator_distance_to_selected_scalar_phase=distance,analytic_bound=float(4*eps)))
 inputs=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/geometry_source_priority_review.md',HERE/'drafts/transport_holonomy_contract_working.md',STAGE/'research_note_924.md',STAGE/'research_note_925.md',STAGE/'research_note_932.md',ROOT/'research_cognition_physics/archive_467_530/research_note_512.md']
 return dict(round=933,date='2026-10-07',all_scientific_checks_passed=True,
  one_bounded_joint_classification_group=True,isolated_edge_checks=link,
  unequal_singular_value_negative_control=filter_rows,arrival_witness=samples,
  smallest_full_H_eigenvalue=float(np.linalg.eigvalsh(h).min()),
  full_H_unitarity_error=norm(ut.conj().T@ut-np.eye(8)),
  blocked_sector_squared_identity_error=blocked_identity_error,
  open_sector_sin4_formula_error=max(open_formula_errors),
  same_H_current_rows=current_rows,continuity_error=continuity_error,
  node_frame_H_covariance_error=gauge_error,loop_conjugacy_error=loop_gauge_error,
  node_frame_actual_probability_error=gauge_prob_error,
  central_connection_factorization_error=central_gauge_error,
  central_loop_scalar_error=central_loop_error,
  central_connection_arrival_effect_spread=central_arrival_spread,
  one_quadrature_can_hide_internal_information=quadrature_spreads,
  finite_precision_rows=approximate,
  all_path_quantifier_and_visible_full_state_family_are_extra_contracts=True,
  edge_rate_and_internal_transport_constrained_together=True,
  noncentral_loop_not_claimed_nonAbelian_group_derivation=True,
  initial_mean_and_variance_equal_but_full_energy_distributions_differ=True,
  full_global_unknown_state_preservation_does_not_imply_receiver_only_preservation=True,
  coordinate_dimension_metric_or_Gauss_dynamics_generated=False,
  graph_and_control_permissions_derived_from_six=False,whole_program_disproved=False,
  candidate_detail_optimization_stopped=True,full_goal_completed=False,
  source_hashes={str(p.relative_to(ROOT)):sha(p) for p in inputs})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();data=run()
 if args.write:
  with TARGET.open('x',encoding='utf-8') as f:json.dump(data,f,ensure_ascii=False,indent=2);f.write('\n')
 else:assert json.loads(TARGET.read_text('utf-8'))==data,'Recomputed result differs'
 print(json.dumps({k:v for k,v in data.items() if k in ['round','all_scientific_checks_passed','arrival_witness','smallest_full_H_eigenvalue','continuity_error','node_frame_actual_probability_error','central_connection_arrival_effect_spread','finite_precision_rows']},ensure_ascii=False,indent=2))
