"""948: actual 946 portal kernel connects the 947 event source to SM vertices.

Tree-level, spacelike quadratic response only. This is NOT a finite-time
matching certificate for the full SM/Einstein parent or a physical data fit.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'neutral_matter_bridge_results.json'

def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))

def run():
    old=read(STAGE/'946/portal_common_process_results.json')
    protocol=read(STAGE/'947/protocol_field_transport_results.json')
    assert old['all_scientific_checks_passed'] and protocol['all_scientific_checks_passed']
    K=np.array(old['parameters']['mass_matrix'])
    lam=old['parameters']['lambda_H'];v=old['parameters']['v']
    alpha=old['parameters']['portal_alpha'];mu=old['parameters']['mu']
    assert np.allclose(K,[[mu**2+2*lam*alpha**2,2*lam*alpha*v],
                         [2*lam*alpha*v,2*lam*v*v]])
    mass2,O=np.linalg.eigh(K)
    # Numerical electroweak/Yukawa values are diagnostic physical inputs.
    gs=.5;g2=.65;yf=.7
    mw=g2*v/2;mf=yf*v/math.sqrt(2)
    cw=2*mw**2/v;cf=mf/v
    # Energy-current convention: source densities absorb the Lorentz signs.
    C=np.array([[gs,0,0],[0,cw,cf]])
    null=np.array([0,cf,-cw])
    # b = (-1)^x in the frozen first-event instrument.
    inputs=protocol['joint_instrument']['input_examples']
    bmeans=[r['first_event_probabilities'][0]-r['first_event_probabilities'][1] for r in inputs]
    delta_b=bmeans[0]-bmeans[1]
    rows=[];inv_error=0.;fd_error=0.;stationarity_error=0.
    for q in (0.,.25,.5,1.,2.,4.):
        A=q*q*np.eye(2)+K
        D=np.linalg.inv(A)
        spectral=(O/(q*q+mass2)[None,:])@O.T
        direct_off=-K[0,1]/np.linalg.det(A)
        inv_error=max(inv_error,norm(D-spectral),abs(D[0,1]-direct_off))
        response=-C.T@D@C
        # Independent elimination uses solve(A, C@currents) at every call.
        def energy(j):
            phi=-np.linalg.solve(A,C@j)
            return float(.5*phi@A@phi+(C@j)@phi)
        j=np.array([.6,.13,-.08]);step=1e-3
        e=np.eye(3)
        for i,k in ((0,1),(0,2),(1,2)):
            mixed=(energy(j+step*e[i]+step*e[k])
                   -energy(j+step*e[i]-step*e[k])
                   -energy(j-step*e[i]+step*e[k])
                   +energy(j-step*e[i]-step*e[k]))/(4*step**2)
            fd_error=max(fd_error,abs(mixed-response[i,k]))
        phi=-np.linalg.solve(A,C@j)
        stationarity_error=max(stationarity_error,float(np.linalg.norm(A@phi+C@j)))
        assert norm(response-response.T)<1e-14
        assert np.linalg.norm(response@null)<1e-14
        assert np.linalg.eigvalsh(response).max()<1e-14
        assert np.linalg.matrix_rank(response,tol=1e-12)==2
        # Best zero-momentum matched light-pole-only kernel.
        one_pole=np.linalg.inv(K)[0,1]*mass2[0]/(q*q+mass2[0])
        rel_error=abs(one_pole-D[0,1])/abs(D[0,1])
        assert abs(rel_error-q*q/mass2[1])<1e-12
        assert abs(response[0,1]/cw-response[0,2]/cf)<1e-14
        # Omitting reverse response cannot be a Hessian of one scalar energy.
        bad=response.copy();bad[1,0]=0.
        rows.append(dict(spacelike_Q=q,mixed_propagator=float(D[0,1]),
            record_W_exchange_coefficient=float(response[0,1]),
            record_fermion_exchange_coefficient=float(response[0,2]),
            common_coefficient_after_vertex_normalization=float(response[0,1]/cw),
            actual_protocol_source_weight_difference=float(delta_b),
            protocol_weighted_W_response_difference=float(delta_b*response[0,1]),
            protocol_weighted_fermion_response_difference=float(delta_b*response[0,2]),
            response_hessian=response.tolist(),
            zero_matched_light_pole=float(one_pole),
            one_pole_relative_error=float(rel_error),
            omitted_reverse_source_antisymmetric_norm=norm(bad-bad.T)))
    assert inv_error<1e-13 and fd_error<1e-9 and stationarity_error<1e-13
    Kzero=np.diag([mu**2,2*lam*v*v])
    assert np.linalg.inv(Kzero)[0,1]==0
    # Couplings also follow by differentiating the same input mass functions.
    def masses(eta):return (g2*(v+eta)/2)**2,yf*(v+eta)/math.sqrt(2)
    eps=1e-4
    dp=np.array(masses(eps));dm=np.array(masses(-eps))
    derivatives=(dp-dm)/(2*eps)
    vertex_error=float(np.max(abs(derivatives-[cw,cf])))
    assert vertex_error<1e-11
    q1=next(r for r in rows if r['spacelike_Q']==1.)
    assert q1['record_W_exchange_coefficient']>.02
    assert q1['one_pole_relative_error']>.23
    assert abs(sum(O[0,a]*O[1,a] for a in range(2)))<1e-14
    files=[Path(__file__),STAGE/'946/portal_common_process_results.json',
           STAGE/'947/protocol_field_transport_results.json']
    return dict(round=948,date='2026-10-07',all_scientific_checks_passed=True,
        parameters=dict(inherited_K=K.tolist(),scalar_mass_squared=mass2.tolist(),
                        g_source=gs,v=v,g2=g2,y_fermion=yf,m_W=mw,m_fermion=mf,
                        W_mass_squared_derivative=cw,fermion_mass_derivative=cf),
        checks=dict(direct_spectral_kernel_error=inv_error,
                    independent_eliminated_action_mixed_derivative_error=fd_error,
                    on_shell_source_stationarity_error=stationarity_error,
                    independent_mass_vertex_derivative_error=vertex_error,
                    orthogonal_mode_mixed_residue_sum=float(np.dot(O[0],O[1])),
                    zero_portal_mixed_response=0.,
                    all_joint_response_hessians_rank_two=True),
        rows=rows,
        scope=dict(same_946_two_mode_kernel=True,
                   same_947_first_event_source_operator=True,
                   standard_electroweak_and_Yukawa_vertices_are_physical_inputs=True,
                   neutral_record_to_gauge_and_fermion_tree_response=True,
                   reverse_response_from_same_eliminated_action=True,
                   added_gauge_charges_to_record_flavors=False,
                   full_SM_Einstein_finite_time_matching_certified=False,
                   new_joint_detector_probability_bound_certified=False,
                   empirical_parameter_fit=False,
                   microscopic_continuity_or_minimum_scale_inferred=False,
                   full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert abs(a-b)<1e-10,(a,b)
    else:assert a==b,(a,b)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','rows')},
                     ensure_ascii=False,indent=2))
    print(json.dumps(next(r for r in out['rows'] if r['spacelike_Q']==1.),ensure_ascii=False,indent=2))
