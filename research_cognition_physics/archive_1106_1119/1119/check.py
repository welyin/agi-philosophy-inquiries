"""1119: affine finite-context and port-fibre compatibility certificate.

Default: read-only recalculation against saved results. No old device execution.
--write-results: save only this round's results.json.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCES=(
 "archive_1106_/1106/proof.md", "archive_1106_/1106/results.json",
 "archive_1106_/1111/proof.md", "archive_1106_/1111/results.json",
 "archive_1106_/1112/proof.md", "archive_1106_/1114/proof.md",
 "archive_1106_/1115/proof.md", "archive_1106_/1115/results.json",
)
ETA=np.diag([1.,-1.,-1.,-1.])
TOL=2e-10

def norm(a): return float(np.linalg.norm(a,2))
def ok(x):
    assert x<TOL,x
    return float(x)
def expm_h(h, t):
    e,u=np.linalg.eigh(h)
    return (u*np.exp(-1j*t*e))@u.conj().T
def boost(v):
    v=np.array(v,float); b2=v@v; assert b2<1
    g=1/np.sqrt(1-b2); a=np.eye(4)
    a[0,0]=g; a[0,1:]=g*v; a[1:,0]=g*v
    if b2: a[1:,1:]+=(g-1)*np.outer(v,v)/b2
    ok(np.linalg.norm(a.T@ETA@a-ETA))
    return a
def clock():
    # Exactly the undressed H_c block in 1106, not a replacement clock model.
    h=1.5*np.eye(4)
    for j in range(3):
        h[j,j+1]=h[j+1,j]=.5*np.sqrt((j+1)*(3-j))
    m0=15+2*np.sqrt(19)
    return h,m0*np.eye(4)+h
def blocks(items):
    n=sum(a.shape[0] for a in items)
    out=np.zeros((n,n),complex); k=0
    for a in items:
        d=len(a);out[k:k+d,k:k+d]=a;k+=d
    return out
def affine(g,x): return g[0]@x+g[1]
def compose(g2,g1): return (g2[0]@g1[0],g2[0]@g1[1]+g2[1])

def group1(old_clock,old_scatter):
    events2=old_clock["group1_material_and_two_port_events"]["events"][0]["events_o_a_b_f"]
    events=[np.array([x[0],x[1],0.,0.]) for x in events2]
    lam=boost([.3,.2,.1]); shift=np.array([.7,-.4,.25,.15])
    moved=[lam@x+shift for x in events]
    cone=[]
    for i in range(len(events)-1):
        d=events[i+1]-events[i]; dd=moved[i+1]-moved[i]
        assert d[0]>0 and dd[0]>0 and d@ETA@d>0
        cone.append({
          "original_interval":float(d@ETA@d),
          "moved_interval":float(dd@ETA@dd),
          "interval_error":ok(abs(d@ETA@d-dd@ETA@dd)),
        })
    arrival=events[1]; complete=events[-1]
    ok(np.linalg.norm(complete-np.array(old_scatter["parameters"]["scatter_complete"])))
    oldgeo=old_scatter["group4_material_recoil_and_actual_geometry"]
    hits=[]
    for row in oldgeo["incoming_pair"][:2]:
        initial=np.array(row["initial_event"]); final=np.array(row["arrival_event"])
        p=np.array(row["P"]); dt=final[0]-initial[0]
        original_end=initial+dt*p/p[0]
        moved_end=lam@initial+shift+dt*(lam@p)/p[0]
        hits.append(ok(np.linalg.norm(moved_end-(lam@final+shift))))
        ok(np.linalg.norm(original_end-final))
    outgoing=[]
    ef=np.array(old_scatter["parameters"]["scatter_complete"])
    for row in oldgeo["rows"]:
        p=np.array(row["P"]); endpoint=np.array(row["endpoint"])
        sigma=float(row["proper_flight"]); m=old_scatter["parameters"]["m"]
        outgoing.append(ok(np.linalg.norm(
          lam@endpoint+shift-(lam@ef+shift+sigma*(lam@p)/m))))
    return {
      "status":"PASS", "actual_source_clock_events": [x.tolist() for x in events],
      "translated_and_boosted_same_events":[x.tolist() for x in moved],
      "timelike_segments":cone,
      "old_two_incoming_bodies_same_actual_collision_residual":max(hits),
      "all_old_output_endpoint_transport_residual":max(outgoing),
      "all_resources_and_event_cuts_transport_required":True,
      "affine_predeployment_is_new_declared_permission":True,
      "sample_not_proof_of_all_affine_permission":True,
    }

def group2():
    hc,h0=clock()
    v=.6; gamma=1.25; tau=np.pi/2; dt=gamma*tau
    energy=gamma*h0; momentum=gamma*v*h0
    xin=[0.,3.]; xout=[x+v*dt for x in xin]
    jin=blocks([expm_h(momentum,x) for x in xin])
    # Output catalog ordering is reversed: a real port permutation is necessary.
    jout=blocks([expm_h(momentum,x) for x in reversed(xout)])
    pi=np.kron(np.array([[0.,1.],[1.,0.]]),np.eye(4))
    u_tau=blocks([expm_h(h0,tau)]*2)
    u_t=blocks([expm_h(energy,dt)]*2)
    lhs=jout@pi@u_tau
    rhs=pi@u_t@jin
    psi=np.zeros(8,complex);psi[0]=1/np.sqrt(2);psi[4]=1j/np.sqrt(2)
    effect=blocks([np.diag([0,0,0,1])]*2)
    original=pi@u_tau@psi
    transformed=jout@original
    moved_effect=jout@effect@jout.conj().T
    prob=float(np.vdot(original,effect@original).real)
    moved_prob=float(np.vdot(transformed,moved_effect@transformed).real)
    inconsistent_prob=float(np.vdot(transformed,effect@transformed).real)
    assert abs(prob-inconsistent_prob)>.01
    return {
      "status":"PASS", "old_complete_clock_mass_spectral_block":np.linalg.eigvalsh(h0).tolist(),
      "input_port_catalog":xin, "output_catalog_reversed":list(reversed(xout)),
      "full_port_fibre_intertwining_residual":ok(norm(lhs-rhs)),
      "complete_code_unitarity_residual":ok(norm(jin.conj().T@jin-np.eye(8))),
      "correct_transported_instrument_probability":moved_prob,
      "old_probability":prob,
      "instrument_transport_error":ok(abs(prob-moved_prob)),
      "wrong_unchanged_effect_probability":inconsistent_prob,
      "omitted_port_permutation_operator_defect":norm(jout@u_tau-u_t@jin),
      "arbitrary_input_reference_claim_uses_operator_identity":True,
      "J_is_representation_not_physical_free_gate":True,
    }

def group3():
    g1=(boost([.25,-.1,.05]),np.array([.4,.3,-.2,.1]))
    angle=.42
    rot=np.eye(4);rot[1:3,1:3]=[[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]]
    g2=(rot@boost([-.15,.2,.1]),np.array([-.2,.1,.4,-.3]))
    g21=compose(g2,g1)
    point=np.array([1.7,.2,-.4,.15])
    seq=affine(g2,affine(g1,point));direct=affine(g21,point)
    inv=(np.linalg.inv(g21[0]),-np.linalg.inv(g21[0])@g21[1])
    _,h0=clock()
    u=boost([.6,0,0])[:,0]
    p=[u[mu]*h0 for mu in range(4)]
    a=np.array([.8,.2,-.1,.3]);b=np.array([-.3,.5,.4,.2])
    def phase(p,x): return expm_h(sum(ETA[mu,mu]*x[mu]*p[mu] for mu in range(4)),1)
    lam=g21[0]
    pp=[sum(lam[mu,nu]*p[nu] for nu in range(4)) for mu in range(4)]
    return {
      "status":"PASS",
      "affine_composition_error":ok(np.linalg.norm(seq-direct)),
      "affine_inverse_error":ok(np.linalg.norm(affine(inv,direct)-point)),
      "full_phase_translation_composition_error":ok(norm(phase(p,a+b)-phase(p,a)@phase(p,b))),
      "four_P_and_event_phase_covariance_error":ok(norm(phase(pp,lam@a)-phase(p,a))),
      "identity_lifts_are_fixed_in_transported_internal_bases":True,
      "negative_parameter_inverse_not_backward_signal_permission":True,
    }

def group4():
    hc,h0=clock()
    tau=np.pi/2;gamma=1.25
    psi=np.eye(4,dtype=complex)[:,0]
    actual=expm_h(h0,tau)@psi
    wrong=expm_h(gamma*h0,tau)@psi
    p=float(abs(actual[3])**2);pw=float(abs(wrong[3])**2)
    expected=np.sin(tau/2)**6
    wrong_expected=np.sin(gamma*tau/2)**6
    assert abs(pw-p)>.1
    u=boost([.6,0,0])[:,0]
    pvec=[u[mu]*h0 for mu in range(4)]
    rest_energy=sum(ETA[mu,mu]*u[mu]*pvec[mu] for mu in range(4))
    return {
      "status":"PASS","proper_window":tau,"coordinate_duration":gamma*tau,
      "old_Mr_probability_via_1106_head_identity":p,
      "wrong_delta_times_P0_probability":pw,
      "finite_record_probability_gap":abs(pw-p),
      "same_old_record_formula_error":ok(abs(p-expected)),
      "wrong_clock_formula_check":ok(abs(pw-wrong_expected)),
      "u_dot_P_equals_complete_H0_residual":ok(norm(rest_energy-h0)),
      "readout_identity_from_frozen_1106_not_new_detector":True,
    }

def group5(old_scatter,old_parent):
    a=old_scatter["group1_same_complete_device"]
    assert a["full_product_dimension"]==1152
    assert a["cross_sector_actual_geometry_not_certified"]
    for key in ["full_product_energy_commutator","full_product_spatial_commutator",
                "dressed_unitary_residual","complete_H_commutator_residual"]:
        assert a[key]<TOL
    for key,value in old_parent["source_sha256"].items():
        target=(ROOT/"archive_1106_"/key).resolve()
        assert hashlib.sha256(target.read_bytes()).hexdigest()==value,(key,"changed old source")
    # Small full two-body mathematical block, solely to check common-position
    # phase commutation. The actual dressed device stays bound by its source.
    px=np.diag([0.,4.,-4.,0.])
    py=np.diag([0.,-4.,4.,0.])
    energy=10*np.eye(4)
    swap=np.eye(4);swap[0,0]=swap[3,3]=0;swap[0,3]=swap[3,0]=1
    w=swap@expm_h(energy,np.pi/8)
    j=expm_h(.37*px-.22*py,1)
    assert np.linalg.norm(swap-np.eye(4))>1
    q=(energy,px,py,np.zeros((4,4)))
    residual=max(norm(w@p-p@w) for p in q)
    return {
      "status":"PASS","old_full_mother_dimension_bound":1152,
      "old_full_mother_energy_commutator":a["full_product_energy_commutator"],
      "old_full_mother_spatial_commutator":a["full_product_spatial_commutator"],
      "small_two_body_complete_boundary_commutator":ok(residual),
      "common_position_J_boundary_commutator":ok(norm(j@w-w@j)),
      "old_parent_frozen_source_bindings_checked":len(old_parent["source_sha256"]),
      "no_old_288_or_1152_device_execution":True,
      "collision_not_reidentified_as_free_wait":True,
      "crossed_sector_preparation_not_added":True,
    }

def calculate():
    hashes={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in SOURCES}
    c=json.loads((ROOT/SOURCES[1]).read_text(encoding="utf-8"))
    s=json.loads((ROOT/SOURCES[3]).read_text(encoding="utf-8"))
    p=json.loads((ROOT/SOURCES[7]).read_text(encoding="utf-8"))
    return {
      "round":1119, "status":"PASS_CONDITIONAL_AFFINE_FINITE_CONTEXT_CERTIFICATE",
      "new_adopted_axioms":0, "scientific_calibration":0, "numerical_check_groups":5,
      "source_sha256":hashes,
      "scope":{
        "same_1115_complete_parent_and_material_context":True,
        "whole_context_affine_predeployment_publicly_added_condition":True,
        "not_unknown_live_clock_boost_permission":True,
        "strong_continuity_after_fixed_type_identification_only":True,
        "not_full_position_Hilbert_generator_or_CCR":True,
        "J_and_port_permutation_both_retained":True,
        "effective_collision_boundary_not_microscopic_wait":True,
        "all_finite_context_resources_and_original_event_windows_retained":True,
        "no_axiom_adoption_or_goal_completion_claimed":True,
      },
      "groups":{
        "1_old_actual_events_intersections_and_future_cones":group1(c,s),
        "2_complete_clock_port_phase_and_permutation":group2(),
        "3_affine_and_phase_finite_composition":group3(),
        "4_wrong_proper_time_P0_record_countercheck":group4(),
        "5_complete_boundary_and_frozen_parent_evidence":group5(s,p),
      }
    }

def compare(a,b,path="$"):
    if isinstance(a,dict):
        assert isinstance(b,dict) and set(a)==set(b),path
        for k in a: compare(a[k],b[k],path+"."+k)
    elif isinstance(a,(tuple,list)):
        assert isinstance(b,(tuple,list)) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,f"{path}[{i}]")
    elif isinstance(a,(bool,str)) or a is None: assert type(a)is type(b) and a==b,(path,a,b)
    elif isinstance(a,(int,float)):
        assert isinstance(b,(int,float)) and not isinstance(b,bool),path
        assert np.isfinite(a) and np.isfinite(b),path
        assert abs(a-b)<2e-10*(1+max(abs(a),abs(b))),(path,a,b)
    else:raise TypeError(path)

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--write-results",action="store_true");args=parser.parse_args()
    r=calculate();path=HERE/"results.json"
    if args.write_results:
        path.write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n",encoding="utf-8");mode="saved"
    else:
        compare(r,json.loads(path.read_text(encoding="utf-8")));mode="saved results matched"
    for name,item in r["groups"].items():print("PASS",name,item["status"])
    print("PASS 5 groups;",mode,"; old full devices not rerun; no axiom adopted.")

if __name__=="__main__":main()