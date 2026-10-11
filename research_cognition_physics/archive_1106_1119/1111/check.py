"""1111: complete-device controlled recoil and actual finite CO endpoints.

Default is read-only. --write-results writes only this round's results.json.
The effective collision is a declared positive-duration local material law;
it is not a microscopic scattering derivation or full-stage certification.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
RESULTS=HERE/"results.json"
TOL=2e-9

def tensor(*args):
    out=np.array([1.0+0j])
    for arg in args: out=np.kron(out,arg)
    return out

def opnorm(x):
    if min(x.shape)<=8: return float(np.linalg.norm(x,2))
    return float(np.linalg.norm(x,"fro"))

def small(x,label):
    value=float(x)
    if value>TOL: raise AssertionError(f"{label}: {value}")
    return value

def evolve(h,t,a):
    w,v=np.linalg.eigh(h)
    return v@(np.exp(-1j*t*w)[:,None]*(v.conj().T@a))

def mink(p,q): return float(p[0]*q[0]-p[1:]@q[1:])

def boost(beta):
    beta=np.asarray(beta,float); b2=beta@beta
    g=1/math.sqrt(1-b2)
    out=np.eye(4); out[0,0]=g
    out[0,1:]=-g*beta; out[1:,0]=-g*beta
    if b2: out[1:,1:]+=((g-1)/b2)*np.outer(beta,beta)
    return out

def model():
    # The same complete 1106 material, not a frozen bare receiver.
    i2=np.eye(2,dtype=complex); i3=np.eye(3,dtype=complex); i72=np.eye(72)
    z=i2[:,0]; one=i2[:,1]; vac=i3[:,0]
    x=np.array([[0,1],[1,0]],complex)
    y=np.array([[0,-1j],[1j,0]],complex)
    b=-y; ap=(i2+b)/2; am=(i2-b)/2
    ur=np.kron(ap,i2)-1j*np.kron(am,y)
    um=np.kron(np.outer(z,z),i2)-1j*np.kron(np.outer(one,one),x)
    swap=np.zeros((9,9),complex)
    for a in range(3):
        for b0 in range(3): swap[b0*3+a,a*3+b0]=1
    gates=[np.kron(np.eye(8),swap),np.kron(ur,np.eye(18)),
           np.kron(np.kron(i2,um),np.eye(9))]
    hc=1.5*np.eye(4)
    hfk=1.5*np.eye(288,dtype=complex)
    products=[i72]
    for j,g in enumerate(gates):
        weight=math.sqrt((j+1)*(3-j))/2
        hc[j+1,j]=hc[j,j+1]=weight
        hfk[(j+1)*72:(j+2)*72,j*72:(j+1)*72]=weight*g
        hfk[j*72:(j+1)*72,(j+1)*72:(j+2)*72]=weight*g.conj().T
        products.append(g@products[-1])
    gamma=np.zeros((288,288),complex)
    for j,p in enumerate(products): gamma[j*72:(j+1)*72,j*72:(j+1)*72]=p
    baseline=15+2*math.sqrt(19)
    hout=baseline*np.eye(288)+hfk
    hin=baseline*np.eye(288)+np.kron(hc,i72)
    clock=np.array([1,1j])/math.sqrt(2)
    ready=np.column_stack([tensor(np.eye(4)[:,0],clock,z,z,i3[:,q],vac)
                           for q in (1,2)])
    arrival=23/15; delta=math.pi/4; T=math.pi
    phase=np.exp(-1j*15*arrival)
    wcap=gamma@evolve(hin,delta,phase*ready)
    def w(t):
        assert t>=delta
        return evolve(hout,t-delta,wcap)
    return locals()

def calculate():
    m=model(); hout=m["hout"]; gamma=m["gamma"]; prod=m["products"]; T=m["T"]
    wt=m["w"](T); ds=math.pi/8; mass=3.0; k=4.0; E=5.0; flight=2.0
    x=np.array([[0,1],[1,0]],complex)
    I2=np.eye(2); I144=np.eye(144)
    pc=np.kron(np.eye(8),np.kron(np.diag([0,0,1]),np.eye(3)))
    pr=np.kron(np.eye(24),np.diag([0,0,1]))
    sin=np.kron(np.eye(72)-pc,I2)+np.kron(pc,x)
    sreceiver=np.kron(np.eye(72)-pr,I2)+np.kron(pr,x)
    sd=np.zeros((576,576),complex)
    for j,p in enumerate(prod):
        pj=np.kron(p,I2)
        sd[j*144:(j+1)*144,j*144:(j+1)*144]=pj@sin@pj.conj().T
    ht=np.kron(hout,I2)+2*E*np.eye(576)
    bare=np.kron(np.eye(4),sreceiver)
    comm=sd@ht-ht@sd
    barecomm=bare@ht-ht@bare
    unitres=small(opnorm(sd.conj().T@sd-np.eye(576)),"dressed unitary")
    commres=small(opnorm(comm),"complete energy")
    headres=small(opnorm(sd[432:,432:]-sreceiver),"head3 actual Qr control")
    # Two distinguishable material factors are individually two dimensional.
    # chi=(xx,yy) is a correlated code, not the complete A tensor B mother.
    # Extend the scattering by identity on xy,yx and keep their actual P.
    fullperm=np.arange(1152)
    active=np.array([4*d+c for d in range(288) for c in (0,3)])
    for d0 in range(288):
        head=d0//72; data=d0%72
        control=(data//3)%3 if head==0 else data%3
        if control==2:
            fullperm[4*d0]=4*d0+3
            fullperm[4*d0+3]=4*d0
    assert np.array_equal(np.sort(fullperm),np.arange(1152))
    assert np.array_equal(fullperm[fullperm],np.arange(1152))
    sfull=np.eye(1152,dtype=complex)[fullperm]
    full_sector=small(opnorm(sfull[np.ix_(active,active)]-sd),"same correlated code restriction")
    hfull=np.kron(hout,np.eye(4))+2*E*np.eye(1152)
    full_H=small(opnorm(hfull[fullperm,:]-hfull[:,fullperm]),"full product mother H")
    pairpx=np.array([0,k,-k,0.]); pairpy=-pairpx
    full_spatial=max(small(float(np.max(abs(np.tile(p,288)[fullperm]-np.tile(p,288)))),
                           "full product mother spatial P") for p in (pairpx,pairpy))
    cross=np.array([4*d+c for d in range(288) for c in (1,2)])
    assert np.array_equal(fullperm[cross],cross)
    assert opnorm(barecomm)>0.1
    small(abs(opnorm(barecomm)-8*math.sqrt(3)),"analytic bare-control failure")
    # Apply the declared full positive collision after the same physical clock runs.
    before=np.kron(m["w"](T-ds),np.array([[1.0],[0.0]]))
    free=np.kron(m["w"](T),np.array([[1.0],[0.0]]))*np.exp(-2j*E*ds)
    # Independently evolve 288x4 rather than exponentiating a 576 matrix.
    before3=before.reshape(288,2,2)
    afterfree=np.empty_like(before3)
    for channel in range(2):
        afterfree[:,channel,:]=evolve(hout,ds,before3[:,channel,:])*np.exp(-2j*E*ds)
    actual=sd@afterfree.reshape(576,2)
    actual_identity=small(opnorm(actual-sd@free),"same uninterrupted clock")
    expected=bare@free
    endpoint=small(opnorm(actual-expected),"actual receiver controlled recoil")
    iso=small(opnorm(actual.conj().T@actual-np.eye(2)),"all unknown inputs/reference")
    hstart=before.conj().T@ht@before
    energy=small(opnorm(actual.conj().T@ht@actual-hstart),"complete energy on source code")
    energy2=small(opnorm((ht@actual).conj().T@(ht@actual)-(ht@before).conj().T@(ht@before)),
                  "complete squared energy")
    # The global output retains phase: Qr alone is intentionally not undisturbed.
    plus=np.array([1,1])/math.sqrt(2)
    out=(actual@plus).reshape(96,3,2)
    rho=np.einsum("aqc,ard->qcrd",out,out.conj()).reshape(6,6)
    idealpair=np.zeros(6,complex); idealpair[2]=idealpair[5]=1/math.sqrt(2)
    jointres=small(opnorm(rho-np.outer(idealpair,idealpair.conj())),"full Qr/path coherence")
    reduced=np.trace(rho.reshape(3,2,3,2),axis1=1,axis2=3)
    idealq=np.array([0,1,1],complex)/math.sqrt(2)
    loss=float(np.linalg.norm(np.linalg.eigvalsh(reduced-np.outer(idealq,idealq.conj())),1)/2)
    small(abs(loss-.5),"discarded recoil dephases Qr")
    # Original actual CO ports of the two material participants.
    center=np.array([m["arrival"]+T,1.,0.,0.])
    begin=center-np.array([ds,0,0,0])
    incoming=[]
    for word_slot in (0,1):
        slot_begin=begin+np.array([2*math.pi*word_slot,0,0,0])
        pair_action=0.; chi_in=0.; chi_out=0.
        for sign in (1,-1):
            p0=np.array([E,sign*k,0.,0.])
            start=np.r_[0.,slot_begin[1:]-p0[1:]/p0[0]*slot_begin[0]]
            small(float(np.max(abs(start[1:]+p0[1:]/p0[0]*slot_begin[0]-slot_begin[1:]))),
                  "predeployed pair actually arrives")
            action=mink(p0,slot_begin-start); pair_action+=action
            chi_in+=float(p0[1:]@start[1:]); chi_out+=float(p0[1:]@slot_begin[1:])
            incoming.append({"word_slot":word_slot,"P":p0.tolist(),"initial_event":start.tolist(),
                             "arrival_event":slot_begin.tolist(),"action":action})
        small(abs(pair_action-(2*E*slot_begin[0]+chi_in-chi_out)),
              "same incoming pair momentum and port phase")
    axis=[np.array([1.,0,0]),np.array([0.,1,0])]
    momenta=[]; endpoints=[]
    for a in axis:
        pm=[np.r_[E,k*a],np.r_[E,-k*a]]
        ee=[center+np.r_[flight,flight*p[1:]/p[0]] for p in pm]
        momenta.append(pm); endpoints.append(ee)
    B=boost([.3,.2,.1]); u=np.array([1.,0,0,0])
    rows=[]; max_mass=0.; max_tau=0.; max_path=0.; max_phase=0.; max_total=0.
    for s in range(2):
        total=np.sum(momenta[s],axis=0)
        max_total=max(max_total,float(np.max(abs(total-2*E*u))))
        for j in range(2):
            p=momenta[s][j]; e=endpoints[s][j]; dx=e-center
            pp=B@p; dd=B@dx
            max_mass=max(max_mass,abs(mink(pp,pp)-mass*mass))
            tau=math.sqrt(mink(dx,dx)); taup=math.sqrt(mink(dd,dd))
            max_tau=max(max_tau,abs(tau-taup))
            max_path=max(max_path,float(np.max(abs(dd[1:]-pp[1:]/pp[0]*dd[0]))))
            max_phase=max(max_phase,abs(mink(p,dx)-mink(pp,dd)))
            assert dd[0]>0 and mink(dd,dd)>0
            rows.append({"channel":s,"body":j,"P":p.tolist(),"endpoint":e.tolist(),
                         "proper_flight":tau,"action":mink(p,dx),
                         "boosted_P":pp.tolist(),"boosted_endpoint":(B@e).tolist()})
    small(max_mass,"individual boosted mass")
    small(max_tau,"actual proper flight")
    small(max_path,"same physical endpoint and momentum")
    small(max_phase,"rail phase in same reference")
    small(max_total,"opposite recoils both retained")
    changes=[momenta[1][j]-momenta[0][j] for j in range(2)]
    small(float(np.max(abs(changes[0]+changes[1]))),"recoil cancellation")
    bx=boost([.6,0,0])
    energychanges=[float((bx@p)[0]) for p in changes]
    small(float(np.max(abs(np.asarray(energychanges)-np.array([3.,-3.])))),
          "moving-reference opposite energy transfer")
    # Endpoints after y scattering cannot be kept equal to preselected x rails.
    wrong=float(np.linalg.norm(endpoints[1][0][1:]-endpoints[0][0][1:]))
    small(abs(wrong-1.6*math.sqrt(2)),"finite wrong endpoint margin")
    assert wrong>2
    # At fixed parameter masks the actual participant's membership of a port is observable.
    # Basis order is the channel x/y; no measurement of unknown Q was inserted.
    portx=np.diag([1.,0.]); porty=np.diag([0.,1.])
    fx=actual.conj().T@np.kron(np.eye(288),portx)@actual
    fy=actual.conj().T@np.kron(np.eye(288),porty)@actual
    readres=max(small(opnorm(fx-np.diag([1,0])),"actual x endpoint effect"),
                small(opnorm(fy-np.diag([0,1])),"actual y endpoint effect"))
    # One added pair: exact coherent finite composition in the same complete device.
    # An early/late pair is a finite future material slot, not a cloned unknown input.
    arr=actual.reshape(288,2,2)
    both=np.zeros((288,2,2,2),complex); both[:,:,0,:]=arr
    # A second dressed operation at a later odd-pi endpoint; the full h keeps running.
    wait=2*math.pi
    flat=both.reshape(288,-1)
    afterwait=evolve(hout,wait,flat).reshape(288,2,2,2)*np.exp(-4j*E*wait)
    # Coordinate momentum-basis phase of both pairs; endpoint rail basis phases
    # are separately known from p.dx. The last ds of wait is the second
    # positive local collision window, not an extra instantaneous operation.
    second=np.empty_like(afterwait)
    for old in range(2):
        second[:,old,:,:]=(sd@afterwait[:,old,:,:].reshape(576,2)).reshape(288,2,2)
    second=second.reshape(1152,2)
    twoiso=small(opnorm(second.conj().T@second-np.eye(2)),"two-pair preserved source")
    second_full=np.zeros((4608,2),complex)
    two_active=np.array([16*d+4*a+b for d in range(288) for a in (0,3) for b in (0,3)])
    second_full[two_active,:]=second
    small(opnorm(second_full.conj().T@second_full-np.eye(2)),"two pairs full product output")
    # Derive only correlations, not an unprovided second capture preparation law.
    probs=[]
    for q0 in range(2):
        state=second[:,q0].reshape(288,2,2)
        probs.append([[float(np.sum(abs(state[:,a,b])**2)) for b in range(2)] for a in range(2)])
        target=np.zeros((2,2)); target[q0,q0]=1
        small(float(np.max(abs(np.asarray(probs[-1])-target))),"finite coherent pair history")
    # Global conservation under Lorentz is scalar multiples of the same total H.
    fourp=max(small(opnorm((B@u)[mu]*comm),"boosted full P") for mu in range(4))
    source=HERE.parent/"1106"/"check.py"
    return {
      "status":"PASS_CONDITIONAL_DRESSED_RECOIL_CO_ENDPOINT_BRIDGE",
      "round":1111,"scientific_calibration":0,"new_adopted_axioms":0,"numerical_check_groups":7,
      "scope":{"same_1106_HFK_continues":True,"positive_local_collision_law_declared":True,
               "only_single_original_source_cell":True,"bare_Qr_not_claimed_undisturbed":True,
               "all_resources_and_both_recoils_retained":True,
               "all_SR_or_exact_C_effective_completion_not_claimed":True,
               "finite_word_predeployment_law_declared":True,
               "port_effect_diagnostic_not_new_detector":True},
      "parameters":{"collision_duration":ds,"scatter_complete":center.tolist(),
                    "scatter_start":begin.tolist(),"m":mass,"k":k,"E":E,"flight":flight},
      "source_sha256":{"../1106/check.py":hashlib.sha256(source.read_bytes()).hexdigest()},
      "group1_same_complete_device":{"full_product_dimension":1152,"correlated_sector_dimension":576,
         "full_product_restriction_residual":full_sector,"full_product_energy_commutator":full_H,
         "full_product_spatial_commutator":full_spatial,
         "cross_sector_scattering_identity":True,"cross_sector_actual_geometry_not_certified":True,
         "dressed_unitary_residual":unitres,
         "complete_H_commutator_residual":commres,"head3_receiver_control_residual":headres,
         "bare_receiver_commutator_frobenius":opnorm(barecomm)},
      "group2_positive_window_complete_transfer":{"uninterrupted_clock_residual":actual_identity,
         "same_Qr_controlled_endpoint_residual":endpoint,"full_isometry_residual":iso,
         "whole_energy_residual":energy,"whole_energy_square_residual":energy2},
      "group3_arbitrary_Q_reference_and_recoil_record":{"joint_coherence_residual":jointres,
         "bare_Qr_plus_input_half_trace_error_if_recoil_ignored":loss,
         "port_dictionary_effect_residual":readres},
      "group4_material_recoil_and_actual_geometry":{"rows":rows,"incoming_pair":incoming,
         "body_recoils":[p.tolist() for p in changes],"total_pair_P_residual":max_total,
         "wrong_unchanged_rail_endpoint_distance":wrong},
      "group5_same_full_Lorentz_transport":{"beta":[.3,.2,.1],"mass_residual":max_mass,
         "proper_time_residual":max_tau,"motion_endpoint_residual":max_path,
         "rail_action_residual":max_phase,"full_four_P_intertwiner_residual":fourp,
         "beta_x_06_individual_energy_changes":energychanges,
         "whole_event_slice_transported_not_resliced":True,
         "second_preparation_beta_declared_as_material_input":True},
      "group6_finite_two_pair_word":{"full_product_output_dimension":4608,"active_code_dimension":1152,
         "isometry_residual":twoiso,"joint_material_channel_probabilities":probs,
         "no_new_global_clock_or_controller":True,
         "second_pair_predeployed_from_actual_t0_ports":True,
         "incoming_pair_port_phase_identity_checked":True},
      "group7_deleted_contracts_and_scope":{"bare_control_fails":opnorm(barecomm)>0.1,
         "discard_recoil_fails_complete_Q_identity":loss>.49,
         "keep_old_endpoint_fails":wrong>2,
         "full_influence_covariance_requires_declared_primitive_class":True,
         "fixed_deadline_general_C_compiler_not_proved":True}
    }

def compare(saved,actual,path="root"):
    if isinstance(actual,dict):
        if not isinstance(saved,dict) or saved.keys()!=actual.keys(): raise AssertionError(path+": keys")
        for k in actual: compare(saved[k],actual[k],path+"."+k)
    elif isinstance(actual,list):
        if not isinstance(saved,list) or len(saved)!=len(actual): raise AssertionError(path+": length")
        for i,(s,a) in enumerate(zip(saved,actual)): compare(s,a,f"{path}[{i}]")
    elif isinstance(actual,float):
        if not math.isclose(saved,actual,rel_tol=3e-8,abs_tol=3e-8): raise AssertionError(f"{path}: {saved} != {actual}")
    elif saved!=actual: raise AssertionError(f"{path}: {saved} != {actual}")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    args=parser.parse_args(); result=calculate()
    if args.write_results:
        RESULTS.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
        print("WROTE results.json")
    else:
        compare(json.loads(RESULTS.read_text(encoding="utf-8")),result)
        print("SAVED_RESULTS_MATCH")
    print("PASS: 7 numerical groups; full product mother 1152; two-pair physical isometry 4608 x 2")
    print("Conditional actual recoil / CO endpoint connection; scientific calibration: 0")
    print("Whole-stage SR, exact physical C and microscopic scattering: not claimed")

if __name__=="__main__": main()

