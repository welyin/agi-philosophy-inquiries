"""859: existing reference probe, full density variation and receiver cones.
Exact rational principal symbols plus the inherited 64-component coefficient.
This is a candidate-specific classical quadratic check, not quantum transport.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,importlib.util,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
TARGET=HERE/'probe_reference_cone_compatibility_results.json'
import shared_probe_reference_probe as background

def mat(rows):return np.array([[F(x) for x in row] for row in rows],dtype=object)
I=mat(np.eye(4,dtype=int));g=mat(np.diag([-1,1,1,1]))

def symbol(C,a):
    H=a*C;mixed=g@H;trace=sum(mixed[i,i] for i in range(4))
    L=(1-trace/2)*I+mixed/2
    G=L@g@L.T
    assert np.array_equal(L@g,g@L.T)
    assert np.array_equal(G@g,L@L)
    return L,G

def rational_checks():
    longitudinal=mat([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,2]])
    shear=mat([[0,0,0,0],[0,0,0,0],[0,0,0,2],[0,0,2,0]])
    conformal=2*g
    tilt=mat([[0,1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]])
    rows=[];k=np.array([F(-5),F(0),F(3),F(4)],dtype=object)
    assert k@g@k==0
    for name,C in [('old_847_longitudinal',longitudinal),('old_847_shear',shear),('conformal_control',conformal),('time_space_diagnostic',tilt)]:
        for a in (F(-1,10),F(1,20),F(1,10)):
            L,G=symbol(C,a)
            omega=-G[0,0];same=np.array_equal(G,omega*g)
            scalar=np.array_equal(g@C,sum((g@C)[i,i] for i in range(4))*I/4)
            assert same==scalar
            nullvalue=k@G@k
            if name=='old_847_shear':assert nullvalue==48*a+25*a*a
            if name=='old_847_longitudinal':assert nullvalue==16*(2*a-a*a)
            if same:assert nullvalue==0
            rows.append(dict(case=name,amplitude=str(a),effective_inverse_metric=[[str(x) for x in row] for row in G],old_null_covector=[-5,0,3,4],new_quadratic_value=str(nullvalue),same_cone=same))
    return rows

def original_coefficient_checks():
    path=HERE.parent/'802/original_bff_vertex.py'
    spec=importlib.util.spec_from_file_location('round802_density',path)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    candidates=[('longitudinal',np.diag([0.,0.,0.,2.])),('shear',np.array([[0.,0.,0.,0.],[0.,0.,0.,0.],[0.,0.,0.,2.],[0.,0.,2.,0.]])),('conformal',2*old.ETA)]
    rows=[]
    for name,C in candidates:
        a=.1;nu,dc,_=old.frame_tangent(np.eye(4),C,np.zeros((4,4,4)))
        exact=np.array([old.ALPHA[mu]-a*(nu*old.ALPHA[mu]+dc[mu]) for mu in range(4)])
        L=(1-a*np.trace(old.ETA@C)/2)*np.eye(4)+a*old.ETA@C/2
        predicted=np.einsum('ij,jab->iab',L,np.array(old.ALPHA))
        err=float(np.max(abs(exact-predicted)));assert err<1e-13
        differences=[]
        for step in (2e-4,1e-4,5e-5):
            vp,cp,_=old.frame_data(np.eye(4),C,np.zeros((4,4,4)),step)
            vm,cm,_=old.frame_data(np.eye(4),C,np.zeros((4,4,4)),-step)
            measured=np.array(old.ALPHA)-a*(vp*cp-vm*cm)/(2*step)
            differences.append(float(np.max(abs(measured-exact))))
        assert differences[-1]<3e-9
        eigM,U=np.linalg.eigh(exact[0]);assert min(eigM)>.69
        invroot=(U*eigM**-.5)@U.conj().T
        speeds=np.linalg.eigvalsh(invroot@exact[3]@invroot)
        predicted_speed={'longitudinal':1/(1-a),'shear':np.sqrt(1+a*a),'conformal':1.}[name]
        assert np.max(abs(abs(speeds)-predicted_speed))<2e-13
        rows.append(dict(case=name,original_components=64,density_variation_error=err,independent_frame_difference_errors=differences,positive_time_coefficient_min=float(min(eigM)),old_fermion_z_characteristic_speed=float(max(speeds)),unchanged_receiver_z_characteristic_speed=1.))
    return rows

def run():
    bg=background.run();assert bg==json.loads(background.TARGET.read_text('utf-8'))
    exact=rational_checks();coeff=original_coefficient_checks()
    return dict(round=859,all_checks_passed=True,new_numbered_scientific_groups=1,
        background_probe_reproduced=True,background_rows=bg['rows'],exact_principal_cases=exact,original_density_checks=coeff,
        scope='One inherited scalar and the color magnetic scalar form local references on newly constrained classical backgrounds. In a small-coupling active nonconformal stress window, reference rank, unchanged inherited read action and exact original/receiver common cones cannot all hold.',
        original_802_density_and_847_shear_used=True,new_species_added=False,
        exact_common_cone_an_explicit_extra_contract=True,nonzero_probe_changes_free_fermion_quadratic_part=True,
        failure_of_all_relational_references=False,full_reduced_canonical_dictionary_completed=False,
        new_quantum_state_transport_completed=False,original_graph_continuum_limit_proved=False,full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in r.items() if k not in ('background_rows','exact_principal_cases')},ensure_ascii=False,indent=2))
