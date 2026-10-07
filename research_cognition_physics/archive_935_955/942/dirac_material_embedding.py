"""942: a declared neutral Dirac multiplet realizes the old internal update.

Reuses 720 FW/source identities and 940 Dirac stress kernels.
New species and engineered mass mixing are physical inputs, not derived matter.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'dirac_material_embedding_results.json'
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
old=load('old925',STAGE/'925/record_resource_exchange.py')
weak=load('old940',STAGE/'940/native_exchange_bridge.py')
I2=np.eye(2);PAULI=weak.PAULI;BETA=weak.G0
ALPHA=[BETA@g for g in weak.GAMMA[1:]]
def norm(a):return float(np.linalg.norm(a,2))
def fun(a,f):
    e,v=np.linalg.eigh(a);return (v*f(e))@v.conj().T
def positive_data(mass,p):
    d=len(mass);eye=np.eye(d);p2=float(p@p)
    E=fun(mass,lambda x:np.sqrt(p2+x*x))
    A=fun(mass,lambda x:np.sqrt((np.sqrt(p2+x*x)+x)/(2*np.sqrt(p2+x*x))))
    B=fun(mass,lambda x:1/np.sqrt(2*np.sqrt(p2+x*x)*(np.sqrt(p2+x*x)+x)))
    sp=sum(p[i]*PAULI[i] for i in range(3))
    W=np.vstack([np.kron(I2,A),np.kron(sp,B)])
    H=sum(p[i]*np.kron(ALPHA[i],eye) for i in range(3))+np.kron(BETA,mass)
    return H,W,np.kron(I2,E)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run():
    z=old.model();h=z['h'];d=z['d'];mass=100*np.eye(d)+h
    p=np.array([12.,-9.,13.]);H,W,E=positive_data(mass,p)
    isometry=norm(W.conj().T@W-np.eye(2*d));intertwiner=norm(H@W-W@E)
    shell=norm(H@H-np.kron(np.eye(4),mass@mass+float(p@p)*np.eye(d)))
    es=np.linalg.eigvalsh(E);all_spectrum=np.linalg.eigvalsh(H)
    spectrum=norm(np.diag(all_spectrum-np.sort(np.r_[-es,es])))
    assert max(isometry,intertwiner,spectrum)<2e-11 and shell<2e-10
    # Actual old unknown-record isometry in the positive rest sector.
    H0,W0,E0=positive_data(mass,np.zeros(3));initial=[];target=[]
    T=np.pi/(2*z['g']);eye=np.eye(d)
    for q in range(3):
        a=z['idx'](0,q,3);b=z['idx'](1,q,3+q)
        initial.append(eye[:,a]);target.append(-1j*np.exp(-1j*T*z['h0'][a,a])*eye[:,b])
    initial=np.kron(np.array([[1.],[0.]]),np.column_stack(initial))
    target=np.kron(np.array([[1.],[0.]]),np.column_stack(target))
    evolved=fun(H0,lambda e:np.exp(-1j*T*e))@W0@initial*np.exp(1j*100*T)
    rest_error=norm(evolved-W0@target)
    # Geometry source is the original local tetrad derivative, not a fitted scalar mass.
    dp=np.array([p[0],-p[1],0.])
    G=sum(dp[i]*np.kron(ALPHA[i],np.eye(d)) for i in range(3))
    expected_geo=np.kron(I2,fun(mass,lambda x:float(p@dp)/np.sqrt(float(p@p)+x*x)))
    geo_error=norm(W.conj().T@G@W-expected_geo)
    lapse_error=norm(W.conj().T@H@W-E)
    # A mass perturbation changes the positive-spinor frame. Reuse the old 720 identity.
    D=z['nref']@z['nref'];field_source=np.kron(BETA,D)
    physical_source=W.conj().T@field_source@W
    delta=2e-4
    def data(t):return positive_data(mass+t*D,p)
    plus2=data(2*delta);plus=data(delta);minus=data(-delta);minus2=data(-2*delta)
    dW=(-plus2[1]+8*plus[1]-8*minus[1]+minus2[1])/(12*delta)
    dE=(-plus2[2]+8*plus[2]-8*minus[2]+minus2[2])/(12*delta)
    connection=W.conj().T@dW
    transported=dE-(E@connection-connection@E)
    source_error=norm(transported-physical_source)
    connection_anti_error=norm(connection+connection.conj().T)
    omitted_frame_gap=norm(dE-physical_source)
    assert rest_error<2e-12 and geo_error<2e-12 and lapse_error<2e-11
    assert source_error<1e-6 and connection_anti_error<1e-7 and omitted_frame_gap>1e-7
    # Exact projectors representing old logical records are generally momentum dependent.
    PB=np.diag(np.r_[np.zeros(d//2),np.ones(d//2)])
    effect=W@np.kron(I2,PB)@W.conj().T
    effect_error=norm(effect@effect-effect)
    Hminus,Wminus,Eminus=positive_data(mass,-p)
    effect_minus=Wminus@np.kron(I2,PB)@Wminus.conj().T
    multiplier_difference=norm(effect-effect_minus)
    assert effect_error<1e-12 and multiplier_difference>.1
    # Conserved on-shell stress of this material, with its actual mass spectrum.
    kernel_error=0.;ward_error=0.;static_masses=np.linalg.eigvalsh(mass)
    for m in static_masses:
        p0=np.array([.7,-.4,-.25]);p1=np.array([.7,-.4,.25])
        q,j,t,_,_=weak.source(p0,p1,float(m));ql=weak.ETA@q
        ward_error=max(ward_error,float(np.max(abs(np.einsum('a,abij->bij',ql,t)))))
        for a,b in ((0,0),(0,1),(1,0),(1,1)):
            tt=t[:,:,a,b];jj=np.zeros(4,complex)
            cov,rad,con,*_=weak.contractions(q,tt,tt,jj,jj)
            scale=max(1.,abs(cov));kernel_error=max(kernel_error,float(abs(cov-rad-con)/scale))
    assert kernel_error<1e-12 and ward_error<1e-10
    inherited=json.loads((STAGE/'929/joint_protocol_selection_results.json').read_text('utf-8'))
    assert inherited['same_six_finite_contracts_jointly_realized']
    paths=[STAGE/'925/record_resource_exchange.py',STAGE/'929/joint_protocol_selection_results.json',
           STAGE/'940/native_exchange_bridge.py',STAGE/'941/composite_operation_bridge_results.json',
           STAGE.parent/'archive_702_741/research_note_720.md',Path(__file__)]
    return dict(round=942,date='2026-10-07',all_scientific_checks_passed=True,
        new_neutral_Dirac_flavours=d,Dirac_symbol_dimension=4*d,positive_one_particle_dimension=2*d,
        mass_minimum=float(static_masses.min()),mass_maximum=float(static_masses.max()),
        positive_isometry_error=isometry,actual_Dirac_intertwiner_error=intertwiner,
        squared_symbol_error=shell,full_positive_and_negative_spectrum_error=spectrum,
        old_unknown_record_rest_isometry_error=rest_error,
        geometry_source_projection_error=geo_error,lapse_source_projection_error=lapse_error,
        changing_mass_source=dict(transported_source_error=source_error,
            frame_connection_antihermiticity_error=connection_anti_error,
            source_gap_if_frame_omitted=omitted_frame_gap,
            original_720_frame_identity_reused=True),
        logical_record=dict(projector_error=effect_error,
            opposite_momentum_multiplier_difference=multiplier_difference,
            point_local_instrument_or_spatial_role_partition_proved=False),
        actual_material_stress=dict(all_60_masses_four_spin_pairs_checked=True,
            maximum_Ward_error=ward_error,relative_covariant_vs_constraint_kernel_error=kernel_error),
        general_929_embedding=dict(required_neutral_Dirac_flavours=inherited['full_hilbert_dimension'],
            internal_spectrum_maximum=inherited['full_H_maximum'],
            original_abstract_six_protocol_isometry_transports_in_rest_sector=True,
            physical_spatial_agents_or_local_access_automatically_realized=False),
        same_new_action_fixes_update_propagation_and_sources=True,
        neutral_vectorlike_representation_introduces_no_SM_gauge_charges=True,
        original_E_rec_action_unchanged=False,
        full_new_action_finite_gravity_evolution_or_matching_verified=False,
        full_goal_completed=False,
        source_hashes={str(q.relative_to(ROOT)):sha(q) for q in paths})
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');a=ap.parse_args();r=run()
    if a.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(r,f,ensure_ascii=False,indent=2);f.write('\n')
    else:
        oldr=json.loads(TARGET.read_text('utf-8'))
        def compare(x,y):
            if isinstance(x,dict):
                assert x.keys()==y.keys()
                for k in x:compare(x[k],y[k])
            elif isinstance(x,list):
                assert len(x)==len(y)
                for v,w in zip(x,y):compare(v,w)
            elif isinstance(x,float):assert abs(x-y)<1e-8
            else:assert x==y,(x,y)
        compare(r,oldr)
    print(json.dumps({k:v for k,v in r.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
