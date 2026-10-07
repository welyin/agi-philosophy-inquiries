"""993: joint static source matching and conditional thermal-interface stability.
The old material's actual Higgs response and cosmological history are NOT computed.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, importlib.util, json, math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'joint_matching_results.json'

def read(p): return json.loads(p.read_text('utf-8-sig'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def exact(x): return dict(exact=str(x),value=float(x))
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module
def entropy(rho):
    e=np.linalg.eigvalsh(rho)
    assert min(e)>-1e-12
    e=e[e>0]
    return float(-e@np.log(e))
def availability(rho,h,temp):
    e=np.linalg.eigvalsh(h)
    return float(np.trace(rho@h).real-temp*entropy(rho)
                 +temp*np.logaddexp.reduce(-e/temp))
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a: compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=2e-10,abs_tol=2e-12),(a,b)
    else:assert a==b,(a,b)

def run():
    old=read(STAGE/'991/dark_mode_screen_results.json')
    m2=F(old['masses_squared']['higgs']['exact'])
    rows=[]
    # j is an explicit calibration density, NOT a matched native-material source.
    for j in (F(-1,1000),F(1,1000)):
        for data in old['rows']:
            # This is a contrast in dark source relative to the fixed old n=0.
            d=F(data['delta_higgs_source']['exact'])
            cross= -j*d/m2
            def energy(x):return -x*x/(2*m2)
            four_term=energy(j+d)-energy(j)-energy(d)+energy(F(0))
            assert four_term==cross
            z=-(j+d)/m2
            direct=m2*z*z/2+(j+d)*z
            assert direct==energy(j+d)
            assert m2*z+j+d==0
            rows.append(dict(n=data['n'],ordinary_calibration_source=exact(j),
                dark_source_contrast=exact(d),cross_energy_density=exact(cross),
                four_source_contrast=exact(four_term)))
    # Coefficient recorded for the analytic mixed derivative; not an extra test.
    polynomial={(2,0):-1/(2*m2),(1,1):-1/m2,(0,2):-1/(2*m2)}
    ordinary_then_dark=polynomial[(1,1)]
    md2=F(old['masses_squared']['dark']['exact'])
    threshold=m2==4*md2
    assert threshold
    # Free |4> -> |2> matrix element of (a+a^dagger)^2 is sqrt(12).
    a=np.diag(np.sqrt(np.arange(1,8,dtype=float)),1)
    qosc=a+a.T
    pair_element=float((qosc@qosc)[2,4])
    assert abs(pair_element**2-12)<1e-12

    _,h,q,f,_,_=load('native968_for993',STAGE/'968/internal_relay.py').material()
    h=h+.05*f
    e,u=np.linalg.eigh(h)
    beta=2*math.log(2);p=np.exp(-beta*(e-e[0]));p/=sum(p)
    rho=(u*p)@u.T
    thermal=read(STAGE/'992/cooling_availability_results.json')
    source=next(x for x in thermal['rows'] if x['a']==2)
    temp=source['temperature'];A0=availability(rho,h,temp)
    assert abs(A0-source['availability'])<1e-12
    # Prospective budgets, not certified matching errors of the candidate.
    eps=delta=theta=.001
    d=4;span=float(e[-1]-e[0])
    continuity=-delta*math.log(delta)-(1-delta)*math.log1p(-delta)+delta*math.log(d-1)
    bound=2*eps+span*delta+(temp+theta)*continuity+theta*math.log(d)
    lower=A0-bound
    assert lower>0
    witness=[]
    # Old material operators give commuting/noncommuting sensitivity probes.
    # Neither operator is asserted to be its physical Higgs matching derivative.
    for name,op in (('charge_Q',q),('exchange_F',f)):
        op=op/np.linalg.norm(op,2)
        for sign in (-1,1):
            hp=h+sign*eps*op
            target=np.outer(u[:,-1],u[:,-1].conj())
            rp=(1-delta)*rho+delta*target
            tp=temp+sign*theta
            actual_delta=float(np.sum(abs(np.linalg.eigvalsh(rp-rho)))/2)
            Ap=availability(rp,hp,tp)
            assert actual_delta<=delta+1e-12
            assert np.linalg.norm(hp-h,2)<=eps+1e-12
            assert abs(Ap-A0)<=bound+1e-12
            witness.append(dict(operator=name,sign=sign,actual_trace_distance=actual_delta,
                actual_availability=Ap,absolute_change=abs(Ap-A0)))
    inputs=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/selection.md',
        HERE/'common_candidate_v1.md',STAGE/'research_note_952.md',
        STAGE/'981/drafts/common_parent_contract_v1.md',STAGE/'research_note_990.md',
        STAGE/'991/dark_mode_screen_results.json',STAGE/'research_note_991.md',
        STAGE/'992/cooling_availability_results.json',STAGE/'research_note_992.md',
        STAGE/'968/internal_relay.py',STAGE/'956/native_material_interface.py']
    return dict(round=993,all_scientific_checks_passed=True,
        kind='joint_candidate_adoption_with_conditional_interface_calibration',
        static_source_rows=rows,mixed_susceptibility=exact(ordinary_then_dark),
        free_pair_frequency_equals_higgs_mass=threshold,
        free_pair_4_to_2_squared_element=pair_element**2,
        interacting_resonance_claimed=False,
        thermal_interface=dict(a=2,temperature=temp,old_availability=A0,energy_span=span,
            prospective_operator_error=eps,prospective_state_trace_error=delta,
            prospective_temperature_error=theta,error_bound=bound,
            conditional_availability_lower_bound=lower,sensitivity_probes=witness),
        new_physical_species_beyond_D991=0,
        native_higgs_matching_computed=False,prospective_budgets_physically_certified=False,
        cross_term_is_new_general_theorem=False,full_common_process_verified=False,
        cosmological_arrow_derived=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in inputs})

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
