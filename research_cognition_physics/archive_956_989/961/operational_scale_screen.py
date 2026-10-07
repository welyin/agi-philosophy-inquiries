"""961: small mechanism screen, not an expanding-universe simulation.

Five frozen transfer/clock contracts. Prescribed sequential exchange pulses
are a kinematic diagnostic; their autonomous controls and growth costs are
not implemented or certified. Unknown two-dimensional payloads are retained.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/"operational_scale_screen_results.json"

def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))
def evolution(h,t):
    e,v=np.linalg.eigh(h);return (v*np.exp(-1j*t*e))@v.conj().T

def echo(n,j):
    # Same five cursor positions for every case: no free Hilbert-space growth.
    dim=10; U=np.eye(dim,dtype=complex); step=math.pi/(2*j)
    for k in list(range(n))+list(reversed(range(n))):
        h=np.zeros((dim,dim))
        for q in range(2):
            a=2*k+q;b=2*(k+1)+q;h[a,b]=h[b,a]=j
        U=evolution(h,step)@U
    W=np.eye(dim)[:,:2]
    return norm(U@W-(-1)**n*W),2*n*step

def run():
    specs=[
        ("base",2,.2,1.,4.),
        ("all_rates_times_0p4",2,.08,.4,1.6),
        ("more_relays_same_local_rates",4,.2,1.,4.),
        ("same_relays_slower_transfer",2,.1,1.,4.),
        ("refinement_with_faster_transfer",4,.4,1.,4.)]
    rows=[]
    for name,n,j,omega,gap in specs:
        error,T=echo(n,j)
        ticks=omega*T/(2*math.pi)
        local_crossing_ticks=omega/(4*j)
        assert error<1e-13
        assert abs(ticks-n*omega/(2*j))<1e-13
        rows.append(dict(name=name,active_links=n,
            exchange_rate=j,clock_angular_frequency=omega,local_gap=gap,
            echo_parameter_time=T,echo_clock_ticks=ticks,
            single_link_clock_ticks=local_crossing_ticks,
            gap_to_clock_ratio=gap/omega,
            all_unknown_payload_echo_isometry_error=error))
    b,u,g,s,r=rows
    assert abs(b["echo_clock_ticks"]-u["echo_clock_ticks"])<1e-13
    assert abs(b["gap_to_clock_ratio"]-u["gap_to_clock_ratio"])<1e-13
    assert abs(g["echo_clock_ticks"]/b["echo_clock_ticks"]-2)<1e-13
    assert abs(g["gap_to_clock_ratio"]-b["gap_to_clock_ratio"])<1e-13
    assert abs(g["echo_clock_ticks"]-s["echo_clock_ticks"])<1e-13
    assert abs(g["single_link_clock_ticks"]-s["single_link_clock_ticks"])>1
    assert abs(r["echo_clock_ticks"]-b["echo_clock_ticks"])<1e-13
    # Algebraic source-bookkeeping example only: an internal relation bit Q.
    # Same total H changes relation and payload; freezing it is a new model.
    X=np.array([[0.,1.],[1.,0.]]);Z=np.diag([1.,-1.]);I=np.eye(2)
    Hq=.3*np.kron(X,I);Hm=.7*np.kron(I,Z)
    V=.2*np.kron(Z,X);H=Hq+Hm+V
    current=lambda A:1j*(H@A-A@H)
    currents=[current(A) for A in (Hq,Hm,V)]
    balance=norm(sum(currents))
    assert balance<1e-14 and all(norm(A)>.05 for A in currents)
    # No interpretation of this two-qubit Q as a metric or a scale factor.
    paths=[Path(__file__),HERE/"drafts/mechanism_priority_entry.md",
        HERE/"drafts/deferred_material_gravity_interface.md",
        STAGE/"960/research_round_960_checks.json",
        STAGE/"957/drafts/unified_operation_hypotheses_v0_2.md"]
    return dict(round=961,date="2026-10-07",all_scientific_checks_passed=True,
        kind="mechanism map with a finite falsification screen, not a new physical model",
        frozen_transfer_cases=rows,
        decisions=dict(common_rate_rescaling_is_not_observable_expansion=True,
            more_nodes_alone_do_not_imply_observable_expansion=True,
            equal_echo_can_hide_different_local_transfer_mechanisms=True,
            clock_calibrated_delay_and_local_material_ratio_both_required=True),
        internal_relation_source_example=dict(total_current_balance_error=balance,
            subsystem_current_norms=[norm(A) for A in currents],
            new_gravity_or_cosmology_claim=False),
        scope=dict(pulse_schedule_is_given=True,clock_frequency_is_given=True,
            no_autonomous_growth_process_certified=True,
            no_redshift_or_3D_metric_or_Friedmann_equation_derived=True,
            concept_map_is_not_a_mathematical_proof=True,
            source_example_is_not_a_common_physics_realization=True,
            full_goal_completed=False,stop_local_optimization=True),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in paths})

def compare(a,b):
    if isinstance(b,dict):
        assert a.keys()==b.keys()
        for k in b:compare(a[k],b[k])
    elif isinstance(b,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(b,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-13),(a,b)
    else:assert a==b,(a,b)

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");a=p.parse_args()
    if a.write:assert not TARGET.exists()
    result=run()
    if a.write:
        with TARGET.open("x",encoding="utf-8") as out:
            json.dump(result,out,ensure_ascii=False,indent=2);out.write("\n")
    else:compare(read(TARGET),result)
    print(json.dumps({k:v for k,v in result.items() if k!="source_hashes"},
                     ensure_ascii=False,indent=2))
