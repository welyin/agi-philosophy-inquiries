"""852 working calibration of a dynamical switch's missing force.

Classical clock and pointer plus a qubit Bloch mean: this is an Ehrenfest
energy-bookkeeping diagnostic, NOT the proposed covariant quantum theory.
No scientific round is completed by this script.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,json
HERE=Path(__file__).resolve().parent;TARGET=HERE/'relational_probe_budget_probe_results.json'
def run():
    rows=[]
    for q,p,qdot,pdot,sx,sy,sz in [
        (F(1,3),F(2,5),F(3,7),F(-1,4),F(1,5),F(2,5),F(1,3)),
        (F(-2,3),F(3,8),F(-4,9),F(2,7),F(-1,4),F(1,3),F(1,2)),
        (F(4,5),F(-2,7),F(1,2),F(-2,9),F(1,3),F(-1,5),F(-2,5))]:
        assert sx*sx+sy*sy+sz*sz<1
        lam=F(2,7);omega=F(3,5);mass2=F(7,11);f=q*q;df=2*q
        clock_force=-lam*df*p*sz
        pointer_force=-mass2*p-lam*f*sz
        dsx=-2*lam*f*p*sy;dsy=2*lam*f*p*sx-omega*sz;dsz=omega*sy
        clock=qdot*clock_force
        pointer=pdot*(pointer_force+mass2*p)
        qubit=omega*dsx/2
        interaction=lam*(df*qdot*p*sz+f*pdot*sz+f*p*dsz)
        balance=clock+pointer+qubit+interaction
        omitted=pointer+qubit+interaction
        assert balance==0 and omitted==-clock and omitted!=0
        assert sx*dsx+sy*dsy+sz*dsz==0
        rows.append(dict(clock_force=str(clock_force),clock_energy_rate=str(clock),pointer_energy_rate=str(pointer),qubit_energy_rate=str(qubit),interaction_energy_rate=str(interaction),total_energy_rate=str(balance),frozen_switch_missing_energy_rate=str(omitted)))
    return dict(working_round=852,formal_round_completed=False,new_scientific_test_groups=0,
        all_diagnostic_identities_passed=True,rows=rows,
        clock_profile='q squared in a local mechanical diagnostic, not a compact spacetime switch',
        original_spacetime_action_or_quantum_probe_realized=False,
        covariance_QME_state_source_joint_extension_proven=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
