"""743: original scalar fluctuations, physical Gauss preparation and CAR jets.

Analytic evolution is the full original fixed-graph Hamiltonian. Numerical
checks evaluate the initial jets and normal-packet moments only, not that
Hamiltonian's full time evolution or a continuum interacting theory.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'round743_drafts'))
import quantum_scalar_record_entry as entry
import joint_record_mass_feedback as old
TARGET=HERE/'joint_native_quantum_record_results.json'

def packet(order):
    # Exactly the R=0 normal Gauss packet from717; no width optimization.
    u,wu=np.polynomial.legendre.leggauss(order)
    radial=.7*(u+1)/2;wr=.7*wu/2
    r,s=np.meshgrid(radial,u,indexing='ij')
    z=np.sqrt(1+(r*r+s*s)/6)
    weight=wr[:,None]*wu[None,:]*r**3/z
    weight*=np.exp(-2/(1-(r/.7)**2)-2/(1-s*s))
    probability=weight/weight.sum()
    moment=lambda A:float(np.sum(probability*A))
    variance=moment(s*s)-moment(s)**2
    F=2/z**2;phi_s=np.sqrt(F)*s
    return dict(normalization=float(weight.sum()),mean_x5=moment(s),variance_x5=variance,
        inherited_sin_readout_force_weight=moment(np.sqrt(F)*np.cos(phi_s)),
        nonGaussian_occupation_Wick_coefficient=abs(old.matter.Y['s'])**2*variance,
        connected_scalar_pair_derivative_abs=abs(old.matter.Y['s'])*variance,
        global_fermion_linear_entropy_coefficient_per_independent_node=2*abs(old.matter.Y['s'])**2*variance)

def gauss_jet_check():
    fine,coarse=packet(80),packet(48)
    inherited=old.packet(0.,80)
    assert abs(fine['normalization']-inherited['Z'])<1e-16
    assert abs(fine['mean_x5'])<1e-15
    assert fine['variance_x5']>0 and fine['inherited_sin_readout_force_weight']>0
    qerr=max(abs(fine[k]-coarse[k]) for k in fine)
    assert qerr<1e-7
    # Actual full64-mode two-node mass and nontrivial hopping on CAR vacuum.
    rng=np.random.default_rng(743)
    pair0=(1<<30)|(1<<31);pair1=(1<<62)|(1<<63)
    residuals=[];normal_jet=[];anomalous_jet=[]
    for _ in range(9):
        x=rng.normal(size=(2,5))*.4
        h=np.zeros((64,64),complex);d=np.zeros_like(h)
        for site in range(2):
            sl=slice(site*32,(site+1)*32);h[sl,sl],d[sl,sl]=old.mass_x(x[site])
        C=old.matter.gauge.group_exp(rng.normal(size=8)*.2,3)
        W=old.matter.gauge.group_exp(rng.normal(size=3)*.2,2)
        R=old.matter.representation(C,W,np.exp(.17j))
        hopping=.23-.11j
        h[:32,32:]=hopping*R;h[32:,:32]=hopping.conjugate()*R.conj().T
        image=old.car.quadratic({0:1.},h,d)
        expected={pair0:old.matter.Y['s']*x[0,4],pair1:old.matter.Y['s']*x[1,4]}
        residuals.append(float(old.car.difference(image,expected)))
        # Both actual occupations and their product have the same t^2 jet.
        norm_a=sum(abs(z)**2 for key,z in image.items() if (key>>30)&1)
        norm_ab=sum(abs(z)**2 for key,z in image.items() if (key&pair0)==pair0)
        target=abs(old.matter.Y['s']*x[0,4])**2
        normal_jet.append(float(max(abs(norm_a-target),abs(norm_ab-target))))
        twice=old.car.old.annihilate(old.car.old.annihilate(image,31),30)
        # a*b acting on a^dagger*b^dagger|0> gives -|0>.
        anomalous_jet.append(float(abs(-1j*twice.get(0,0)-1j*old.matter.Y['s']*x[0,4])))
    assert max(residuals+normal_jet+anomalous_jet)<2e-13
    return dict(original_two_node_CAR_modes=64,full_Yukawa_Majorana_and_nontrivial_hopping_retained=True,
        vacuum_first_jet_error=max(residuals),occupation_second_jet_error=max(normal_jet),
        anomalous_first_jet_error=max(anomalous_jet),normal_Gauss_packet=fine,
        quadrature_48_to_80_max_change=qerr,
        old_717_packet_reused=True,
        full_Hb_remains_in_analytic_proof=True,
        no_full_boson_time_evolution_simulated=True,
        alternate_finite_graph_preparation_not_round730_continuum_reference=True)

def run():
    first=entry.run();assert first==json.loads(entry.TARGET.read_text('utf8'))
    second=gauss_jet_check()
    deps=('research_note_598.md','research_note_623.md','research_note_717.md','research_note_719.md',
          'research_note_735.md','research_note_742.md','joint_record_mass_feedback.py',
          'round743_drafts/quantum_scalar_record_entry.py','round743_drafts/quantum_scalar_record_entry_results.json')
    return dict(round=743,tests_run=2,failures=0,errors=0,
        checks=['original_conditional_mass_fluctuation_jet','full_original_Gauss_vacuum_pair_jet'],
        results=dict(conditional_reference=first,normal_Gauss_preparation=second),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original quantum scalar fluctuations generate a nonGaussian fermion four-point cumulant in normal strict-Gauss finite-graph preparations under the full original H. Same Y and scalar variance fix scalar-pair correlation and initial entanglement. The717 readout/source bridge is reused. This is not exact realization of633 occupation instrument, autonomous final readout, the730 continuum reference, a continuum interacting SM or dynamic quantum gravity.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r['results']['normal_Gauss_preparation'],ensure_ascii=False,indent=2))
