"""A bounded nuclear formation interface; not a BBN history or QCD derivation."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/'nuclear_formation_results.json'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def num(x):return dict(exact=str(x),value=float(x))


def calculate():
    # NIST CODATA 2022 tabulated central values. Fractions certify arithmetic,
    # not experimental exactness; covariance/uncertainty is not propagated.
    mp,mn,md=map(F,('938.27208943','939.56542194','1875.61294500'))
    W=mp+mn;B=W-md
    photon=(W*W-md*md)/(2*W)
    recoil=W-photon-md
    threshold=(W*W-md*md)/(2*md)
    assert B>0 and photon>0 and recoil>0
    assert (md+recoil)**2-photon**2==md**2
    assert photon+md+recoil==W
    assert threshold+md>W and md*md+2*md*threshold==W*W
    # Bad bookkeeping: count an extra positive released B while retaining
    # already reduced deuteron mass and the actual photon/recoil.
    double_count=md+recoil+photon+B-W
    assert double_count==B
    # If photon is incorrectly assigned all of B, its on-shell recoil costs
    # additional energy; rationalize the square root to avoid cancellation.
    bad_recoil=float(B)**2/(math.hypot(float(md),float(B))+float(md))
    eta=6.1e-10;zeta3=1.2020569031595943;a=0.13;b=0.87
    prefactor=eta*2*zeta3/math.pi**2*0.75*(2*math.pi)**1.5
    rows=[]
    for temp in (1.0,0.1,0.08,0.06):
        R=prefactor*(float(md/(mn*mp))*temp)**1.5*math.exp(float(B)/temp)
        # Stable smaller root: y=R(a-y)(b-y), 0<=y<=a.
        middle=1+R*(a+b)
        disc=middle*middle-4*R*R*a*b
        y=2*R*a*b/(middle+math.sqrt(disc))
        yn,yp=a-y,b-y
        assert 0<y<a and abs(yn+yp+2*y-1)<2e-15
        assert math.isclose(y,R*yn*yp,rel_tol=2e-13)
        # Direct Maxwell-Boltzmann densities give the same chemical equilibrium.
        nb=eta*2*zeta3/math.pi**2*temp**3
        # Thermal fugacity z_i=exp((mu_i-M_i)/T).
        zn=nb*yn/(2*(float(mn)*temp/(2*math.pi))**1.5)
        zp=nb*yp/(2*(float(mp)*temp/(2*math.pi))**1.5)
        zd=nb*y/(3*(float(md)*temp/(2*math.pi))**1.5)
        residual=math.log(zd)-math.log(zn)-math.log(zp)-float(B)/temp
        assert max(zn,zp,zd)<1e-8 and abs(residual)<5e-13
        rows.append(dict(temperature_MeV=temp,R=R,Yd=y,Yn=yn,Yp=yp,
                         neutron_fraction_bound=y/a,baryon_fraction_bound=2*y,
                         max_thermal_fugacity=max(zn,zp,zd),
                         chemical_equilibrium_log_residual=residual))
    assert rows[0]['neutron_fraction_bound']<1e-10
    assert rows[-1]['neutron_fraction_bound']>0.9
    files=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/adoption_decision.md',
           STAGE/'1000/research_round_1000_checks.json',
           STAGE/'981/drafts/common_parent_contract_v1.md',
           STAGE/'999/overall_operation_hypothesis_v1.md']
    return dict(round=1001,date='2026-10-07',all_scientific_checks_passed=True,
        kind='conditional_nuclear_binding_kinematics_and_restricted_equilibrium',
        mass_inputs_MeV={k:num(v) for k,v in zip(('proton','neutron','deuteron'),(mp,mn,md))},
        binding_energy_MeV=num(B),capture_photon_MeV=num(photon),
        deuteron_recoil_MeV=num(recoil),photodissociation_threshold_MeV=num(threshold),
        matched_energy_residual_MeV=num(F(0)),double_count_residual_MeV=num(double_count),
        photon_equals_binding_on_shell_energy_excess_MeV=bad_recoil,
        thermal_inputs=dict(eta=eta,neutron_constituent_fraction=a,
                            proton_constituent_fraction=b,zeta3=zeta3,
                            spins_degeneracies=dict(n=2,p=2,d=3),
                            chemical_potential_photon=0),equilibrium_rows=rows,
        sources=dict(masses='https://physics.nist.gov/cuu/Constants/Table/allascii.txt',
            equilibrium='https://arxiv.org/html/1801.08023#A1.SS4',
            description='CODATA 2022 central values; Pitrou et al. A.4 equations 187-191'),
        capture_rate_derived=False,cosmic_history_simulated=False,
        physical_abundances_predicted=False,mass_uncertainty_propagated=False,
        work_extraction_certified=False,joint_record_lifecycle_certified=False,
        cognitive_derivation_of_nuclear_force=False,full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})


def compare(actual,saved):
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys()
        for k in actual:compare(actual[k],saved[k])
    elif isinstance(actual,list):
        assert len(actual)==len(saved)
        for a,b in zip(actual,saved):compare(a,b)
    elif isinstance(actual,float):assert math.isclose(actual,saved,rel_tol=2e-12,abs_tol=5e-14)
    else:assert actual==saved,(actual,saved)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();out=calculate()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as dest:
            json.dump(out,dest,ensure_ascii=False,indent=2);dest.write('\n')
    else:compare(out,json.loads(TARGET.read_text('utf-8')))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','sources')},ensure_ascii=False,indent=2))
