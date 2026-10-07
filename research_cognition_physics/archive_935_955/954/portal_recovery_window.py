"""954: a joint tree-level portal response window, not a full SM certificate."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
ROOT=HERE.parents[2]
TARGET=HERE/"portal_recovery_window_results.json"
def read(p):return json.loads(p.read_text("utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    portal=read(STAGE/"946/portal_common_process_results.json")
    bridge=read(STAGE/"948/neutral_matter_bridge_results.json")
    weak=read(STAGE/"949/joint_weak_domain_results.json")
    assert all(r["all_scientific_checks_passed"] for r in (portal,bridge,weak))
    p=portal["parameters"]; bp=bridge["parameters"]
    lam,v,alpha,mu=(p[k] for k in ("lambda_H","v","portal_alpha","mu"))
    a=mu**2; A=2*lam*alpha**2; b=2*lam*alpha*v; c=2*lam*v*v
    gs=bp["g_source"]; cw=bp["W_mass_squared_derivative"]; cf=bp["fermion_mass_derivative"]
    assert abs(A*c-b*b)<1e-15
    checks=dict(schur_identity_error=0.,spectral_resolvent_error=0.,
        pole_residue_sum_error=0.,joint_response_tradeoff_error=0.,
        analytic_uniform_bounds_verified=True,window_endpoint_equation_error=0.)
    rows=[]
    for zeta in (1.,.3,.1,.03,0.):
        K=np.array([[a+A*zeta*zeta,b*zeta],[b*zeta,c]])
        vals,vec=np.linalg.eigh(K)
        residues=vec[1]**2
        checks["pole_residue_sum_error"]=max(checks["pole_residue_sum_error"],abs(float(residues.sum())-1))
        entries=[]
        for q in (0.,.1,.5,1.,math.sqrt(2),4.):
            t=q*q
            D=np.linalg.inv(t*np.eye(2)+K)
            baseline=1/(t+c)
            den=(t+a)*(t+c)+A*zeta*zeta*t
            relative=b*b*zeta*zeta/den
            loss=A*zeta*zeta*t/den
            checks["schur_identity_error"]=max(checks["schur_identity_error"],
                abs(float(D[1,1]-baseline-D[0,1]**2/D[0,0])))
            checks["spectral_resolvent_error"]=max(checks["spectral_resolvent_error"],
                abs(float(D[1,1]-sum(residues/(t+vals)))))
            rSS=gs*gs*D[0,0]; rSW=-gs*cw*D[0,1]
            correction=cw*cw*(D[1,1]-baseline)
            checks["joint_response_tradeoff_error"]=max(checks["joint_response_tradeoff_error"],
                abs(float(rSW*rSW-rSS*correction)))
            assert relative <= A*zeta*zeta/a+1e-15
            assert loss <= A*zeta*zeta/(math.sqrt(a)+math.sqrt(c))**2+1e-15
            entries.append(dict(Q=q,relative_matter_correction=relative,
                relative_record_self_response_loss=loss,record_W_response=float(rSW),
                record_fermion_response=float(-gs*cf*D[0,1]),
                record_self_response=float(rSS),matter_baseline_response=cw*cw*baseline))
        rows.append(dict(zeta=zeta,poles_mass_squared=vals.tolist(),Higgs_residues=residues.tolist(),
            Higgs_like_mass_squared_shift=float(vals[1]-c),rows=entries,
            all_spacelike_relative_matter_bound=A*zeta*zeta/a,
            all_spacelike_record_self_loss_bound=A*zeta*zeta/(math.sqrt(a)+math.sqrt(c))**2))
    # Fixed diagnostic task in coefficient units, not an empirical precision.
    delta=.003  # matter relative tolerance over ALL real spacelike Q
    smin=.002   # lower mixed W coefficient at Q=1
    q=1.;t=q*q;D0=(t+a)*(t+c);d=A*t;k=abs(gs*cw*b)
    zmax=math.sqrt(delta*a/A)
    discriminant=k*k-4*smin*smin*d*D0
    assert discriminant>0
    zmin=2*smin*D0/(k+math.sqrt(discriminant)) # stable smaller root
    ztest=.1
    signal=k*ztest/(D0+d*ztest*ztest)
    assert 0<zmin<ztest<zmax
    checks["window_endpoint_equation_error"]=abs(k*zmin/(D0+d*zmin*zmin)-smin)
    assert signal>smin
    # A second task that cannot meet both criteria within this alpha-only family.
    impossible_signal=.003
    max_signal=k*zmax/(D0+d*zmax*zmax)
    assert zmax < math.sqrt(D0/d) # monotone part of the cross response
    assert max_signal<impossible_signal
    # Reuse the old epsilon family: changes to absolute couplings do not
    # remove the normalized propagator distortion at its unchanged K.
    family_evidence=[]
    for row in weak["rows"]:
        assert np.max(abs(np.array(row["scalar_masses_squared"])-np.array(rows[0]["poles_mass_squared"])))<1e-12
        family_evidence.append(dict(epsilon=row["epsilon"],
            inherited_mass_squared=row["scalar_masses_squared"],
            matter_relative_correction_Q1=1/10.25))
    assert max(checks[k] for k in ("schur_identity_error","spectral_resolvent_error",
        "pole_residue_sum_error","joint_response_tradeoff_error","window_endpoint_equation_error"))<1e-13
    files=[Path(__file__),STAGE/"946/portal_common_process_results.json",
        STAGE/"948/neutral_matter_bridge_results.json",STAGE/"949/joint_weak_domain_results.json"]
    return dict(round=954,date="2026-10-07",all_scientific_checks_passed=True,
        parameters=dict(a=a,A=A,b=b,c=c,g_source=gs,cW=cw,cf=cf,
            family="alpha -> zeta*alpha; all other displayed physical parameters fixed"),
        checks=checks,rows=rows,inherited_949_evidence=family_evidence,
        declared_coefficient_window=dict(all_spacelike_matter_tolerance=delta,Q_for_signal=q,
            record_W_minimum=smin,zeta_minimum=zmin,zeta_maximum=zmax,
            chosen_zeta=ztest,chosen_signal=signal,chosen_uniform_matter_bound=A*ztest*ztest/a,
            incompatible_signal_request=impossible_signal,
            maximum_signal_under_same_uniform_tolerance=max_signal),
        scope=dict(joint_tree_level_response_window_proved=True,
            same_source_schur_tradeoff_exact=True,
            original_internal_h_B_masses_and_clock_kept=True,
            full_interacting_947_instrument_bound_reproved=False,
            all_SM_experimental_predictions_or_loops_certified=False,
            Higgs_pole_neighborhood_uniform_error_claimed=False,
            window_tolerance_is_empirical=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})

if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--write",action="store_true");args=parser.parse_args()
    if args.write:assert not TARGET.exists()
    out=run()
    if args.write:
        with TARGET.open("x",encoding="utf-8") as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write("\n")
    else:
        before=read(TARGET)
        assert before["source_hashes"]==out["source_hashes"] and before["scope"]==out["scope"]
        for key in ("zeta_minimum","zeta_maximum","chosen_signal"):
            assert abs(out["declared_coefficient_window"][key]-before["declared_coefficient_window"][key])<1e-13
    print(json.dumps({k:v for k,v in out.items() if k not in ("source_hashes","rows","inherited_949_evidence")},ensure_ascii=False,indent=2))
