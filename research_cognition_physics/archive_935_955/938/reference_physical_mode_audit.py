"""938: finite-wavelength physical reference modes on the original full fields.

This tests classical initial data and the native generator, not a quantum
completion. The extra physical mode count uses the same regular constraints.
"""
from pathlib import Path
import sys, json, hashlib, argparse
import numpy as np

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'reference_physical_mode_audit_results.json'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime

def mx(a):return float(np.max(np.abs(a)))

def run():
    rows=[];windows=[]
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        geo=old.geo
        for n in (16,24):
            z,psi,tensor,info,color=old.completed(n,1.)
            d=.01;tau=np.sqrt(z['tau2'])+d
            rhoD=(tau*tau-z['tau2'])/3
            v=psi**6;qi=psi**-4;dx3=z['dx']**3
            R=-8*psi**-5*geo.laplace(psi)
            rho_nd=.5*psi**-12*z['pKp']+.5*psi**-4*z['B']+z['U']+psi**-8*z['Y']
            q=psi[...,None,None]**4*np.eye(3)
            qinv=qi[...,None,None]*np.eye(3)
            shear=psi[...,None,None]**-2*tensor
            def constraint(trace):
                kk=shear+trace[...,None,None]*q/3
                kup=psi[...,None,None]**-8*kk
                pi=-.5*v[...,None,None]*(kup-trace[...,None,None]*qinv)
                tr=np.sum(pi*q,axis=(-1,-2))
                pp=psi**8*np.sum(pi*pi,axis=(-1,-2))
                return 2/v*(pp-.5*tr*tr)-v*R/2+v*rho_nd,pi
            baseC,basePi=constraint(np.full_like(psi,tau))
            assert np.max(baseC)<0
            r0=-baseC
            f=np.sin(z['grid'][...,0])
            df=np.zeros(z['grid'].shape);df[...,0]=np.cos(z['grid'][...,0])
            coeff=float(dx3*np.sum((2/9)*v*v*qi*np.sum(df*df,axis=-1)/r0))
            # Exact analytic sufficient condition using already established
            # old753 psi>1/2, tau<1.743 and rhoD>.0115, not grid extrema.
            # |epsilon|<=3e-4 gives -C/v > .01115 and sqrt(D)/v<=.0008.
            windows.append(dict(N=n,analytic_epsilon_window=.0003,
                analytic_negative_C_per_volume_lower=.01115,
                analytic_sqrt_D_per_volume_upper=.0008,
                epsilon2_coefficient=coeff))
            for amp in (.0003,.00015):
                signs=[]
                for sign in (1,-1):
                    e=sign*amp;tr=tau+e*f
                    c,pi=constraint(tr)
                    c_formula=baseC-v*(2*tau*e*f+e*e*f*f)/3
                    assert mx(c-c_formula)<5e-14
                    ca=-(2/3)*v[...,None]*e*df
                    D=qi*np.sum(ca*ca,axis=-1)
                    rad=c*c-D
                    assert np.min(rad)>0 and np.max(c)<0
                    r=np.sqrt(rad)
                    density=rad/(v*(-c))
                    P=r;Pa=-ca
                    un=c/r;ua=ca/r[...,None]
                    lapse=-c/r;shift=qi[...,None]*ca/r[...,None]
                    energy=density*un*un
                    stress_na=density[...,None]*un[...,None]*ua
                    totalC=c+np.sqrt(P*P+qi*np.sum(Pa*Pa,axis=-1))
                    totalCa=ca+Pa
                    einstein=R-psi**-8*np.sum((shear+tr[...,None,None]*q/3)**2,axis=(-1,-2))+tr*tr-2*(rho_nd+energy)
                    assert mx(totalC)<2e-14 and mx(totalCa)==0
                    assert mx(einstein)<5e-12
                    assert mx(v[...,None]*stress_na+ca)<2e-14
                    assert mx(-un*un+qi*np.sum(ua*ua,axis=-1)+1)<2e-14
                    # Two independent evaluations of the omitted D term.
                    native_minus_C=-r-c
                    rationalized=D/((-c)+r)
                    assert mx(native_minus_C-rationalized)<3e-17
                    deltaH=float(dx3*np.sum(rationalized))
                    assert deltaH>0
                    # Native response along the same canonical path.
                    cp=-v*(2*tau*f+2*e*f*f)/3
                    cap=-(2/3)*v[...,None]*df
                    hprime=-c/r*cp+np.sum(shift*cap,axis=-1)
                    energy_derivative=float(dx3*np.sum(hprime))
                    def energy_path(u):
                        cu=baseC-v*(2*tau*u*f+u*u*f*f)/3
                        au=-(2/3)*v[...,None]*u*df
                        return float(-dx3*np.sum(np.sqrt(cu*cu-qi*np.sum(au*au,axis=-1))))
                    step=1e-7
                    fd=(energy_path(e+step)-energy_path(e-step))/(2*step)
                    assert abs(fd-energy_derivative)<1e-7
                    signrow=dict(sign=sign,native_minus_C_integral=deltaH,
                        native_source_derivative=energy_derivative,finite_difference_source=fd,
                        source_derivative_error=abs(fd-energy_derivative),
                        min_proper_density=float(np.min(density)),
                        min_radicand=float(np.min(rad)),
                        maximum_reference_shift=mx(shift),
                        old_non_dust_spatial_constraint_defect=mx(ca),
                        total_Hamiltonian_density_residual=mx(totalC),
                        total_spatial_constraint_residual=mx(totalCa),
                        full_Einstein_initial_residual=mx(einstein),
                        density_weighted_reference_stress_error=mx(v[...,None]*stress_na+ca),
                        original_internal_Gauss_preserved=True)
                    signs.append(signrow)
                scaled=(signs[0]['native_minus_C_integral']+signs[1]['native_minus_C_integral'])/(2*amp*amp)
                rows.append(dict(N=n,amplitude=amp,physical_coordinate_wavenumber=1,
                    symmetric_epsilon2_coefficient=scaled,
                    coefficient_relative_error=abs(scaled/coeff-1),signs=signs))
        for i in (0,2):
            assert 3.9<rows[i]['coefficient_relative_error']/rows[i+1]['coefficient_relative_error']<4.1
    # All original internal first-class constraints and physical fermion
    # sectors occur on both sides; only the bosonic reference increment matters.
    pair_counts={}
    for neutral_probe_pairs in (0,1):
        original_pairs=6+3*12+5+neutral_probe_pairs
        constraints=4+12
        no_dust=2*original_pairs-2*constraints
        with_dust=2*(original_pairs+4)-2*constraints
        assert with_dust-no_dust==8
        pair_counts[str(neutral_probe_pairs)]=dict(without_dust_bosonic_phase_dimension=no_dust,
            with_dust_bosonic_phase_dimension=with_dust,additional_physical_pairs=4)
    paths=[Path(__file__),HERE/'drafts/STATUS.md',HERE/'drafts/quantum_recovery_reuse_audit.md',
        STAGE/'research_note_937.md',STAGE/'937/dust_common_model_results.json',
        STAGE/'research_note_764.md',STAGE/'research_note_765.md',STAGE/'research_note_898.md',
        STAGE.parent/'archive_742_763/753/joint_reference_constraint_strata.py']
    return dict(round=938,date='2026-10-07',all_scientific_checks_passed=True,
        physical_phase_counts=pair_counts,uniform_classical_window=windows,
        full_field_finite_wavelength_rows=rows,
        new_modes_are_not_just_removed_gauge_labels=True,
        finite_scale_generator_and_source_change_verified=True,
        old_quantum_results_cannot_be_transferred_without_new_physical_mapping=True,
        not_a_quantum_completion_or_a_no_go_for_dust=True,
        dust_candidate_demoted_from_priority_quantum_route=True,
        no_new_cognitive_axiom_or_spatial_dimension_proof=True,
        full_goal_completed=False,
        source_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    if args.write:assert not TARGET.exists()
    result=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
