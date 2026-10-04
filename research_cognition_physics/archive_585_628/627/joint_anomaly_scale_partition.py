"""627: original Weyl modules constrain retained/eliminated EFT partitions.

Exact anomaly arithmetic is a necessary test for the declared continuum
target, not an anomaly calculation on the finite CAR Hilbert space.
The mass diagnostic keeps the original 598 scalar potential and Yukawas.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_gauge_matter_constraints as old
import joint_quantum_response_matching as mass

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_anomaly_scale_partition_results.json'
NAMES=('Q','uc','dc','L','ec','nuc')
LABELS=('SU3_cubed','SU3_squared_Q','SU2_squared_Q','Q_cubed','gravity_squared_Q')
EXPECTED={
    'Q':(2,2,3,6,6),'uc':(-1,-4,0,-192,-12),
    'dc':(-1,2,0,24,6),'L':(0,0,-3,-54,-6),
    'ec':(0,0,0,216,6),'nuc':(0,0,0,0,0)}


def component_rows(name):
    # Explicit diagonal generators on the 16 internal left-handed components.
    # Color uses diag(1,1,-2); weak uses diag(1,-1).
    charge=6*old.fields(3)[name][1]
    assert charge.denominator==1
    color=(1,1,-2) if name=='Q' else (-1,-1,2) if name in ('uc','dc') else (0,)
    weak=(1,-1) if name in ('Q','L') else (0,)
    return [(c,w,int(charge)) for c in color for w in weak]


def coefficients(name):
    rows=component_rows(name)
    # Tr fundamental color generator^3=-6, ^2=6; weak ^2=2.
    sums=(sum(c**3 for c,w,q in rows)//(-6),
          sum(c*c*q for c,w,q in rows)//6,
          sum(w*w*q for c,w,q in rows)//2,
          sum(q**3 for c,w,q in rows),
          sum(q for c,w,q in rows))
    assert sums==EXPECTED[name]
    return sums


def subset_check():
    coeff={n:coefficients(n) for n in NAMES}
    allrows=[];passed=[]
    for flags in itertools.product((0,1),repeat=6):
        selected=[n for n,x in zip(NAMES,flags) if x]
        local=[sum(x*coeff[n][j] for n,x in zip(NAMES,flags)) for j in range(5)]
        doublets=3*flags[0]+flags[3]
        valid=not any(local) and doublets%2==0
        if valid:passed.append(selected)
        # Independent closed classification from three anomaly equations.
        condition=len(set(flags[:5]))==1
        assert valid==condition
        allrows.append(dict(modules=selected,local_coefficients=local,
                            usual_SU2_doublet_parity=doublets%2,
                            passes_listed_necessary_anomaly_tests=valid))
    assert len(passed)==4
    quarks=[sum(coeff[n][j] for n in ('Q','uc','dc')) for j in range(5)]
    leptons=[sum(coeff[n][j] for n in ('L','ec','nuc')) for j in range(5)]
    assert quarks==[0,0,3,-162,0] and leptons==[0,0,-3,162,0]
    return dict(normalization='Q=6Y; non-Abelian quadratic fundamental index normalized to 1',
                component_count=sum(len(component_rows(n)) for n in NAMES),
                table={n:dict(zip(LABELS,coeff[n])) for n in NAMES},
                all_64_subsets=allrows,passing_subsets=passed,
                quark_removed_coefficients=quarks,lepton_retained_coefficients=leptons,
                both_sectors_have_odd_usual_SU2_doublet_count=True,
                classification_is_necessary_not_a_full_quotient_group_anomaly_certificate=True)


def actual_mass_domain_check():
    def block_eigen(phi,ids):
        h,d=mass.matter.mass_matrices(phi)
        h=h[np.ix_(ids,ids)];d=d[np.ix_(ids,ids)]
        B=np.block([[h,d],[d.conj().T,-h.T]])
        return np.linalg.eigvalsh(B)
    rows=[];max_error=0.
    for scale in (1.,.5,.1,.01,0.):
        phi=np.array([0.,scale*mass.HIGGS,0.,0.,mass.S0])
        F=float(mass.matter.original.F(phi))
        q=block_eigen(phi,list(range(24)))
        expected=np.array(sorted(
            [-abs(mass.matter.Y['u'])*scale*mass.HIGGS/np.sqrt(F)]*12+
            [-abs(mass.matter.Y['d'])*scale*mass.HIGGS/np.sqrt(F)]*12+
            [abs(mass.matter.Y['u'])*scale*mass.HIGGS/np.sqrt(F)]*12+
            [abs(mass.matter.Y['d'])*scale*mass.HIGGS/np.sqrt(F)]*12))
        error=float(np.max(abs(q-expected)));max_error=max(max_error,error)
        assert error<1e-12 and F>0
        U=float(mass.matter.original.node_potential(phi))
        assert np.isfinite(U)
        rows.append(dict(Higgs_scale=scale,F=F,original_potential=U,
                         quark_mass_min=float(q[24]),quark_mass_max=float(q[-1]),
                         BdG_formula_error=error))
    phi=np.array([0.,mass.HIGGS,0.,0.,mass.S0])
    q=block_eigen(phi,list(range(24)));le=block_eigen(phi,list(range(24,32)))
    qmin=float(q[24]);lmax=float(le[-1])
    assert qmin<lmax and rows[-1]['quark_mass_max']==0.
    return dict(path='h=t*h0, s=s0, 0<=t<=1; same original F and U',
                rows=rows,maximum_spectral_formula_error=max_error,
                vacuum_quark_min=qmin,vacuum_lepton_max=lmax,
                no_single_mass_cut_leaves_all_original_leptons_below_all_quarks=True,
                whole_configuration_quark_gap_has_infimum_zero=True,
                full_thermal_bad_region_probability_not_computed=True,
                no_physical_pole_mass_or_full_graph_gap_claim=True)


def run():
    dependencies=('research_note_531.md','research_note_536.md','research_note_537.md',
                  'research_note_598.md','research_note_599.md','research_note_603.md',
                  'research_note_605.md','research_note_613.md','research_note_614.md',
                  'research_note_615.md','research_note_626.md',
                  'joint_fermion_gauss_completion.py','joint_gauge_matter_constraints.py')
    return dict(round=627,tests_run=2,failures=0,errors=0,
                retained_module_classification=subset_check(),
                original_mass_and_state_domain=actual_mass_domain_check(),
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest()
                                   for n in dependencies},
                scope='Necessary continuum chiral anomaly tests for complete binary subsets of the original one-generation Weyl modules, with no compensating sector: only all or none of the five charged modules, optional sterile singlet. A removed anomalous sector must leave matching gauge variation; full Wess-Zumino/global completion is not constructed. Original finite-graph mass-domain diagnostic and inherited faithful thermal support preclude a uniform positive pointwise quark mass gap on the full original state. Not a proof of the full continuum limit, all quotient-group anomalies, GR, or Standard Model uniqueness.')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(round=627,tests_run=2,passing_subsets=result['retained_module_classification']['passing_subsets'],
                         mass_domain=result['original_mass_and_state_domain']),ensure_ascii=False,indent=2))
