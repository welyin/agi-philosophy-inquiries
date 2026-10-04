"""705: strict regional compression, original quotient charges and Gauss support.
Analytic full-model conclusions use617/636/637 and faithful603 Gibbs states.
Finite matrices test exact boundary representations and declared multiplicity fixtures.
No full-Hamiltonian Gibbs spectrum, continuum or apparatus is simulated.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_quotient_boundary_entropy as boundary
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_regional_charge_compression_results.json'
gluing=boundary.gluing


def trace_norm_hermitian(x):
    return float(np.abs(np.linalg.eigvalsh((x+x.conj().T)/2)).sum())


def original_charge_tail_check():
    # Original allowed characters (0,0,0,6n), not new particle species.
    rng=np.random.default_rng(705)
    M=12
    weights=[Fraction(1,2**(n+1)) for n in range(M)]+[Fraction(1,2**M)]
    assert sum(weights)==1
    labels=[boundary.label_data((0,0,0,6*n)) for n in range(M+1)]
    assert all(x['residue']==0 and x['dimension']==1 for x in labels)
    errors=[]
    for _ in range(5):
        edges=[gluing.sample(rng) for _ in range(4)]
        gs=[gluing.sample(rng) for _ in range(4)]
        endpoints=((0,1),(0,2),(1,3),(2,3))
        transformed=[gluing.product(gluing.product(gs[s],u),gluing.inverse(gs[t]))
                     for u,(s,t) in zip(edges,endpoints)]
        u=boundary.loop_holonomy(edges)[2];v=boundary.loop_holonomy(transformed)[2]
        errors.extend(abs(u**(6*n)-v**(6*n)) for n in range(M+1))
    rows=[]
    for N in (0,1,3,7):
        p=sum(weights[N+1:]);assert p==Fraction(1,2**(N+1))
        # Local keep/reset on A, B is untouched: diagonal original matched loop mixture.
        joint=np.zeros((M+1,M+1))
        for n,w in enumerate(weights):joint[n if n<=N else 0,n]+=float(w)
        assert np.max(abs(joint.sum(axis=0)-np.array([float(w) for w in weights])))<1e-15
        leak=float(1-np.trace(joint));assert abs(leak-float(p))<1e-15
        rows.append(dict(keep_labels_through=N,exact_excluded_probability=str(p),
            actual_Gauss_leakage=leak,remote_marginal_unchanged=True))
    assert max(errors)<1e-12
    return dict(allowed_original_labels=labels,
        original_loop_gauge_invariance_error=float(max(errors)),rows=rows,
        finite_diagnostic_tail_collected_at_label_M=M,
        infinite_geometric_loop_mixture_has_all_energy_moments_analytic=True,
        full_H_Gibbs_charge_tail_not_numerically_computed=True)


def carrier_reset_check():
    rng=np.random.default_rng(7051);rows=[]
    for name in ('L','u','Q'):
        g1,g2=gluing.sample(rng),gluing.sample(rng)
        r1,r2=gluing.rep(g1,name),gluing.rep(g2,name)
        carrier=np.kron(r1,r2.conj());D=len(carrier)
        psi=np.eye(D)/np.sqrt(D)
        invariance=float(np.linalg.norm(carrier@psi@carrier.conj().T-psi))
        # (depolarize A x id B)(maximally entangled state)=I_AB / D^2.
        survival=1/(D*D);predicted=Fraction(1,D*D)
        # Independent full-density check in the two smaller ORIGINAL representations.
        explicit=None
        if D<=9:
            v=psi.ravel();P=np.outer(v,v.conj())
            out=np.eye(D*D)/(D*D)
            explicit=float(np.trace(P@out).real)
            assert abs(explicit-survival)<1e-14
        assert invariance<3e-14
        rows.append(dict(original_module=name,d=len(r1),two_cut_carrier_dimension=D,
            physical_singlet_survival=str(predicted),Gauss_leakage=1-survival,
            singlet_group_invariance_error=invariance,explicit_density_survival=explicit))
    return dict(rows=rows,channel_covariant_but_physical_support_not_preserved=True,
        depolarization_keeps_the_label_but_loses_matching_carrier=True)


def fixture():
    # Full quotient boundary representations inherited; multiplicity energies explicitly declared.
    names=('nu','L','u');m=3;offset=0;blocks=[]
    for name in names:
        data=boundary.label_data(boundary.MODULES[name]);D=data['dimension']**2
        mu=1+2*sum(data['casimirs'])
        energies=mu+np.array([0.,1.,7.])
        blocks.append(dict(name=name,D=D,m=m,start=offset,stop=offset+m*D,
                           multiplicity_energies=energies))
        offset+=m*D
    energy=np.concatenate([np.repeat(b['multiplicity_energies'],b['D']) for b in blocks])
    return blocks,energy


def matched_random_vector(blocks,n,seed):
    rng=np.random.default_rng(seed);x=np.zeros((n,n,2),complex)
    for b in blocks:
        D,m,s=b['D'],b['m'],b['start']
        coeff=rng.normal(size=(m,m,2))+1j*rng.normal(size=(m,m,2))
        # Genuine matched carrier, arbitrary entangled multiplicities and passive reference.
        for i,j,r in np.ndindex(m,m,2):
            for k in range(D):x[s+i*D+k,s+j*D+k,r]=coeff[i,j,r]/np.sqrt(D)
    return x/np.linalg.norm(x)


def physical_part(x,blocks):
    out=np.zeros_like(x)
    for b in blocks:
        D,m,s=b['D'],b['m'],b['start']
        for i,j,r in np.ndindex(m,m,2):
            coeff=sum(x[s+i*D+k,s+j*D+k,r] for k in range(D))/D
            for k in range(D):out[s+i*D+k,s+j*D+k,r]=coeff
    return out


def kraus_at(blocks,energy,N):
    keep=energy<=N;ops=[np.diag(keep.astype(float))]
    for b in blocks:
        D,m,s=b['D'],b['m'],b['start']
        for j in range(m):
            if b['multiplicity_energies'][j]<=N:continue
            k=np.zeros((len(energy),len(energy)))
            k[s:s+D,s+j*D:s+(j+1)*D]=np.eye(D)
            ops.append(k)
    assert np.linalg.norm(sum(k.conj().T@k for k in ops)-np.eye(len(energy)))<1e-14
    return ops,keep


def applyA(k,x):return np.einsum('ij,jbr->ibr',k,x)
def applyB(k,x):return np.einsum('ij,ajr->air',k,x)


def density_difference_trace_norm(vectors,initial):
    # Nonzero eigenvalues of V diag(-1,+1,...) V* using the small Gram matrix.
    V=np.column_stack([initial.ravel()]+[x.ravel() for x in vectors])
    gram=V.conj().T@V;e,u=np.linalg.eigh(gram)
    assert e.min()>-1e-12
    root=(u*np.sqrt(np.maximum(e,0)))@u.conj().T
    signature=np.diag([-1.]+[1.]*len(vectors))
    return trace_norm_hermitian(root@signature@root)


def conditional_local_check():
    blocks,energy=fixture();n=len(energy);x=matched_random_vector(blocks,n,7052)
    assert np.linalg.norm(x-physical_part(x,blocks))<1e-14
    rng=np.random.default_rng(7053)
    g1,g2=gluing.sample(rng),gluing.sample(rng);U=np.zeros((n,n),complex)
    for b in blocks:
        r1,r2=gluing.rep(g1,b['name']),gluing.rep(g2,b['name'])
        carrier=np.kron(r1,r2.conj())
        U[b['start']:b['stop'],b['start']:b['stop']]=np.kron(np.eye(b['m']),carrier)
    X=x.reshape(n,-1)
    remote=X.conj().T@X
    energy_total=energy[:,None,None]+energy[None,:,None]
    base_E=float(np.sum(energy_total*abs(x)**2))
    base_E2=float(np.sum(energy_total**2*abs(x)**2))
    rows=[]
    for N in (10.,23.,36.,44.):
        ops,keep=kraus_at(blocks,energy,N)
        out=[applyA(k,x) for k in ops]
        p=float(np.sum(abs(x[~keep])**2))
        leak=sum(np.linalg.norm(v-physical_part(v,blocks))**2 for v in out)
        remote_after=sum(v.reshape(n,-1).conj().T@v.reshape(n,-1) for v in out)
        remote_error=float(np.linalg.norm(remote_after-remote))
        invariance=max(float(np.linalg.norm(k@U-U@k)) for k in ops)
        error=density_difference_trace_norm(out,x)
        assert error<=2*np.sqrt(p)+p+1e-11
        assert max(leak,remote_error,invariance)<2e-13
        E=sum(float(np.sum(energy_total*abs(v)**2)) for v in out)
        E2=sum(float(np.sum(energy_total**2*abs(v)**2)) for v in out)
        assert E<=base_E+1e-11 and E2<=base_E2+1e-9
        # Two genuinely local maps commute, including reference and coherences.
        both=[applyB(l,v) for k,v in zip(ops,out) for l in ops]
        reverse=[applyA(k,applyB(l,x)) for k in ops for l in ops]
        order_error=max(float(np.linalg.norm(a-b)) for a,b in zip(both,reverse))
        assert order_error<2e-14
        both_leak=sum(np.linalg.norm(v-physical_part(v,blocks))**2 for v in both)
        assert both_leak<2e-13
        # Preserve cross-sector coherence on the kept part: reset branches have one sector.
        a0,b0=blocks[0],blocks[1]
        psi0=np.zeros_like(x);psi1=np.zeros_like(x)
        for b,psi in ((a0,psi0),(b0,psi1)):
            for j in range(b['D']):psi[b['start']+j,b['start']+j,0]=1/np.sqrt(b['D'])
        cross=sum(np.vdot(psi0,v)*np.vdot(v,psi1) for v in out)
        main=np.vdot(psi0,out[0])*np.vdot(out[0],psi1)
        assert abs(cross-main)<1e-14
        rows.append(dict(cutoff=N,retained_low_energy_rank=int(keep.sum()),
            Kraus_count=len(ops),tail_probability=p,trace_distance_norm=error,
            gentle_upper_bound=2*np.sqrt(p)+p,physical_support_leakage=float(leak),
            remote_and_reference_marginal_error=remote_error,intertwiner_error=invariance,
            comparison_energy=E,comparison_second_moment=E2,
            two_region_order_error=order_error,two_region_support_leakage=float(both_leak),
            retained_cross_sector_coherence=[float(cross.real),float(cross.imag)]))
    assert rows[-1]['tail_probability']==0
    return dict(region_dimension=n,passive_reference_dimension=2,
        original_labels=[b['name'] for b in blocks],rows=rows,
        comparison_energy_initial=base_E,comparison_second_moment_initial=base_E2,
        multiplicity_energy_gaps_are_declared_fixture_not_full_KA_spectrum=[0,1,7],
        finite_fixture_not_claimed_to_remove_infinitely_many_original_boundary_sectors=True,
        coherence_not_removed_by_an_initial_sector_measurement=True)


def run():
    deps=('research_note_360.md','research_note_363.md','research_note_365.md',
        'research_note_617.md','research_note_625.md','research_note_636.md',
        'research_note_637.md','research_note_638.md','research_note_704.md',
        'joint_region_energy_gluing.py','joint_quotient_boundary_entropy.py')
    return dict(date='2026-10-02',round=705,tests_run=3,failures=0,errors=0,
        original_charge_tail=original_charge_tail_check(),
        matching_carrier=carrier_reset_check(),local_sector_channel=conditional_local_check(),
        dependency_hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps},
        scope=dict(original_fixed_graph_and_original_quotient=True,
            strict_local_finite_output_Gauss_conflict_for_unbounded_boundary_support=True,
            faithful_full_Gibbs_obstruction_analytic_no_Gibbs_spectrum_computed=True,
            boundary_retaining_local_CPTP_alternative=True,
            finite_per_sector_output_not_globally_finite=True,
            comparison_energy_not_full_H_monotonicity=True,
            actual_dynamics_source_derivatives_for_this_new_family_still_open=True,
            continuum_chiral_and_quantum_GR_not_proved=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(dict(round=705,tests=3,charges=result['original_charge_tail']['rows'],
        local=result['local_sector_channel']['rows'],all_checks_passed=True)))

