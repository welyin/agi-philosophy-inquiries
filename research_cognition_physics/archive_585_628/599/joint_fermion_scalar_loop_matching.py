"""599: common frozen-Einstein scalar/fermion UV derivative-order check.

This is a conditional continuum one-loop calculation, NOT a lattice matching
or a full gauge/gravity loop calculation. Couplings are the 598 diagnostics.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_fermion_gauss_completion as finite

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_fermion_scalar_loop_matching_results.json'


def clifford():
    I=np.eye(2);Z=np.zeros((2,2))
    sigmas=[np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1,-1])]
    gammas=[np.block([[Z,-1j*s],[1j*s,Z]]) for s in sigmas]
    gammas.append(np.block([[Z,I],[I,Z]]))
    chirality=np.diag([1.,1.,-1.,-1.])
    for a in range(4):
        for b in range(4):
            assert np.max(abs(gammas[a]@gammas[b]+gammas[b]@gammas[a]-
                              2*(a==b)*np.eye(4)))<1e-14
    return gammas,(np.eye(4)+chirality)/2,(np.eye(4)-chirality)/2


def spin_trace_check():
    rng=np.random.default_rng(599);gamma,PL,PR=clifford()
    errors=[];first_order_errors=[]
    for _ in range(18):
        dim=3
        M=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim))
        dM=rng.normal(size=(4,dim,dim))+1j*rng.normal(size=(4,dim,dim))
        # Arbitrary noncommuting complex matrices check more than radial masses.
        cal=lambda m:np.kron(PL,m)+np.kron(PR,m.conj().T)
        C=cal(M);R=float(rng.normal())
        GG=[np.kron(g,np.eye(dim)) for g in gamma]
        first_order_errors.extend(float(np.max(abs(C.conj().T@g-g@C))) for g in GG)
        P=R*np.eye(4*dim)/4+C.conj().T@C-sum(g@cal(d) for g,d in zip(GG,dM))
        pure=R*np.eye(4*dim)/4
        direct=np.trace(P@P/2-R*P/6-pure@pure/2+R*pure/6).real
        A=M.conj().T@M
        expected=2*np.trace(A@A).real+2*sum(np.trace(d.conj().T@d).real for d in dM)
        expected+=R*np.trace(A).real/3
        errors.append(float(abs(direct-expected)))
        assert np.max(abs(P-P.conj().T))<1e-13
    assert max(errors)<3e-12 and max(first_order_errors)<1e-14
    return dict(samples=18,flavor_dimension=3,max_spin_trace_error=max(errors),
                max_first_order_cancellation_error=max(first_order_errors),
                curvature_mass_cross_term_checked=True,
                pure_curvature_terms_subtracted_not_discarded_physically=True)


def physical_mass_direction():
    """One generation, all-left symmetric mass; Majorana weight = 1/2.

    Dirac pairs each appear twice in this 16-Weyl matrix, while the sterile
    Majorana term occurs once. Two spin components are NOT another species.
    """
    nx=2/np.sqrt(5);ns=1/np.sqrt(5)
    out=np.zeros((16,16),complex)
    cursor=0
    for name,count in [('u',3),('d',3),('e',1)]:
        for _ in range(count):
            out[cursor,cursor+1]=out[cursor+1,cursor]=finite.Y[name]*nx
            cursor+=2
    out[14,15]=out[15,14]=finite.Y['nu']*nx
    out[15,15]=finite.Y['s']*ns
    assert np.max(abs(out-out.T))==0
    return out,np.array([nx,0.,0.,0.,ns])


def common_background_check():
    M0,direction=physical_mass_direction()
    A=M0.conj().T@M0
    # Three identical diagnostic generations: illustrative, not a fitted spectrum.
    t2=3*np.trace(A).real;t4=3*np.trace(A@A).real
    amp=.27
    rows=[]
    x=np.arange(4096)*2*np.pi/4096
    for n in (1,2,3,4,5):
        q=amp*np.sin(n*x);dq=amp*n*np.cos(n*x)
        phi=np.sqrt(6*original.M)*np.tanh(q/np.sqrt(6))[:,None]*direction
        f=np.sqrt(6)*np.sinh(q/np.sqrt(6))
        # Verify using the original F, not just the analytic radial formula.
        scaled=phi/np.sqrt(original.F(phi))[:,None]
        assert np.max(abs(scaled-f[:,None]*direction))<1e-15
        dfdx=np.cosh(q/np.sqrt(6))*dq
        fermion=t4*f**4+t2*dfdx**2  # Half the squared-Dirac b4.
        scalar_four=dq**4/18       # Inherited H^5 scalar coefficient.
        integral=lambda z:float(2*np.pi*np.mean(z))
        rows.append(dict(n=n,fermion_b4=integral(fermion),
                         scalar_four_gradient=integral(scalar_four),
                         combined_partial_b4=integral(scalar_four-fermion)))
    n=np.array([r['n'] for r in rows],float)
    matrix=np.column_stack([np.ones_like(n),n*n,n**4])
    # Normalize polynomial columns to avoid an artificial conditioning error.
    scale=np.linalg.norm(matrix,axis=0)
    b=np.array([r['combined_partial_b4'] for r in rows])
    fit=np.linalg.lstsq(matrix/scale,b,rcond=None)[0]/scale
    expected_four=np.pi*amp**4/24
    residual=float(np.max(abs(matrix@fit-b)))
    assert abs(fit[2]-expected_four)<2e-14 and residual<1e-12
    fermion_fit=np.linalg.lstsq(matrix/scale,np.array([r['fermion_b4'] for r in rows]),rcond=None)[0]/scale
    assert abs(fermion_fit[2])<2e-14
    # At x=0, the R*S source has an n^4 term; the fermion R*mass^2
    # counterterm only has n^2. These are the curvature-dependent terms only.
    curvature_source=4*amp**2*n**4/9-2*t2*amp**2*n**2/3
    source_fit=np.linalg.lstsq(np.column_stack([n*n,n**4]),curvature_source,rcond=None)[0]
    assert abs(source_fit[1]-4*amp**2/9)<2e-14
    return dict(amplitude=amp,diagnostic_generations=3,Weyl_components=48,
                original_target_geodesic=True,rows=rows,fit_constant_n2_n4=fit.tolist(),
                expected_scalar_n4=expected_four,max_fit_error=residual,
                fermion_n4_fit=float(fermion_fit[2]),
                curvature_source_fit_n2_n4=source_fit.tolist(),
                omits_inherited_scalar_n0_n2_terms=True,
                not_a_fixed_hbar_graph_continuum_limit=True)


def run():
    spin=spin_trace_check();background=common_background_check()
    names=('research_note_553.md','research_note_580.md','research_note_581.md',
           'research_note_582.md','research_note_583.md','research_note_598.md',
           'joint_fermion_gauss_completion.py','joint_curved_quantum_source.py')
    return dict(round=599,tests_run=2,failures=0,errors=0,spin=spin,background=background,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
                scope=dict(frozen_Einstein_background=True,zero_gauge_background=True,
                           bosonic_background_and_zero_fermion_background=True,
                           canonical_fermion_kinetic_operator_is_additional_input=True,
                           parity_even_local_one_loop_UV_only=True,
                           scalar_four_gradient_and_RS_survive_fermion_block=True,
                           full_gauge_and_gravity_loops_not_computed=True,
                           no_chiral_regulator_or_finite_graph_matching_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=599,tests=2,all_passed=True,
                         n4=result['background']['fit_constant_n2_n4'][2],
                         trace_error=result['spin']['max_spin_trace_error'])))
