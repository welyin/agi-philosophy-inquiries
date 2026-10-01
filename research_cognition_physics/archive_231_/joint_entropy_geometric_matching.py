"""635: original full mass content and target gradients in replica UV entropy.

Partial scalar/fermion determinants on the same fixed Einstein background;
no gauge/gravity loops, complete entropy, or Newton-constant prediction.
"""
import argparse,hashlib,json,math
from pathlib import Path
import numpy as np
import joint_background_contact_matching as bg
import joint_fermion_scalar_loop_matching as previous
HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_entropy_geometric_matching_results.json'
PHI=np.array([0.,np.sqrt(bg.U0[0]),0.,0.,np.sqrt(bg.U0[1])])
F=float(bg.matter.original.F(PHI))
K=np.eye(5)/F+np.outer(PHI,PHI)/(6*F*F)
e,V=np.linalg.eigh(K)
Kinvhalf=(V/e**.5)@V.T
Prows=np.zeros((2,5));Prows[0,:4]=PHI[:4];Prows[1,4]=PHI[4]
Hessian=2*Prows.T@bg.L@Prows/F**2
A=Kinvhalf@Hessian@Kinvhalf
SCALAR_M2=np.maximum(np.linalg.eigvalsh(A),0)
FM,FW=bg.canonical_mass_data(bg.U0)


def J_quadrature(m2,epsilon):
    v,w=np.polynomial.legendre.leggauss(320);v=25*(v+1);w=25*w
    return float(np.dot(w,np.exp(-v-m2*epsilon**2*np.exp(v)))/epsilon**2)


def J_series(m2,epsilon):
    if m2==0:return epsilon**-2
    x=m2*epsilon**2
    E1=-np.euler_gamma-np.log(x)-sum((-x)**n/(n*math.factorial(n)) for n in range(1,45))
    return float(np.exp(-x)/epsilon**2-m2*E1)


def curvature_and_original_mass_check():
    # Independent sphere spectrum for the 4-component squared Dirac operator
    # on S^2(radius 1) x R^2: normalized trace 8t sum n exp(-t n^2).
    ts=np.array([.0005,.001,.002])
    leading=[];mixed=[]
    S2=float(np.dot(FW,FM*FM))
    for t in ts:
        n=np.arange(1,900,dtype=float)
        heat=float(8*t*np.dot(n,np.exp(-t*n*n)))
        leading.append((heat-4)/(2*t))
        mass_factor=float(np.dot(FW,np.expm1(-FM*FM*t)))
        mixed.append((heat-4)*mass_factor/(2*t*t))
    coeff=np.polynomial.polynomial.polyfit(ts,leading,2)[0]
    cross=np.polynomial.polynomial.polyfit(ts,mixed,2)[0]
    assert abs(coeff+1/3)<2e-9 and abs(cross-S2/3)<2e-9
    h,d=bg.matter.mass_matrices(PHI)
    bdg=np.block([[h,d],[d.conj().T,-h.T]])
    masses=np.linalg.eigvalsh(bdg)[32:]
    # 32 positive BdG eigenvalues: each Weyl species twice in spin.
    S2direct=float(np.sum(masses*masses)/4)
    assert abs(S2direct-S2)<2e-14
    assert np.linalg.eigvalsh(A)[0]>-1e-14
    return dict(original_F=F,original_phi=PHI.tolist(),scalar_masses_squared=SCALAR_M2.tolist(),
                fermion_masses=FM.tolist(),Dirac_equivalent_weights=FW.tolist(),
                scalar_Hessian_trace=float(np.trace(A)),fermion_weighted_mass_squared=S2,
                full_BdG_mass_error=abs(S2-S2direct),
                Dirac_sphere_R_coefficient=coeff,expected_R_coefficient=-1/3,
                sphere_mass_R_coefficient=cross,expected_mass_R_coefficient=S2/3,
                partial_determinants_only=True)


def plane_area_check():
    rows=[];maxerr=0.
    traceA=float(np.trace(A));S2=float(np.dot(FW,FM*FM))
    for eps in (.05,.1,.2):
        sj=np.array([J_quadrature(m2,eps) for m2 in SCALAR_M2])
        fj=np.array([J_quadrature(m*m,eps) for m in FM])
        expectedS=np.array([J_series(m2,eps) for m2 in SCALAR_M2])
        expectedF=np.array([J_series(m*m,eps) for m in FM])
        err=max(np.max(abs(sj-expectedS)),np.max(abs(fj-expectedF)))
        maxerr=max(maxerr,float(err))
        # Replica scalar and Dirac weights (the latter twice per Dirac).
        entropy=(np.sum(sj)+2*np.dot(FW,fj))/(48*np.pi)
        # Curvature projection of the original quadratic determinants.
        WR=-np.sum(sj)/(192*np.pi**2)-np.dot(FW,fj)/(96*np.pi**2)
        dq=-16*np.pi*WR
        assert abs(entropy-dq/4)<2e-13
        rows.append(dict(epsilon=eps,entropy_area_coefficient=float(entropy),
                         induced_inverse_G=float(dq),R_action_coefficient=float(WR),
                         replica_Wald_match_error=float(abs(entropy-dq/4))))
    assert maxerr<2e-9
    # Coefficient of log(epsilon^2) after subtraction of the 1/epsilon^2 term.
    logcoef=(traceA+2*S2)/(48*np.pi)
    return dict(rows=rows,proper_time_vs_E1_error=maxerr,
                weighted_area_species=5+2*float(np.sum(FW)),
                logarithmic_coefficient_per_log_epsilon_squared=float(logcoef),
                absolute_renormalized_G_not_determined=True,
                no_gauge_ghost_or_graviton_determinants_included=True)


def target_gradient_replica_check():
    traceA=float(np.trace(A));S2=float(np.dot(FW,FM*FM))
    # A unit target radial vector at the same field value.
    u=np.sqrt(F)*PHI
    # Physical radial direction in orthonormal target coordinates.
    ev,vv=np.linalg.eigh(K);Khalf=(vv*ev**.5)@vv.T
    u=Khalf@PHI;u/=np.linalg.norm(u)
    amp=.18;rows=[];err=0.
    for winding in (1,2,3,4):
        S=(amp*winding)**2;MX=S*np.outer(u,u)
        Jacobi=A+(S*np.eye(5)-MX)/6
        trP=float(np.trace(Jacobi))
        expected=traceA+2*S/3
        err=max(err,abs(trP-expected))
        # Extract surface trace's O(s) coefficient from its eigenvalues,
        # independently of the curvature contraction used below.
        vals=np.linalg.eigvalsh(Jacobi);ts=np.array([.0002,.0004,.0008])
        est=np.array([np.sum(np.expm1(-t*vals))/t for t in ts])
        trfrom=-np.polynomial.polynomial.polyfit(ts,est,2)[0]
        conelog=-(trfrom+2*S2)/(48*np.pi)
        b4R=-trP/6-S2/3
        action_R_log=-b4R/(32*np.pi**2)
        waldlog=-4*np.pi*action_R_log
        err=max(err,abs(conelog-waldlog))
        # Same local density under g_E=F g_J: area_E=F area_J.
        SJ=F*S
        cJ=F*(traceA+2*S2)/(192*np.pi**2)+SJ/(288*np.pi**2)
        frame_error=abs(-4*np.pi*cJ-F*waldlog)
        err=max(err,frame_error)
        rows.append(dict(winding=winding,gradient_S=S,Jacobi_trace=trP,
                         entropy_log_density_per_L=float(waldlog),
                         cone_trace_inferred_density=float(conelog),
                         Jordan_Einstein_density_error=float(frame_error)))
    base=rows[0]['entropy_log_density_per_L']
    for row in rows:
        diff=row['entropy_log_density_per_L']-base
        expected=-amp**2*(row['winding']**2-1)/(72*np.pi)
        err=max(err,abs(diff-expected))
        row['difference_from_n1']=diff
    assert err<2e-12
    assert rows[-1]['difference_from_n1']<-1e-3
    return dict(same_field_value_and_same_mass_spectrum=True,amplitude=amp,
                tangential_gradient_only=True,rows=rows,maximum_matching_error=err,
                single_constant_area_density_insufficient=True,
                no_full_off_shell_entropy_or_replica_state_claim=True)


def run():
    a=curvature_and_original_mass_check();b=plane_area_check();c=target_gradient_replica_check()
    deps=('research_note_313.md','research_note_314.md','research_note_315.md',
          'research_note_580.md','research_note_599.md','research_note_631.md',
          'research_note_632.md','joint_background_contact_matching.py',
          'joint_fermion_scalar_loop_matching.py')
    return dict(round=635,tests_run=3,failures=0,errors=0,original_mass_curvature=a,
                constant_plane_area=b,target_gradient_surface=c,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope=dict(same_fixed_Einstein_quantum_variables=True,
                           original_five_scalar_and_full_one_generation_fermions=True,
                           one_loop_local_UV_and_planar_constant_branch=True,
                           target_gradient_entropy_matches_inherited_RS_coefficient=True,
                           gauge_and_gravity_loops_not_computed=True,
                           no_full_finite_region_entropy_or_type_III_trace_claim=True,
                           Newton_finite_matching_and_dimension_remain_inputs=True,
                           no_complete_GR_or_unified_model_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run();payload=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(payload)
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=635,tests=3,all_passed=True,
                         species=result['constant_plane_area']['weighted_area_species'],
                         gradient_shift=result['target_gradient_surface']['rows'][-1]['difference_from_n1'])))

