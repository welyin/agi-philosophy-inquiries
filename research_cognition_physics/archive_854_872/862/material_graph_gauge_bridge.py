"""862: original color, gauge-orbit Wilson witness and transported Gauss.
Analytic continuum statements are in research_note_862.md. The matrix checks
calibrate those statements; no finite graph quantum equivalence is inferred.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,sys
import numpy as np
import material_link_path_probe as inherited
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'material_graph_gauge_bridge_results.json'
sys.path.insert(0,str(HERE.parent/'861'))
import magnetic_coordinate_canonical_probe as magnetic


def gauss_shift_jet():
    gen=magnetic.generators()
    A=[F(i+2,7)*gen[i]+F(1,11)*gen[(i+1)%3] for i in range(3)]
    cov=magnetic.comm
    C=[[cov(A[i],A[j]) for j in range(3)] for i in range(3)]
    metric=np.array([[F(2),F(1,5),F(0)],[F(1,5),F(3,2),F(1,7)],[F(0),F(1,7),F(5,4)]],dtype=object)
    # The proof allows any positive metric here; determinant-one normalization
    # is unnecessary for the antisymmetric double-divergence identity.
    R=[[sum((metric[i,k]*metric[j,l]*C[k][l] for k in range(3) for l in range(3)),magnetic.zero((3,3))) for j in range(3)] for i in range(3)]
    f=F(3,5);df=[F(1,3),F(-2,7),F(1,11)]
    ddf=np.array([[F(1),F(2,3),F(1,4)],[F(2,3),F(-1),F(1,7)],[F(1,4),F(1,7),F(2)]],dtype=object)
    shift=[-F(3,2)*sum((df[i]*R[i][j]+f*cov(A[i],R[i][j]) for i in range(3)),magnetic.zero((3,3))) for j in range(3)]
    # Original constant connection, variable f: evaluate D_j D_i(f R^ij)
    # including all first/second derivatives. This is an exact continuum jet,
    # not a finite-difference covariant derivative with a failed Leibniz rule.
    divergence=magnetic.zero((3,3))
    for i in range(3):
        for j in range(3):
            divergence-=F(3,2)*(ddf[j,i]*R[i][j]+df[i]*cov(A[j],R[i][j])+df[j]*cov(A[i],R[i][j])+f*cov(A[j],cov(A[i],R[i][j])))
    assert np.count_nonzero(divergence)==0
    norm2=sum(magnetic.inner(Q,Q) for Q in shift);assert norm2>0
    return dict(full_shift_covariant_divergence='0',nonzero_shift_norm_squared=str(norm2),
                all_coefficient_derivatives_retained=True,exact_continuum_jet_not_lattice_Gauss=True)


def wilson_orbit():
    from research_layout import Layout,ResearchRuntime
    with ResearchRuntime(Layout()).installed():
        import joint_reference_constraint_strata as old
        data=old.color(1.)
    Ay=-1j*data['A'][1];Az=-1j*data['A'][2]
    Fyz=Ay@Az-Az@Ay;rho=.4
    coefficient=-rho*float(np.trace(Fyz@Fyz).real)/30
    assert coefficient>0
    exp=inherited.exp_anti
    pts,ws=np.polynomial.legendre.leggauss(64)
    rows=[]
    for ell in (.8,.4,.2):
        rest=exp(Az*ell)@exp(Ay*ell)@exp(-Az*ell)
        derivative=np.zeros((3,3),complex)
        for u,w in zip((pts+1)/2,ws/2):
            f=rho*ell*u*u*(1-u)**2
            derivative-=w*exp(-Ay*ell*(1-u))@(-Fyz)@exp(-Ay*ell*u)*ell*f
        dW=float(np.trace(rest@derivative).real)
        assert dW>0
        rows.append(dict(side=ell,closed_Wilson_derivative=dW,divided_by_side_fourth=dW/ell**4,
                         relative_small_loop_coefficient_error=abs(dW/ell**4/coefficient-1)))
    assert rows[-1]['relative_small_loop_coefficient_error']<rows[0]['relative_small_loop_coefficient_error']/10
    # Explicit local diffeomorphism theta_t(y,z)=(y,z+t f(y)) on the
    # bottom-edge tube, smoothly cut off before the other edges. f=f'=0
    # at both endpoints, so vertices and their first jets are fixed.
    ell=.8;n=512;u=np.arange(n+1)/n
    f=rho*ell*u*u*(1-u)**2;df=np.diff(f)
    rest=exp(Az*ell)@exp(Ay*ell)@exp(-Az*ell)
    baseline=rest@exp(-Ay*ell)
    def bottom(t,inverse_material=False):
        U=np.eye(3,dtype=complex)
        for delta in df:
            dz=-t*delta if inverse_material else 0.
            # theta_t^* A contributes t*df to the z increment.
            U=exp(-(Ay*ell/n+Az*(dz+t*delta)))@U
        return U
    eps=.001
    actual=(np.trace(rest@bottom(eps)).real-np.trace(rest@bottom(-eps)).real)/(2*eps)
    pred=rows[0]['closed_Wilson_derivative'];assert abs(actual-pred)<1e-9
    invariant=[]
    for t in (-.2,.1,.3):
        closed=rest@bottom(t,True)
        error=float(np.linalg.norm(closed-baseline));assert error<3e-12
        moved=float(np.trace(rest@bottom(t)).real-np.trace(baseline).real)
        assert abs(moved)>1e-8
        invariant.append(dict(parameter=t,fully_transported_material_loop_error=error,
                              fixed_embedded_loop_trace_change=moved))
    # Use a representation of the ACTUAL quotient group: d_R=(3,1)_{-2}.
    omega=np.exp(2j*np.pi/3);z=np.exp(1j*np.pi/3)
    center=float(abs(omega*z**-2-1));assert center<1e-14
    return dict(original_753_color_used=True,representation='d_R=(3,1)_{-2}',
        diagonal_Z6_center_error=center,color_fundamental_alone_not_claimed_a_quotient_observable=True,
        small_loop_leading_coefficient=coefficient,small_loop_rows=rows,
        finite_difference_closed_trace_derivative=float(actual),derivative_error=abs(actual-pred),
        pure_spatial_pullback_rows=invariant,
        matrix_calibration_abelian_factor_is_identity=True,
        full_original_abelian_and_weak_fields_retained_in_analytic_argument=True,
        original_Einstein_constraints_preserved_by_diffeomorphism='analytic covariance, not a new constraint solve',
        actual_global_material_grid_numerically_reconstructed=False)


def run():
    assert inherited.run()==json.loads(inherited.TARGET.read_text('utf-8'))
    return dict(round=862,fresh_test_groups=1,all_checks_passed=True,
        inherited_path_probe_reproduced=True,gauss_shift_jet=gauss_shift_jet(),wilson_gauge_orbit=wilson_orbit(),
        analytic_scope='The original local classical reference reduction admits relational holonomies and transported Gauss/continuum flux brackets. A nonzero quotient-group Wilson-loop variation on a pure diffeomorphism orbit excludes identifying untouched fixed-embedding graph observables with reduced observables.',
        finite_graph_TstarG_Poisson_reduction_proved=False,
        full_quantum_graph_measure_or_state_equivalence_proved=False,
        reduced_Hamiltonian_self_adjointness_or_finite_coupling_proved=False,
        full_goal_completed=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
