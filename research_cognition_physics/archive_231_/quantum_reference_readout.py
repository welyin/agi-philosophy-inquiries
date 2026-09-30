"""A quantum clock zero mode and an actual probe on the supplied FLRW background.

The background and homogeneous-mode truncation are inputs. This is not a
self-consistent quantum Einstein solution or an autonomous total-universe model.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import material_reference_geometry as bg

HERE=Path(__file__).resolve().parent
TARGET=HERE/"quantum_reference_readout_results.json"
D,Q,B,AF=3,0.6,0.4,2.0
TA,TB=0.12,0.55
WIDTH=TB-TA
LIP=np.pi/WIDTH
TF=float(bg.clock_of_a(AF,D,Q,B))


def f(x):
    z=(np.asarray(x)-TA)/WIDTH
    return np.where((z>0)&(z<1),np.sin(np.pi*z)**2,0.)


def fp(x):
    z=(np.asarray(x)-TA)/WIDTH
    return np.where((z>0)&(z<1),np.pi/WIDTH*np.sin(2*np.pi*z),0.)


def quad(fun,lo,hi,n=128):
    x,w=np.polynomial.legendre.leggauss(n)
    nodes=(lo+hi)/2+(hi-lo)*x/2
    return float((hi-lo)/2*np.dot(w,fun(nodes)))


GATE_AREA=quad(lambda z:f(z)*bg.a_of_clock(z,D,Q,B)**D/Q,TA,TB)
COUPLING=np.pi/(4*GATE_AREA)
END_TIME=quad(lambda a:1/bg.adot(a,D,Q,B),1,AF)


def analytic_bounds(volume):
    variance=TF/(2*volume*Q)
    width=lambda a:np.sqrt(variance+
        (bg.clock_of_a(a,D,Q,B)/(2*volume*Q*np.sqrt(variance)))**2)
    delta=COUPLING*LIP*quad(lambda a:width(a)/bg.adot(a,D,Q,B),1,AF)
    impulse=COUPLING*LIP*END_TIME
    p0=volume*Q
    pvar=1/(4*variance)
    # Elementary time and gate-area bounds supplement the evaluated integral.
    time_upper=(AF-1)/np.sqrt(B*B/(D-1))
    amin=float(bg.a_of_clock(TA+WIDTH/4,D,Q,B))
    lambda_upper=np.pi*Q/(WIDTH*amin**D)
    delta_upper=lambda_upper*LIP*time_upper*np.sqrt(TF/(volume*Q))
    return dict(initial_Q_variance=variance,duhamel_integral_evaluated=delta,
        elementary_duhamel_upper=delta_upper,impulse_upper=impulse,
        displacement_upper=impulse*TF/(volume*Q),
        initial_clock_kinetic_energy=(p0*p0+pvar)/(2*volume),
        free_variance_density_at_end=pvar/(2*volume*volume*AF**(2*D)),
        recoil_density_upper_at_end=(2*np.sqrt(p0*p0+pvar)*impulse+impulse**2)
                    /(2*volume*volume*AF**(2*D)),
        interaction_density_upper_at_end=COUPLING/(volume*AF**D))


def pointer_instrument(psi,spacing):
    """Construct actual Kraus diagonals from clock vectors and an X probe read."""
    pplus,pminus=psi
    a0=(pplus+pminus)/2
    a1=1j*(pplus-pminus)/2
    kraus=[]
    choi=np.zeros((8,8),complex)  # output=(record,S), input=R
    completeness=np.zeros(2)
    wrong=[]
    for eta in (1,-1):
        v=np.column_stack([(a0+eta*a1)/np.sqrt(2),
                           (a0-eta*a1)/np.sqrt(2)])*np.sqrt(spacing)
        kraus.append(v)
        completeness+=np.sum(abs(v)**2,axis=0)
        rec=0 if eta==1 else 1
        wrong.append(float(np.sum(abs(v[:,1 if eta==1 else 0])**2)))
        w=np.zeros((len(v),8),complex)
        w[:,rec*4]=v[:,0]/np.sqrt(2)
        w[:,rec*4+3]=v[:,1]/np.sqrt(2)
        choi+=w.T@w.conj()
    ideal=np.zeros((8,8),complex)
    ideal[0,0]=0.5
    ideal[7,7]=0.5
    cdist=float(np.sum(abs(np.linalg.eigvalsh(choi-ideal)))/2)
    assert np.max(abs(completeness-1))<1e-10
    assert np.min(np.linalg.eigvalsh(choi))>-1e-12
    assert abs(np.trace(choi)-1)<1e-10
    return dict(normalized_choi_trace_distance=cdist,
                wrong_outcome_probabilities=wrong,
                completeness_residual=float(np.max(abs(completeness-1))),
                choi_min_eigenvalue=float(np.min(np.linalg.eigvalsh(choi))))


def evolve(volume,nx=1024,steps=2048,box=1.):
    """Strang propagation in y=Q-T(a); time dependence of this frame is passive."""
    y=(np.arange(nx)-nx/2)*box/nx
    dy=box/nx
    p=2*np.pi*np.fft.fftfreq(nx,d=dy)
    sigma2=TF/(2*volume*Q)
    initial=np.exp(-y*y/(4*sigma2)).astype(complex)
    initial/=np.sqrt(dy*np.sum(abs(initial)**2))
    psi=np.repeat(initial[None,:],2,axis=0)
    eigen=np.array([1.,-1.])[:,None]
    p0=volume*Q
    momentum=p+p0
    def p_moments(state):
        ft=np.fft.fft(state,axis=1,norm="ortho")
        return (dy*np.sum(abs(ft)**2*momentum,axis=1),
                dy*np.sum(abs(ft)**2*momentum**2,axis=1))
    pstart,p2start=p_moments(psi)
    initial_energy=p2start/(2*volume)+eigen[:,0]*COUPLING*dy*np.sum(
                   abs(psi)**2*f(y),axis=1)
    work=np.zeros(2)
    p2old=p2start.copy()
    da=(AF-1)/steps
    norm_error=0.
    edge_mass=0.
    for index in range(steps):
        left=1+index*da
        right=left+da
        mid=(left+right)/2
        rate=bg.adot(mid,D,Q,B)
        clock=float(bg.clock_of_a(mid,D,Q,B))
        kinetic=np.exp(-1j*da*p*p/(4*volume*mid**D*rate))
        psi=np.fft.ifft(np.fft.fft(psi,axis=1,norm="ortho")*kinetic,
                       axis=1,norm="ortho")
        psi*=np.exp(-1j*eigen*da*COUPLING*f(y+clock)/rate)
        psi=np.fft.ifft(np.fft.fft(psi,axis=1,norm="ortho")*kinetic,
                       axis=1,norm="ortho")
        pmean,p2new=p_moments(psi)
        # Exact a^-d coefficient increment, trapezoid only on the evolving <P^2>.
        work+=(p2old+p2new)*(left**(-D)-right**(-D))/(4*volume)
        p2old=p2new
        if index%64==0 or index==steps-1:
            norm_error=max(norm_error,float(np.max(abs(dy*np.sum(abs(psi)**2,axis=1)-1))))
            edge_mass=max(edge_mass,float(np.max(dy*np.sum(abs(psi[:,abs(y)>box*.4])**2,axis=1))))
    s=TF/(volume*Q)
    free=np.fft.ifft(np.fft.fft(initial,norm="ortho")*np.exp(-1j*p*p*s/2),
                    norm="ortho")
    ideal=np.exp(-1j*eigen*np.pi/4)*free[None,:]
    # Orthogonal B eigenspaces make this the operator norm of the two isometries.
    iso=float(np.max(np.sqrt(dy*np.sum(abs(psi-ideal)**2,axis=1))))
    pmean,p2end=p_moments(psi)
    means=dy*np.sum(abs(psi)**2*y,axis=1)
    energy=p2end/(2*volume*AF**D)+eigen[:,0]*COUPLING*dy*np.sum(
                       abs(psi)**2*f(y+TF),axis=1)
    energy_res=float(np.max(abs(energy-initial_energy+work)))
    var_free=float(dy*np.sum(abs(free)**2*y*y))
    expected_free_var=sigma2+(TF/(2*volume*Q*np.sqrt(sigma2)))**2
    answer=dict(volume=volume,grid=nx,steps=steps,box=box,
                joint_isometry_error=iso,
                momentum_changes=[float(x) for x in pmean-pstart],
                Q_displacements_from_free=[float(x) for x in means],
                clock_density_changes_from_free=[float(x) for x in
                    (p2end-p2start)/(2*volume*volume*AF**(2*D))],
                energy_work_residual=energy_res,
                energy_work_residual_per_volume=energy_res/volume,
                norm_residual=norm_error,max_edge_probability=edge_mass,
                free_variance_error=abs(var_free-expected_free_var),
                probe_instrument=pointer_instrument(psi,dy))
    bounds=analytic_bounds(volume)
    assert iso<=bounds["duhamel_integral_evaluated"]
    assert max(abs(x) for x in answer["momentum_changes"])<=bounds["impulse_upper"]
    assert max(abs(x) for x in answer["Q_displacements_from_free"])<=bounds["displacement_upper"]
    assert max(abs(x) for x in answer["clock_density_changes_from_free"])<=bounds["recoil_density_upper_at_end"]
    assert answer["probe_instrument"]["normalized_choi_trace_distance"]<=iso
    assert edge_mass<1e-12 and norm_error<1e-10 and answer["free_variance_error"]<1e-12
    assert energy_res/volume<1e-8
    return answer,psi,y


def run():
    assert 0<TA<TB<TF<bg.clock_parameters(D,Q,B)[2]
    area2=quad(lambda z:f(z)*bg.a_of_clock(z,D,Q,B)**D/Q,TA,TB,256)
    assert abs(area2-GATE_AREA)<1e-12
    assert abs(COUPLING*GATE_AREA-np.pi/4)<1e-14
    rows=[]
    waves={}
    for volume in (1024,4096,16384):
        item,psi,y=evolve(volume)
        item["analytic_bounds"]=analytic_bounds(volume)
        rows.append(item)
        waves[volume]=psi
    # Independent temporal and spatial/domain refinements; no projection/reset.
    temporal,psi_t,yt=evolve(4096,steps=4096)
    expanded,psi_x,yx=evolve(4096,nx=2048,steps=2048,box=2.)
    err_t=float(np.max(np.sqrt(np.sum(abs(psi_t-waves[4096])**2,axis=1)/1024)))
    grid_refined,psi_g,yg=evolve(4096,nx=2048,steps=2048,box=1.)
    err_g=float(np.max(np.sqrt(np.sum(abs(psi_g[:,::2]-waves[4096])**2,axis=1)/1024)))
    # Same lattice spacing, central 1024 points represent the same Q values.
    err_x=float(np.max(np.sqrt(np.sum(abs(psi_x[:,512:1536]-waves[4096])**2,axis=1)/1024)))
    assert err_t<1e-6 and err_x<1e-8 and err_g<1e-6
    # The canonical reference-energy variance is a source, not zero.
    controls=[]
    for volume in (1024,4096,16384):
        bound=analytic_bounds(volume)
        assert bound["free_variance_density_at_end"]>0
        assert bound["duhamel_integral_evaluated"]<=bound["elementary_duhamel_upper"]
        controls.append(dict(volume=volume,**bound))
    scaling=float(controls[0]["duhamel_integral_evaluated"]/
                  controls[2]["duhamel_integral_evaluated"])
    assert abs(scaling-4)<1e-12
    return dict(date="2026-09-30",diagnostic_tests=6,failures=0,errors=0,
        numbered_round_created=False,numbered_test_increment=0,
        dependency_sha256={"material_reference_geometry.py":
                hashlib.sha256((HERE/"material_reference_geometry.py").read_bytes()).hexdigest()},
        inputs=dict(d=D,q=Q,b=B,a_initial=1.,a_final=AF,clock_final=TF,
                    proper_time_final=END_TIME,window=[TA,TB],coupling=COUPLING,
                    coupling_lipschitz=LIP,gate_area=GATE_AREA),
        diagnostics=dict(
            background_and_gate=dict(area_quadrature_difference=abs(area2-GATE_AREA),
                                      target_angle=COUPLING*GATE_AREA),
            joint_propagation=rows,
            finite_instrument=[r["probe_instrument"] for r in rows],
            momentum_displacement_energy=[dict(volume=r["volume"],
                  momentum=r["momentum_changes"],displacement=r["Q_displacements_from_free"],
                  energy_work_residual=r["energy_work_residual"]) for r in rows],
            discretization=dict(time_step_state_difference=err_t,
                larger_box_same_spacing_state_difference=err_x,
                fixed_box_finer_grid_state_difference=err_g,
                refined_instrument=temporal["probe_instrument"]),
            resource_scaling=dict(bounds=controls,
                 duhamel_ratio_1024_to_16384=scaling)),
        scope=dict(supplied_background=True,homogeneous_mode_truncation=True,
                   complete_reference_target_probe_unitary=True,
                   arbitrary_target_and_passive_reference_bound=True,
                   fully_autonomous_total_model=False,
                   quantum_reference_stress_equals_classical_exactly=False,
                   self_consistent_einstein_backreaction=False,
                   field_localization_or_spatial_dimension_derived=False,
                   full_unification_complete=False))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as out:
            out.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    if args.check:
        assert json.loads(TARGET.read_text("utf8"))==json.loads(json.dumps(result))
    print(json.dumps(result,ensure_ascii=False,indent=2))
