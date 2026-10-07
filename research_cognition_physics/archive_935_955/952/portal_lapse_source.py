"""952: the same portal exchange and its static lapse insertion.

Quadratic scalar sector, fixed proper source densities and spatial metric.
No new joint detector probability, full stress tensor or Einstein solution.
"""
from pathlib import Path
import argparse, hashlib, json, math
import numpy as np
HERE=Path(__file__).resolve().parent;STAGE=HERE.parent;ROOT=HERE.parents[2]
TARGET=HERE/'portal_lapse_source_results.json'
def read(p):return json.loads(p.read_text('utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def norm(a):return float(np.linalg.norm(a,2))

def run():
    old=read(STAGE/'948/neutral_matter_bridge_results.json')
    p=old['parameters'];K=np.array(p['inherited_K']);I2=np.eye(2)
    C=np.array([[p['g_source'],0,0],
                [0,p['W_mass_squared_derivative'],p['fermion_mass_derivative']]])
    cw,cf=p['W_mass_squared_derivative'],p['fermion_mass_derivative']
    def D(q):
        q=np.asarray(q,float);return np.linalg.inv((q@q)*I2+K)
    def vertex(qout,qin):
        qout,qin=np.asarray(qout,float),np.asarray(qin,float)
        core=(qout@qout+qin@qin-qout@qin)*I2+K
        return -C.T@D(qout)@core@D(qin)@C
    rows=[]
    uniform_error=0.;reverse_error=0.
    for qo,qi in (([1,0,0],[1,0,0]),([2,0,0],[1,0,0]),([1,1,0],[1,0,0])):
        V=vertex(qo,qi)
        endpoints=-C.T@(D(qo)+D(qi))@C
        field=C.T@D(qo)@((np.dot(qo,qi))*I2+K)@D(qi)@C
        assert norm(V-endpoints-field)<1e-14
        reverse_error=max(reverse_error,norm(V-vertex(qi,qo).T))
        assert abs(V[0,1]/cw-V[0,2]/cf)<1e-14
        if qo==qi:uniform_error=max(uniform_error,norm(V+C.T@D(qi)@C))
        rows.append(dict(q_out=qo,q_in=qi,lapse_vertex=V.tolist(),
            record_W=float(V[0,1]),record_fermion=float(V[0,2]),
            endpoints_only_record_W=float(endpoints[0,1]),
            field_insertion_record_W=float(field[0,1])))
    assert uniform_error<1e-14 and reverse_error<1e-14
    flat=next(r for r in old['rows'] if r['spacelike_Q']==1.)
    assert abs(rows[0]['record_W']-flat['record_W_exchange_coefficient'])<1e-14
    assert abs(rows[0]['endpoints_only_record_W']/rows[0]['record_W']-2)<1e-14

    # A finite Galerkin test, with fields independent of y,z. It does NOT
    # assert that this finite Fourier band is an exact continuum model.
    modes=np.arange(-3,4);dim=len(modes)*2
    L=np.zeros((dim,dim),complex);An=L.copy();nmat=L.copy()
    for a,kout in enumerate(modes):
        for b,kin in enumerate(modes):
            sl=(slice(2*a,2*a+2),slice(2*b,2*b+2))
            if a==b:L[sl]=kin*kin*I2+K
            n=.5 if abs(kout-kin)==1 else 0.
            nmat[sl]=n*I2
            An[sl]=n*(kout*kin*I2+K)
    D0=np.linalg.inv(L)
    Tprime=nmat@D0+D0@nmat-D0@An@D0
    # Compare A_n to the gradient and mass bilinear in real space,
    # independently of its Fourier matrix formula.
    x=np.arange(128)*2*np.pi/128
    basis=np.exp(1j*np.outer(x,modes))
    derivative=basis*(1j*modes)[None,:]
    realspace_An=(np.kron(derivative.conj().T@(np.cos(x)[:,None]*derivative)/len(x),I2)
                 +np.kron(basis.conj().T@(np.cos(x)[:,None]*basis)/len(x),K))
    spatial_bilinear_error=norm(An-realspace_An)
    assert spatial_bilinear_error<1e-13
    J=np.zeros((dim,3),complex)
    for channel,frequency in enumerate((1,2,2)):
        for k in (-frequency,frequency):
            a=int(np.flatnonzero(modes==k)[0])
            J[2*a:2*a+2,channel]=C[:,channel]/math.sqrt(2)
    predicted=-J.conj().T@Tprime@J
    assert abs(predicted[0,1]-.5*rows[1]['record_W'])<1e-14
    assert abs(predicted[0,2]-.5*rows[1]['record_fermion'])<1e-14
    eye=np.eye(dim)
    def kernel(lam,incorrect=False):
        N=eye+lam*nmat
        return N@np.linalg.inv(L if incorrect else L+lam*An)@N
    step=2e-5
    fd=(-kernel(2*step)+8*kernel(step)-8*kernel(-step)+kernel(-2*step))/(12*step)
    variation_error=norm(fd-Tprime)
    assert variation_error<1e-10
    lam=.1
    T=kernel(lam);A=L+lam*An
    residual=T-D0-lam*Tprime
    B=nmat-An@D0
    exact_residual=lam**2*B.conj().T@np.linalg.solve(A,B)
    remainder_identity_error=norm(residual-exact_residual)
    assert remainder_identity_error<1e-13
    assert np.linalg.eigvalsh(residual).min()>-1e-13
    kmin=float(np.linalg.eigvalsh(K).min())
    universal_bound=4*lam**2/((1-abs(lam))*kmin)
    assert norm(residual)<universal_bound
    # Evaluate the actual stationary field energy, not only the inverse formula.
    current=np.array([.6,.13,-.08])
    source=J@current
    N=eye+lam*nmat
    phi=-np.linalg.solve(A,N@source)
    explicit_energy=float((.5*np.vdot(phi,A@phi)+np.vdot(N@source,phi)).real)
    eliminated_energy=float((-.5*np.vdot(source,T@source)).real)
    assert abs(explicit_energy-eliminated_energy)<1e-14
    energy_derivative=float((.5*np.vdot(phi,An@phi)+np.vdot(nmat@source,phi)).real)
    def energy(x):
        return float((-.5*np.vdot(source,kernel(x)@source)).real)
    fd_energy=(energy(lam-2*step)-8*energy(lam-step)+8*energy(lam+step)-energy(lam+2*step))/(12*step)
    energy_variation_error=abs(energy_derivative-fd_energy)
    assert energy_variation_error<1e-11
    # A spatially varying lapse has a nonzero source-mixing response;
    # constant lapse is used ONLY as a normalization/coordinate check.
    R0=-J.conj().T@D0@J;R=-J.conj().T@T@J
    Rbad=-J.conj().T@kernel(lam,incorrect=True)@J
    assert abs(R0[0,1])<1e-15
    assert abs(R[0,1])>1e-4 and abs(Rbad[0,1]-R[0,1])>1e-4
    assert norm(R-R.conj().T)<1e-14
    assert abs(R[0,1]/cw-R[0,2]/cf)<1e-14
    # A positive bound independent of a Fourier cutoff follows from
    # |<phi,A_n phi>| <= ||n||_infty <phi,L phi>, here ||n||_infty=1.
    val,vec=np.linalg.eigh(L)
    invroot=(vec/np.sqrt(val))@vec.conj().T
    form_norm=norm(invroot@An@invroot)
    assert form_norm<=1+1e-13
    files=[Path(__file__),STAGE/'948/neutral_matter_bridge_results.json',
           STAGE/'951/body_sector_bridge_results.json',
           STAGE.parent/'archive_301_341/research_note_303.md',
           STAGE.parent/'archive_585_628/research_note_588.md']
    return dict(round=952,date='2026-10-07',all_scientific_checks_passed=True,
        parameters=dict(K=K.tolist(),C=C.tolist(),lapse_profile='N(x)=1+0.1*cos(x)',
            spatial_Galerkin_modes=modes.tolist(),spatial_metric_fixed=True),
        checks=dict(uniform_lapse_identity_error=uniform_error,
            source_exchange_reciprocity_error=reverse_error,
            independent_real_space_bilinear_error=spatial_bilinear_error,
            independent_resolvent_variation_error=variation_error,
            exact_positive_remainder_identity_error=remainder_identity_error,
            direct_stationary_energy_error=abs(explicit_energy-eliminated_energy),
            stationary_energy_lapse_variation_error=energy_variation_error,
            relative_quadratic_form_norm=form_norm),
        momentum_vertices=rows,
        finite_lapse=dict(lambda_value=lam,baseline_record_W=float(R0[0,1].real),
            correct_record_W=float(R[0,1].real),correct_record_fermion=float(R[0,2].real),
            endpoints_only_record_W=float(Rbad[0,1].real),
            wrong_minus_correct_record_W=float((Rbad[0,1]-R[0,1]).real),
            leading_record_W=float(lam*predicted[0,1].real),
            kernel_linearization_remainder_norm=norm(residual),
            kernel_remainder_bound=universal_bound,
            exact_remainder_min_eigenvalue=float(np.linalg.eigvalsh(exact_residual).min()),
            stationary_energy=explicit_energy,stationary_lapse_source=energy_derivative),
        scope=dict(same_948_record_W_and_fermion_kernel=True,
            quadratic_static_exchange_lapse_source_connected=True,
            field_and_endpoint_variations_both_included=True,
            constant_lapse_is_only_normalization_check=True,
            finite_Galerkin_result_is_not_continuum_error_certificate=True,
            full_Einstein_solution_or_full_stress_certified=False,
            new_947_finite_time_probability_or_preparation_bound=False,
            full_SM_or_quantum_gravity_restored=False,full_goal_completed=False),
        source_hashes={str(p.relative_to(ROOT)):sha(p) for p in files})
def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=2e-12),(a,b)
    else:assert a==b,(a,b)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=run()
    if args.write:
        with TARGET.open('x',encoding='utf-8') as f:json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare(out,read(TARGET))
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_hashes','momentum_vertices')},ensure_ascii=False,indent=2))

