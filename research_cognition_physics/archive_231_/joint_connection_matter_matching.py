"""664: match the original scalar/CAR model to declared connection branches.

The independent-connection comparison is a MINIMAL JORDAN-FRAME input,
not a consequence of spin, cognition, or the original fixed-graph model.
Finite CAR checks use a stated vacuum normal ordering and the exact original
eight-lepton-mode onsite sector; no continuum quantum torsion integral.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_fermion_gauss_completion as matter
import joint_gauss_fermion_influence as car

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_connection_matter_matching_results.json'
ETA=np.array([-1.,1.,1.,1.])
PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
CHI=dict(Q=-1,u=1,d=1,L=-1,e=1,nu=1)


def inner(a,b):return float(np.dot(ETA*a,b))
def cartan_metric(phi):return np.eye(5)/original.F(phi)


def christoffel(metric,phi):
    n=len(phi);g=metric(phi);inv=np.linalg.inv(g);dg=np.empty((n,n,n))
    # Complex step obtains first derivatives without subtractive cancellation.
    for k in range(n):
        p=phi.astype(complex);p[k]+=1e-24j
        dg[k]=metric(p).imag/1e-24
    out=np.zeros((n,n,n))
    for a in range(n):
        for b in range(n):
            for c in range(n):
                out[a,b,c]=sum(inv[a,d]*(dg[b,d,c]+dg[c,d,b]-dg[d,b,c]) for d in range(n))/2
    return out


def target_curvature(metric,phi,step):
    n=len(phi);gamma=christoffel(metric,phi);dgamma=np.zeros((n,n,n,n))
    for k in range(n):
        delta=np.eye(n)[k]*step
        dgamma[k]=(christoffel(metric,phi+delta)-christoffel(metric,phi-delta))/(2*step)
    ric=np.zeros((n,n))
    for b in range(n):
        for c in range(n):
            ric[b,c]=sum(dgamma[a,a,b,c]-dgamma[c,a,b,a] for a in range(n))
            ric[b,c]+=sum(gamma[a,a,d]*gamma[d,b,c]-gamma[a,c,d]*gamma[d,b,a]
                          for a in range(n) for d in range(n))
    return float(np.sum(np.linalg.inv(metric(phi))*ric))


def scalar_connection_check():
    rng=np.random.default_rng(664);errors=[];stationarity=[]
    for _ in range(24):
        phi=rng.normal(size=5)*.35;f=original.F(phi)
        df=rng.normal(size=4);axial=rng.normal(size=4)
        v=-3*df/(2*f);a=3*axial/f
        def lag(v,a):return -f*inner(v,v)/3+f*inner(a,a)/48-inner(v,df)-inner(a,axial)/8
        expected=3*inner(df,df)/(4*f)-3*inner(axial,axial)/(16*f)
        errors.append(abs(lag(v,a)-expected))
        for k in range(4):
            e=np.eye(4)[k]*1e-5
            stationarity.extend([abs((lag(v+e,a)-lag(v-e,a))/2e-5),
                                 abs((lag(v,a+e)-lag(v,a-e))/2e-5)])
        gradf=-phi/3
        missing=1.5*np.outer(gradf,gradf)/f**2
        errors.append(float(np.max(abs(original.metric(phi)-cartan_metric(phi)-missing))))
        # Same 619 Weyl weight gives a constant Einstein-frame axial coefficient.
        errors.append(abs(f**-2*(-3/(16*f))*f**3+3/16))
    assert max(errors)<3e-13 and max(stationarity)<2e-9
    rows=[]
    for phi in (np.zeros(5),np.array([.6,.2,-.4,.3,.8]),np.array([1.08,.36,-.72,.54,1.44])):
        f=original.F(phi);r2=float(phi@phi)
        expected=[-10/3,-20/3-7*r2/(9*f)]
        checks=[]
        for step in (2e-4,1e-4):
            actual=[target_curvature(original.metric,phi,step),target_curvature(cartan_metric,phi,step)]
            error=float(np.max(abs(np.array(actual)-expected)))
            assert error<2e-6
            checks.append(dict(step=step,curvatures=actual,max_error=error))
        rows.append(dict(phi=phi.tolist(),F=float(f),analytic_curvatures=expected,checks=checks))
    distances=[]
    for frac in (.5,.9,.99,.9999):
        distances.append(dict(radius_fraction=frac,metric_distance=float(np.sqrt(6)*np.arctanh(frac)),
                              Cartan_distance=float(np.sqrt(6)*np.arcsin(frac))))
    return dict(original_M=float(original.M),minimal_Jordan_connection_branch_is_extra_input=True,
                algebraic_elimination_error=float(max(errors)),stationarity_error=float(max(stationarity)),
                target_curvature_rows=rows,radial_distances=distances,
                Cartan_F_zero_distance=float(np.pi*np.sqrt(6)/2),
                original_hyperbolic_metric_is_not_preserved_by_minimal_connection_replacement=True)


def current_matrices():
    mats=np.zeros((4,32,32),complex)
    for name,sl in matter.SLICES.items():
        n=sl.stop-sl.start
        mats[0,sl,sl]=CHI[name]*np.eye(n)
        for a in range(3):mats[a+1,sl,sl]=np.kron(np.eye(n//2),PAULI[a])
    return mats


def contact_operator(mats):
    zero=np.zeros_like(mats[0]);js=[car.fock(c,zero) for c in mats]
    product=sum(eta*(j@j) for eta,j in zip(ETA,js))
    contraction=car.fock(sum(eta*(c@c) for eta,c in zip(ETA,mats)),zero)
    return product-contraction,product,js


def original_current_check():
    mats=current_matrices();rng=np.random.default_rng(6641);errors=[]
    for _ in range(12):
        rep=matter.representation(matter.gauge.group_exp(rng.normal(size=8),3),
                                  matter.gauge.group_exp(rng.normal(size=3),2),np.exp(1j*rng.normal()))
        errors.append(float(np.max(abs(mats@rep-rep@mats))))
    contraction=sum(eta*(c@c) for eta,c in zip(ETA,mats))
    errors.append(float(np.max(abs(contraction-2*np.eye(32)))))
    assert max(errors)<2e-13
    # Quark vacuum is an invariant conditional onsite sector. Keep every lepton.
    lm=mats[:,24:,24:];q,prod,_=contact_operator(lm)
    number=car.fock(np.eye(8),np.zeros((8,8)))
    order_error=float(np.max(abs(prod-q-2*number)))
    sep=np.zeros_like(q)
    for sl in (slice(0,4),slice(4,6),slice(6,8)):
        part=np.zeros_like(lm);part[:,sl,sl]=lm[:,sl,sl]
        sep+=contact_operator(part)[0]
    singlet=np.zeros(256,complex)
    singlet[(1<<4)+(1<<7)]=1/np.sqrt(2)
    singlet[(1<<5)+(1<<6)]=-1/np.sqrt(2)
    total=float(np.vdot(singlet,q@singlet).real)
    separate=float(np.vdot(singlet,sep@singlet).real)
    assert order_error<2e-13 and abs(total+8)<1e-12 and abs(separate)<1e-12
    # Independent diagonal witness: no quadratic CAR Hamiltonian has this term.
    vacuum=0;up=1<<6;down=1<<7;pair=up+down
    connected=float((q[pair,pair]-q[up,up]-q[down,down]+q[vacuum,vacuum]).real)
    assert abs(connected+8)<1e-12
    return dict(original_modes=32,conditional_lepton_modes=8,original_gauge_commutator_error=max(errors),
                normal_order_contraction_per_number=2,normal_order_identity_error=order_error,
                e_nu_spin_singlet_total_contact=total,separate_module_contact=separate,
                neutral_diagonal_connected_quartic=connected,
                conditional_sector_is_not_a_full_Gauss_thermal_state=True), (q,prod,sep)


def thermal_geometry_check(contacts):
    q,prod,sep=contacts;phi=car.PHI[0].copy();h,d=matter.mass_matrices(phi)
    h0=car.fock(h[24:,24:],d[24:,24:]);theta=.13;beta=.7
    def data(theta,operator):
        volume=np.exp(6*theta);v=3*operator/(16*volume)
        ham=h0+v;ev,vec=np.linalg.eigh(ham);weights=np.exp(-beta*(ev-ev.min()))
        logz=float(-beta*ev.min()+np.log(weights.sum()))
        diag=np.real(np.diag(vec.conj().T@v@vec))
        mean=float(weights@diag/weights.sum())
        return logz,mean,v
    full,mean,v=data(theta,q);bare=data(theta,np.zeros_like(q))[0]
    naive=data(theta,prod)[0];split=data(theta,sep)[0]
    source=6*beta*mean;step=2e-5
    finite=(data(theta+step,q)[0]-data(theta-step,q)[0])/(2*step)
    err=abs(source-finite)
    assert err<2e-6 and abs(full-bare)>1e-3 and abs(full-naive)>1e-3 and abs(full-split)>1e-3
    # A matching contact counterterm cancels the operator, hence every source,
    # without fitting the mass or separately adjusting a vacuum normalization.
    cancelled=float(np.max(abs((h0+v)-v-h0)))
    assert cancelled<2e-14
    return dict(phi=phi.tolist(),F=float(original.F(phi)),theta=theta,beta=beta,cell_volume=float(np.exp(6*theta)),
                original_Dirac_block_norm=float(np.linalg.norm(h[24:,24:])),
                original_Majorana_block_norm=float(np.linalg.norm(d[24:,24:])),
                logZ_original_quadratic=bare,logZ_shared_normal_ordered_contact=full,
                logZ_unordered_product=naive,logZ_separate_module_contacts=split,
                geometric_logZ_source=source,finite_difference_source=finite,source_error=err,
                frozen_contact_source_error=abs(source),counterterm_operator_error=cancelled,
                no_quantum_continuum_coefficient_or_torsion_measure_claim=True)


def run():
    currents,contacts=original_current_check()
    deps=('joint_curved_quantum_source.py','joint_fermion_gauss_completion.py','joint_gauss_fermion_influence.py',
          'research_note_358.md','research_note_375.md','research_note_580.md','research_note_598.md',
          'research_note_599.md','research_note_600.md','research_note_619.md','research_note_663.md')
    return dict(date='2026-10-02',round=664,tests_run=3,failures=0,errors=0,
                scalar_connection=scalar_connection_check(),original_currents=currents,
                conditional_thermal_source=thermal_geometry_check(contacts),
                dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
                scope='Original metric scalar/CAR branch compared to a declared minimal Jordan independent-spin-connection branch. Tree-level bulk elimination, target-metric invariant mismatch, original universal axial current, finite-cell normal ordering and conditional eight-lepton-mode geometry source. Matching counterterms restore shared classical metric/matter action, not the quantum measure. No Einstein selection, full Gauss thermal integral, continuum matching or quantum GR.',
                all_checks_passed=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args()
    result=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','all_checks_passed')}))
