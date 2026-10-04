"""727: conditional original-matter ground reference and its retained sources.

Finite CAR spectral tests, not an interacting Gauss ground-state solver.
The same original masses, physical-right convention and604 Weyl edges are used.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as matter
import joint_quantum_response_matching as response
import joint_chiral_source_matching as chiral
import joint_recorded_classical_wall as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_matter_ground_source_results.json'
CENTER=old.source()['center']
IDS=[24,25,30,31]


def block(items):
    out=np.zeros((sum(len(a) for a in items),)*2,complex);pos=0
    for a in items:
        out[pos:pos+len(a),pos:pos+len(a)]=a;pos+=len(a)
    return out


def bdg(h,d):
    return np.block([[h,d],[d.conj().T,-h.T]])


def radial(h,s,ids=None):
    phi=np.array([0.,h,0.,0.,s])
    hm,dm=matter.mass_matrices(phi);F=float(matter.original.F(phi))
    hn=hm*np.sqrt(F);dn=dm*np.sqrt(F)
    ds=np.zeros_like(dm)
    ds[30:32,30:32]=matter.Y['s']*np.array([[0,1],[-1,0]])
    if ids is not None:
        ix=np.ix_(ids,ids);hm,dm,hn,dn,ds=[a[ix] for a in (hm,dm,hn,dn,ds)]
    # Fixed h derivative retains F and every Dirac block.
    dhds=hn*s/(6*F**1.5)
    ddds=dn*s/(6*F**1.5)+ds/np.sqrt(F)
    return (hm,dm),(dhds,ddds)


def data(B,G):
    e,v=np.linalg.eigh(B);neg=e<0;assert min(abs(e))>1e-9
    P=(v[:,neg])@v[:,neg].conj().T
    gg=v.conj().T@G@v
    transition=gg[np.ix_(~neg,neg)]
    de=e[~neg,None]-e[None,neg]
    mean=float(np.trace(P@G).real/2)
    noise=float(np.sum(abs(transition)**2)/2)
    metric=float(np.sum(abs(transition/de)**2)/2)
    return dict(E=float(e[neg].sum()/2),P=P,mean=mean,noise=noise,metric=metric,
                gap=float(min(abs(e))),excitation_min=float(de.min()),
                excitation_max=float(de.max()),e=e,v=v)


def neutral(h,s):
    pair,dp=radial(h,s,IDS)
    H=response.fock(pair);G=response.fock(dp)
    e,v=np.linalg.eigh(H);chi=v[:,0];g=v.conj().T@G@chi
    metric=float(np.sum(abs(g[1:]/(e[1:]-e[0]))**2))
    noise=float(np.sum(abs(g[1:])**2))
    return dict(H=H,G=G,e=e,chi=chi,E=float(e[0]),
                mean=float(g[0].real),metric=metric,noise=noise)


def aligned_vector(h,s,reference):
    chi=neutral(h,s)['chi'];z=np.vdot(reference,chi)
    return chi*np.exp(-1j*np.angle(z))


def local_mass_check():
    h,s=CENTER;pair,dp=radial(h,s);B=bdg(*pair);G=bdg(*dp)
    d=data(B,G);small=neutral(h,s)
    F=float(matter.original.F(np.array([0,h,0,0,s])))
    masses={k:abs(matter.Y[k])*h/np.sqrt(F) for k in ('u','d','e','nu')}
    majorana=abs(matter.Y['s'])*s/np.sqrt(F)
    exact=-6*(masses['u']+masses['d'])-2*masses['e']-np.sqrt(majorana**2+4*masses['nu']**2)
    assert abs(d['E']-exact)<1e-12
    nb=data(bdg(*radial(h,s,IDS)[0]),bdg(*radial(h,s,IDS)[1]))
    errors=[abs(small[key]-nb[key]) for key in ('E','mean','noise','metric')]
    # Charged masses change only eigenvalues on this radial slice.
    errors.extend([abs(d['noise']-small['noise']),abs(d['metric']-small['metric'])])
    step=2e-5
    pm=[data(bdg(*radial(h,s+sign*step)[0]),G) for sign in (1,-1)]
    errors.append(abs((pm[0]['E']-pm[1]['E'])/(2*step)-d['mean']))
    dP=(pm[0]['P']-pm[1]['P'])/(2*step)
    errors.append(abs(float(np.trace(dP@dP).real/4)-d['metric']))
    # Direct Fock eigenvector derivatives and original curved radial kinetic form.
    dh=(aligned_vector(h+step,s,small['chi'])-aligned_vector(h-step,s,small['chi']))/(2*step)
    ds=(aligned_vector(h,s+step,small['chi'])-aligned_vector(h,s-step,small['chi']))/(2*step)
    derivatives=np.array([dh,ds])
    connections=derivatives@small['chi'].conj()
    perpendicular=derivatives-connections[:,None]*small['chi']
    metric=np.real(perpendicular.conj()@perpendicular.T)
    q=np.array([h,s]);K=F*(np.eye(2)-np.outer(q,q)/(6*matter.original.M))
    w=old.source()['w'];a=K/w
    coefficient=float(np.sum(a*metric)/2)
    psi=.73+.21j;grad=np.array([.4-.3j,-.2+.1j])
    direct=grad[:,None]*small['chi']+psi*derivatives
    covariant=grad+psi*connections
    lhs=float(np.einsum('ij,ia,ja->',a,direct.conj(),direct).real)
    rhs=float((np.vdot(covariant,a@covariant)+abs(psi)**2*np.sum(a*metric)).real)
    errors.extend([abs(lhs-rhs),abs(metric[1,1]-small['metric'])])
    assert max(errors)<2e-8
    assert d['noise']>1e-5 and coefficient>0 and d['E']<0
    assert d['excitation_min']**2*d['metric']<=d['noise']+1e-12
    assert d['noise']<=d['excitation_max']**2*d['metric']+1e-12
    return dict(original_center=CENTER.tolist(),CAR_modes=32,
                full_ground_energy=d['E'],exact_mass_formula=exact,
                singlet_source_mean=d['mean'],singlet_source_variance=d['noise'],
                ground_metric_ss=d['metric'],neutral_Fock_dimension=16,
                maximum_independent_identity_error=max(errors),
                radial_Born_Huang_coefficient_without_hbar_squared=coefficient,
                radial_kinetic_product_identity_error=abs(lhs-rhs),
                quantum_metric_radial=metric.tolist(),
                original_bosonic_source_not_held_fixed_by_reference_change=True)


def star(gamma=0.,gauge_matrices=None):
    src=old.source();N=src['N'];original=old.wall.old.fields(N)
    index=[(0,N//4,N//8),(1,N//4,N//8),(0,N//4+1,N//8),(0,N//4,N//8+1)]
    phis=[original['q']['phi'][p] for p in index]
    if gauge_matrices is None:gauge_matrices=[np.eye(32)]*4
    hs=[];ds=[]
    for phi,R in zip(phis,gauge_matrices):
        h,d=matter.mass_matrices(phi);hs.append(R@h@R.conj().T);ds.append(R@d@R.T)
    h=block(hs);d=block(ds);hop=np.zeros_like(h)
    for i,alpha in enumerate(chiral.kinetic_matrices()):
        J=-1j*np.exp(-gamma)*alpha[:32,:32]/(2*src['eps'])
        J=gauge_matrices[0]@J@gauge_matrices[i+1].conj().T
        sl=slice(32*(i+1),32*(i+2));hop[:32,sl]=J;hop[sl,:32]=J.conj().T
    return bdg(h+hop,d),bdg(-hop,np.zeros_like(hop))


def propagation_source_check():
    B,G=star();r=data(B,G);step=2e-5
    up=data(star(step)[0],G);dn=data(star(-step)[0],G)
    dP=(up['P']-dn['P'])/(2*step)
    errors=[abs((up['E']-dn['E'])/(2*step)-r['mean']),
            abs(float(np.trace(dP@dP).real/4)-r['metric']),
            float(np.linalg.norm(r['P']@r['P']-r['P'])),
            float(np.linalg.norm(B@r['P']-r['P']@B))]
    energy_noise=data(B,B)['noise']
    rng=np.random.default_rng(72701);Rs=[]
    for _ in range(4):
        Rs.append(matter.representation(matter.gauge.group_exp(rng.normal(size=8),3),
            matter.gauge.group_exp(rng.normal(size=3),2),np.exp(1j*rng.normal())))
    R=block(Rs);RN=block([R,R.conj()]);Bg,Gg=star(gauge_matrices=Rs);rg=data(Bg,Gg)
    errors.extend([float(np.max(abs(Bg-RN@B@RN.conj().T))),
                   float(np.max(abs(rg['P']-RN@r['P']@RN.conj().T))),
                   abs(rg['E']-r['E']),abs(rg['noise']-r['noise'])])
    assert max(errors)<3e-7 and energy_noise<1e-22 and r['noise']>1e-4
    assert r['excitation_min']**2*r['metric']<=r['noise']+1e-11
    assert r['noise']<=r['excitation_max']**2*r['metric']+1e-11
    return dict(nodes=4,physical_modes=128,Nambu_dimension=256,
                all_three_original_Weyl_edge_directions=True,
                nonuniform_original651_node_fields=True,ground_energy=r['E'],
                quasiparticle_gap=r['gap'],geometry_source_mean=r['mean'],
                geometry_source_variance=r['noise'],ground_metric_geometry=r['metric'],
                ground_energy_variance=energy_noise,max_covariance_and_derivative_error=max(errors),
                covariance_not_a_standalone_Gauss_state_or_global_gap_proof=True)


def recorded_fiber_check():
    h,s=CENTER;center=neutral(h,s);src=old.source();lam=src['lam']
    x,weights=np.polynomial.hermite.hermgauss(40);rows=[]
    for hbar in (1e-4,2.5e-5,6.25e-6):
        delta=np.sqrt(hbar)*x;keep=abs(delta)<.18
        values=s+delta[keep];dx=delta[keep]
        F=matter.original.M-(h*h+values**2)/6
        cutoff=np.exp(-dx**2/(.18**2-dx**2))
        prob=weights[keep]*h**3*np.sqrt(matter.original.M)/F**3*cutoff**2;prob/=prob.sum()
        states=[neutral(h,z) for z in values]
        ene=np.array([z['E'] for z in states])
        mean=np.array([z['mean'] for z in states]);noise=np.array([z['noise'] for z in states])
        sigma=hbar**.4;f=values-lam*h;f0=s-lam*h
        # Actual y outcome kernel on a common grid; no latent phase label.
        width=10*(sigma+np.sqrt(hbar));y=np.linspace(f0-width,f0+width,801)
        kernel=np.exp(-(y[:,None]-f)**2/(2*sigma*sigma))/(np.sqrt(2*np.pi)*sigma)
        density=kernel@prob
        e_res=kernel@(prob*(ene-center['E'])**2)
        g_res=kernel@(prob*(noise+(mean-center['mean'])**2))
        normalization=float(np.trapezoid(density,y))
        energy_residual=float(np.trapezoid(e_res,y))
        source_residual=float(np.trapezoid(g_res,y))
        position=float(np.trapezoid(density*(y-f0)**2,y))
        assert abs(normalization-1)<1e-11
        assert abs(energy_residual-prob@((ene-center['E'])**2))<1e-12
        assert abs(source_residual-prob@(noise+(mean-center['mean'])**2))<1e-12
        assert abs(position-(sigma*sigma+prob@(dx*dx)))<1e-11
        rows.append(dict(hbar=hbar,sigma=sigma,normalization=normalization,
                         actual_recorded_mass_square_residual=energy_residual,
                         actual_recorded_source_square_residual=source_residual,
                         position_square_error=position,
                         conditional_ground_source_floor=center['noise']))
    assert rows[-1]['actual_recorded_mass_square_residual']<rows[0]['actual_recorded_mass_square_residual']/10
    assert abs(rows[-1]['actual_recorded_source_square_residual']-center['noise'])<1e-5
    return dict(rows=rows,neutral_conditional_fiber_fixture_only=True,
                original_actual_Gaussian_y_records=True,
                full_bosonic_derivative_and_Gauss_packet_extension_is_analytic=True)


def run():
    results=dict(original_mass_ground=local_mass_check(),
                 retained_spatial_propagation=propagation_source_check(),
                 actual_recorded_reference=recorded_fiber_check())
    deps=('research_note_558.md','research_note_574.md','research_note_591.md',
          'research_note_602.md','research_note_603.md','research_note_604.md',
          'research_note_630.md','research_note_632.md','research_note_633.md',
          'research_note_720.md','research_note_726.md','joint_fermion_gauss_completion.py',
          'joint_quantum_response_matching.py','joint_chiral_source_matching.py',
          'joint_recorded_classical_wall.py')
    return dict(round=727,tests_run=3,failures=0,errors=0,results=results,
                dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope='Original-matter conditional isolated ground band gives a new energy surface and geometric kinetic terms. Fixed-graph actual recorded semiclassical energy may concentrate about Hb+e0; source variance need not vanish. Original three-direction nonuniform Weyl star retains all128 physical modes, gauge covariance and source noise. No full interacting Gauss ground calculation, universal band gap, autonomous preparation, spatial continuum, self-consistent original Einstein source or adiabatic time error is claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(r,ensure_ascii=False,indent=2)+'\n')
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r['results'],ensure_ascii=False,indent=2))
