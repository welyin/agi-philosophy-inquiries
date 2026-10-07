"""803: original full fermion density vs its common-Cauchy-space generator.

A homogeneous diagnostic background, not the solved 753 spacetime, tests the
metric/frame/inner-product dictionary with all 64 original Nambu coefficients.
The continuum argument and original response scope are in research_note_803.
"""
from pathlib import Path
import argparse, json, sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'802'))
import original_bff_vertex as old
import quasifree_pairing_probe as pairing
TARGET=HERE/'cauchy_vertex_transport_results.json'

E0=np.array([[1.,.18,0.,0.],[0.,1.02,.1,0.],[.09,0.,1.01,0.],[0.,0.,0.,.97]])
G0=E0.T@old.ETA@E0
GI0=np.linalg.inv(G0)
K=np.array([[.13,.08,-.12,.03],[.08,-.11,.04,.07],[-.12,.04,.09,-.06],[.03,.07,-.06,.05]])
PHI=np.array([.07,.55,-.04,.09,.4])
DPHI=np.array([.03,-.05,.07,.02,.08])
MOMENTUM=np.array([.37,-.29,.21])
RNG=np.random.default_rng(830)
GAUGE=old.gauge_connection(RNG)
DGAUGE=old.gauge_connection(RNG)/3


def positive_root_and_derivative(m,dm):
    values,vectors=np.linalg.eigh(m)
    assert values.min()>.5
    roots=np.sqrt(values)
    t=(vectors*roots)@vectors.conj().T
    r=(vectors/roots)@vectors.conj().T
    dt=vectors@((vectors.conj().T@dm@vectors)/(roots[:,None]+roots[None,:]))@vectors.conj().T
    return t,r,dt


def data(lam,f,df):
    """Background parameter lam, local bump value f and time derivative df."""
    g=G0+lam*f*K
    dg=np.zeros((4,4,4));dg[0]=lam*df*K
    b=old.root_near_one(GI0@g)
    db=old.sylvester(b,GI0@dg[0])
    frame=E0@b;dframe=E0@db
    inv=np.linalg.inv(frame);dinv=-inv@dframe@inv
    volume=abs(np.linalg.det(frame));dvolume=volume*np.trace(inv@dframe)
    c=np.array([sum((inv[mu,a]*old.ALPHA[a] for a in range(4)),np.zeros((64,64),complex))
                for mu in range(4)])
    dc=np.array([sum((dinv[mu,a]*old.ALPHA[a] for a in range(4)),np.zeros((64,64),complex))
                 for mu in range(4)])
    christoffel=old.christoffel(g,dg)
    omega=np.array([frame@christoffel[mu]@inv-(dframe@inv if mu==0 else 0)
                    for mu in range(4)])
    spin=np.array([old.spin_lift(w) for w in omega])
    connection=spin+GAUGE+lam*f*DGAUGE
    z=sum((.5j*(c[mu]@connection[mu]-connection[mu].conj().T@c[mu])
           for mu in range(4)),np.zeros((64,64),complex))-old.mass(PHI+lam*f*DPHI)
    m=volume*c[0];dm=dvolume*c[0]+volume*dc[0]
    t,r,dt=positive_root_and_derivative(m,dm)
    norm_term=.5j*(dt@r-r@dt)
    h=sum((MOMENTUM[i]*r@(volume*c[i+1])@r for i in range(3)),np.zeros((64,64),complex))
    h-=r@(volume*z)@r
    h+=norm_term
    assert np.max(abs(h-h.conj().T))<2e-12
    return dict(v=volume,c=c,z=z,m=m,t=t,r=r,dt=dt,h=h,norm_term=norm_term)


def on_shell_density(candidate,base):
    """Coefficient of Psi^dagger ... Psi, before the common Nambu 1/2."""
    r=base['r'];d=-1j*r@base['h']
    matrix=.5j*(r@(candidate['v']*candidate['c'][0])@d-d.conj().T@(candidate['v']*candidate['c'][0])@r)
    matrix-=sum((MOMENTUM[i]*r@(candidate['v']*candidate['c'][i+1])@r for i in range(3)),np.zeros((64,64),complex))
    matrix+=r@(candidate['v']*candidate['z'])@r
    return matrix


def bump(t):
    if t<=0 or t>=1:return 0.,0.
    f=np.exp(4-1/(t*(1-t)))
    return f,f*(1-2*t)/(t*t*(1-t)**2)


def run():
    base=data(0.,0.,0.)
    eps=2e-5
    plus,minus=data(eps,1.,0.),data(-eps,1.,0.)
    h1=(plus['h']-minus['h'])/(2*eps)
    dt=(plus['t']-minus['t'])/(2*eps)
    anti=dt@base['r']-base['r']@dt
    time_plus,time_minus=data(eps,0.,1.),data(-eps,0.,1.)
    h2=(time_plus['h']-time_minus['h'])/(2*eps)
    l1=(on_shell_density(plus,base)-on_shell_density(minus,base))/(2*eps)
    l2=(on_shell_density(time_plus,base)-on_shell_density(time_minus,base))/(2*eps)
    h0=base['h']
    errors=dict(free_on_shell_density=float(np.max(abs(on_shell_density(base,base)))),
        time_generator_identity=float(np.max(abs(h2-.5j*anti))),
        action_generator_identity=float(np.max(abs(l1+h1+(h0@anti-anti@h0)/2))),
        time_spin_density_cancellation=float(np.max(abs(l2))))
    energy,vectors=np.linalg.eigh(h0)
    freq=energy[:,None]-energy[None,:]
    nodes,weights=np.polynomial.legendre.leggauss(160)
    fourier=np.zeros((64,64),complex);derivative=np.zeros_like(fourier)
    for x,weight in zip((nodes+1)/2,weights/2):
        f,df=bump(float(x));phase=np.exp(1j*freq*x)
        fourier+=weight*f*phase
        derivative+=weight*df*phase
    h1e=vectors.conj().T@h1@vectors;h2e=vectors.conj().T@h2@vectors
    l1e=vectors.conj().T@l1@vectors
    action_kernel=fourier*l1e
    evolution_kernel=-fourier*h1e-derivative*h2e
    errors['bump_integration_by_parts']=float(np.max(abs(derivative+1j*freq*fourier)))
    errors['integrated_original_action_vs_evolution']=float(np.max(abs(action_kernel-evolution_kernel)))
    omitted=float(np.linalg.norm(action_kernel+fourier*h1e))
    # Independently solve the linear response of the unitary evolution with the
    # exact Frechet derivative of each midpoint matrix exponential (730 tool).
    midpoint_errors=[]
    for steps in (48,96):
        u=np.eye(64,dtype=complex);du=np.zeros_like(u);step=1/steps
        for j in range(steps):
            f,df=bump((j+.5)*step)
            evolution,devolution=old.original.step(h0,f*h1+df*h2,step)
            du=devolution@u+evolution@du
            u=evolution@u
        kernel=-1j*u.conj().T@du
        kern_e=vectors.conj().T@kernel@vectors
        midpoint_errors.append(float(np.max(abs(kern_e-evolution_kernel))))
    assert max(errors.values())<4e-9,errors
    assert omitted>1e-4 and np.linalg.norm(h2)>1e-4
    assert midpoint_errors[1]<midpoint_errors[0]/3 and midpoint_errors[1]<1e-5,midpoint_errors
    q=pairing.run()
    assert q==json.loads(pairing.TARGET.read_text('utf-8'))
    return dict(round=803,all_checks_passed=True,fresh_test_groups=2,
        original_nambu_dimension=64,
        cauchy_density_transport=dict(maximum_residuals=errors,
            norm_connection_size=float(np.linalg.norm(h2)),
            freeze_time_inner_product_integrated_defect=omitted,
            midpoint_steps=[48,96],midpoint_kernel_errors=midpoint_errors,
            all_original_mass_parameters_retained=True,
            diagnostic_is_homogeneous_background_not_753=True),
        ordered_covariance_and_selection=q,
        original_continuum_vertex_cauchy_mapping_proven_analytically=True,
        continuous_original_W_response_value_computed=False,
        original_spacetime_or_apparatus_simulated=False,
        finite_coupling_or_graph_continuum_proven=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))
