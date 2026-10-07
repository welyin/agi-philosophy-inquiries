"""802: first variation of the original full fermion quadratic density.

Original 32 CAR / 64 Nambu matrices and Yukawa constants are imported read-only.
Local jets below test coefficient identities, not the original continuum state.
No new spacetime solution, finite-coupling signal, or apparatus is simulated.
"""
from pathlib import Path
import argparse
import json
import sys
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from research_layout import Layout, ResearchRuntime
with ResearchRuntime(Layout()).installed():
    import joint_dynamic_continuum_reference as original

I4 = np.eye(4)
ETA = np.diag([-1., 1., 1., 1.])
I64 = np.eye(64, dtype=complex)
ALPHA = [I64, *original.GAMMA]
RESULT = HERE/'original_bff_vertex_results.json'


def mass(phi):
    return original.old.bdg(*original.old.matter.mass_matrices(phi))


# Numerator is linear in the five original REAL Higgs/singlet coordinates.
NUMERATORS = [mass(np.eye(5)[a])*np.sqrt(11/6) for a in range(5)]
star = mass(np.array([0., .5, 0., 0., .4]))
eigen, eigvec = np.linalg.eigh(star)
BETA = (eigvec*np.sign(eigen))@eigvec.conj().T
GAMMA = [1j*BETA@a for a in ALPHA]


def dmass(phi, dphi):
    f = original.old.matter.original.F(phi)
    dx = dphi/np.sqrt(f)+phi*np.dot(phi, dphi)/(6*f**1.5)
    return sum((n*x for n, x in zip(NUMERATORS, dx)), np.zeros((64,64),complex))


def root_near_one(matrix):
    x = matrix-I4
    assert np.linalg.norm(x, 2) < .7
    value, power, coeff = I4.copy(), I4.copy(), 1.
    for n in range(1, 200):
        power = power@x
        coeff *= (1.5-n)/n
        term = coeff*power
        value += term
        if np.linalg.norm(term) < 1e-18:
            return value
    raise AssertionError('square root failed')


def sylvester(b, rhs):
    op = np.kron(I4,b)+np.kron(b.T,I4)
    return np.linalg.solve(op, rhs.reshape(-1,order='F')).reshape((4,4),order='F')


def christoffel(g, dg):
    inv = np.linalg.inv(g)
    return np.array([[[sum(inv[r,s]*(dg[mu,s,nu]+dg[nu,s,mu]-dg[s,mu,nu])
                           for s in range(4))/2 for nu in range(4)]
                      for r in range(4)] for mu in range(4)])


def spin_lift(omega):
    lower = ETA@omega
    return sum((lower[a,b]*GAMMA[a]@GAMMA[b]/4
                for a in range(4) for b in range(4)), np.zeros((64,64),complex))


def frame_data(e0, h, dh, t):
    # E0 is locally constant in this diagnostic. The analytic note is general.
    g0 = e0.T@ETA@e0
    gi0 = np.linalg.inv(g0)
    g, dg = g0+t*h, t*dh
    b = root_near_one(gi0@g)
    frame = e0@b
    dframe = np.array([e0@sylvester(b,gi0@dg[mu]) for mu in range(4)])
    inv = np.linalg.inv(frame)
    connection = christoffel(g,dg)
    omega = np.array([frame@connection[mu]@inv-dframe[mu]@inv for mu in range(4)])
    c = np.array([sum((inv[mu,a]*ALPHA[a] for a in range(4)),np.zeros((64,64),complex))
                  for mu in range(4)])
    return abs(np.linalg.det(frame)), c, np.array([spin_lift(w) for w in omega])


def frame_tangent(e0, h, dh):
    g0 = e0.T@ETA@e0
    gi0, inv0 = np.linalg.inv(g0), np.linalg.inv(e0)
    de = e0@gi0@h/2
    di = -inv0@de@inv0
    dc = np.array([sum((di[mu,a]*ALPHA[a] for a in range(4)),np.zeros((64,64),complex))
                   for mu in range(4)])
    # Flat background jet; delta Gamma is linear in the supplied dh.
    dgamma = christoffel(g0,dh)
    domega = np.array([e0@dgamma[mu]@inv0-e0@gi0@dh[mu]@inv0/2 for mu in range(4)])
    return np.trace(gi0@h)/2, dc, np.array([spin_lift(w) for w in domega])


def kinetic(c, omega, left, dleft, right, dright):
    # Polarized Hermitian Weyl/Nambu density, without the overall Nambu 1/2.
    value = 0j
    for mu in range(4):
        value += .5j*(left.conj()@c[mu]@dright[mu]-dleft[mu].conj()@c[mu]@right)
        value += .5j*left.conj()@(c[mu]@omega[mu]-omega[mu].conj().T@c[mu])@right
    return value


def bilinear(c, omega, b, left, dleft, right, dright):
    return kinetic(c,omega,left,dleft,right,dright)-left.conj()@b@right


def gauge_connection(rng):
    # Original electroweak representation, including all charged/sterile blocks.
    generators = original.gauge_h(rng.normal(size=(3,3))*.1, rng.normal(size=3)*.1)
    generators.append(sum(generators)/5)
    return np.array([np.block([[1j*g,np.zeros_like(g)],
                               [np.zeros_like(g),-1j*g.conj()]]) for g in generators])


def coefficient_check():
    rng = np.random.default_rng(802)
    errors = dict(first_variation=0., mass_derivative=0., spin_connection_contraction=0., clifford=0.)
    missing_density, missing_metric, missing_mass = [], [], []
    for a in range(4):
        for b in range(4):
            errors['clifford'] = max(errors['clifford'],float(np.linalg.norm(
                GAMMA[a]@GAMMA[b]+GAMMA[b]@GAMMA[a]-2*ETA[a,b]*I64)))
    for _ in range(12):
        e0 = I4+rng.normal(size=(4,4))*.035
        h = rng.normal(size=(4,4))*.12; h += h.T
        dh = rng.normal(size=(4,4,4))*.08; dh += dh.transpose(0,2,1)
        phi = rng.normal(size=5)*.3
        dphi = rng.normal(size=5)*.2
        u = (rng.normal(size=64)+1j*rng.normal(size=64))/8
        v = (rng.normal(size=64)+1j*rng.normal(size=64))/8
        du = (rng.normal(size=(4,64))+1j*rng.normal(size=(4,64)))/8
        dv = (rng.normal(size=(4,64))+1j*rng.normal(size=(4,64)))/8
        connection, da = gauge_connection(rng), gauge_connection(rng)
        volume,c,spin = frame_data(e0,h,dh,0.)
        nu,dc,ds = frame_tangent(e0,h,dh)
        b,db = mass(phi),dmass(phi,dphi)
        density = nu*bilinear(c,connection,b,u,du,v,dv)
        metric = kinetic(dc,connection,u,du,v,dv)
        con = sum((.5j*u.conj()@(c[mu]@(ds[mu]+da[mu])-
                  (ds[mu]+da[mu]).conj().T@c[mu])@v for mu in range(4)),0j)
        mass_part = -u.conj()@db@v
        analytic = .5*volume*(density+metric+con+mass_part)
        def evaluate(t):
            vol, cc, ss = frame_data(e0,h,dh,t)
            return .5*vol*bilinear(cc,ss+connection+t*da,mass(phi+t*dphi),u,du,v,dv)
        step = 2e-5
        numeric = (evaluate(step)-evaluate(-step))/(2*step)
        errors['first_variation'] = max(errors['first_variation'],float(abs(numeric-analytic)))
        errors['mass_derivative'] = max(errors['mass_derivative'],float(np.max(abs(
            (mass(phi+step*dphi)-mass(phi-step*dphi))/(2*step)-db))))
        spin_term = sum((c[mu]@ds[mu]-ds[mu].conj().T@c[mu] for mu in range(4)),np.zeros((64,64),complex))
        errors['spin_connection_contraction'] = max(errors['spin_connection_contraction'],float(np.max(abs(spin_term))))
        missing_density.append(float(abs(.5*volume*density)))
        missing_metric.append(float(abs(.5*volume*metric)))
        missing_mass.append(float(abs(.5*volume*mass_part)))
    assert max(errors.values()) < 2e-8,errors
    assert min(missing_density)>1e-6 and min(missing_metric)>1e-5 and min(missing_mass)>1e-5
    return dict(samples=12,original_nambu_dimension=64,maximum_residuals=errors,
        omitted_density_minimum_error=min(missing_density),
        omitted_metric_minimum_error=min(missing_metric),
        omitted_mass_minimum_error=min(missing_mass),
        coefficients_not_continuum_state=True)


def on_shell_density_and_shear():
    rng = np.random.default_rng(803)
    phi = np.array([.13,.55,-.09,.04,.4])
    connection = gauge_connection(rng)
    b = mass(phi)
    z = sum((.5j*(ALPHA[mu]@connection[mu]-connection[mu].conj().T@ALPHA[mu])
             for mu in range(4)),np.zeros((64,64),complex))-b
    def jets(u,spatial):
        dt = -sum((ALPHA[i+1]@spatial[i] for i in range(3)),np.zeros(64,complex))+1j*z@u
        return np.array([dt,*spatial])
    maximum = 0.
    for _ in range(10):
        u = (rng.normal(size=64)+1j*rng.normal(size=64))/8
        v = (rng.normal(size=64)+1j*rng.normal(size=64))/8
        du = jets(u,(rng.normal(size=(3,64))+1j*rng.normal(size=(3,64)))/8)
        dv = jets(v,(rng.normal(size=(3,64))+1j*rng.normal(size=(3,64)))/8)
        maximum=max(maximum,float(abs(bilinear(np.array(ALPHA),connection,b,u,du,v,dv))))
    # Physical metric shear: dh || dt, ds in span(dt,dx) in an adapted frame.
    h=np.zeros((4,4));h[2,3]=h[3,2]=.2
    nu,dc,ds=frame_tangent(I4,h,np.zeros((4,4,4)))
    # +1 eigenvector of the original Gamma_z, so the k slope is fixed analytically.
    vals,vecs=np.linalg.eigh(ALPHA[3]);u=vecs[:,np.argmax(vals)]
    def value(k):
        spatial=np.zeros((3,64),complex);spatial[1]=1j*k*u
        du=jets(u,spatial)
        return float((.5*kinetic(dc,connection,u,du,u,du)).real)
    observed_slope=value(2)-value(1)
    assert maximum<2e-12 and abs(observed_slope-.05)<2e-12
    assert nu==0 and max(np.max(abs(s)) for s in ds)==0
    # This validates the coefficient on fermion equation jets. The metric shear
    # is not asserted to solve the coupled free boson equations.
    return dict(polarized_on_shell_samples=10,density_bilinear_maximum=maximum,
        shear_amplitude=.2,with_original_nambu_half_weight_expected_k_slope=.05,
        measured_k_slope=observed_slope,
        boson_direction_satisfies_linear_reference_slice=True,
        boson_direction_claimed_on_shell=False,
        nonzero_vertex_is_nonzero_record_signal=False)


def run():
    return dict(working_round=802,all_checks_passed=True,
        full_original_coefficient=coefficient_check(),
        original_fermion_equation_cancellation=on_shell_density_and_shear(),
        original_yukawa_parameters={k:[v.real,v.imag] for k,v in original.old.matter.Y.items()},
        same_continuum_W_evaluated=False,complete_joint_response_nonzero_proven=False,
        original_background_solved_here=False,finite_coupling_or_apparatus_proven=False,
        new_numbered_test_groups=0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    result=run()
    if args.write:
        assert not RESULT.exists(),'Do not overwrite saved evidence.'
        RESULT.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:
        saved=json.loads(RESULT.read_text(encoding='utf-8'))
        assert result==saved
    print(json.dumps(result,ensure_ascii=False,indent=2))
