"""587: original curved matter and quotient-gauge energy currents.

Operator claims are on the invariant smooth compact core. The continuum
diagnostic is a classical initial-data symbol comparison, not quantum GR.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_record_source_compression as records
import joint_quotient_gauge_completion as gauge

HERE = Path(__file__).resolve().parent
TARGET = HERE/'joint_matter_energy_current_results.json'
PAR = gauge.parameters()


def dag(a):
    return np.swapaxes(a.conj(), -1, -2)


def trace(a):
    return np.trace(a, axis1=-2, axis2=-1)


def edge_gradients(x, y):
    """Covector derivatives of half the exact target distance squared."""
    fx, fy = original.F(x), original.F(y)
    root = np.sqrt(fx*fy)
    z = (np.sum((x-y)**2, axis=-1)/12+(np.sqrt(fx)-np.sqrt(fy))**2/2)/root
    c = 1+z
    ratio = np.ones_like(c)
    np.divide(2*np.arcsinh(np.sqrt(z/2)), np.sqrt(z*(z+2)), out=ratio, where=z>1e-14)
    gx = ratio[..., None]*(c[..., None]*x/fx[..., None]-y/root[..., None])
    gy = ratio[..., None]*(c[..., None]*y/fy[..., None]-x/root[..., None])
    return gx, gy


def magnetic_direction(C, W, z, dC, dW, dz):
    """Derivative of the FULL faithful 569 quotient potential, not a trace replacement."""
    tc, tw, dc, dw = trace(C), trace(W), trace(dC), trace(dW)
    dchi = PAR['wq']*(dc*tw*z+tc*dw*z+tc*tw*dz+dc.conj()*(z**-4+z**2)
                         +tc.conj()*(-4*z**-5+2*z)*dz)
    dchi += PAR['wl']*(dw*z**-3-3*tw*z**-4*dz+6*z**5*dz)
    delta, coeff = gauge.coefficients(PAR)
    return -delta*dchi.real-2*coeff[0]*(tc.conj()*dc).real-2*coeff[1]*(tw.conj()*dw).real-coeff[2]*(6*z**5*dz).real


def compose_face(links, variations):
    C = links[0][0]@links[1][0]@dag(links[2][0])@dag(links[3][0])
    W = links[0][1]@links[1][1]@dag(links[2][1])@dag(links[3][1])
    z = links[0][2]*links[1][2]/links[2][2]/links[3][2]
    answers = []
    for k in range(4):
        dc, dw = [], []
        for j in range(4):
            cc, ww = variations[j][:2] if j == k else links[j][:2]
            dc.append(dag(cc) if j >= 2 else cc)
            dw.append(dag(ww) if j >= 2 else ww)
        dC, dW = dc[0]@dc[1]@dc[2]@dc[3], dw[0]@dw[1]@dw[2]@dw[3]
        dz = z*variations[k][2]/links[k][2]*(1 if k < 2 else -1)
        answers.append(magnetic_direction(C, W, z, dC, dW, dz))
    return answers


def incidence_check():
    path = HERE/'round587_drafts/local_current_entry.py'
    spec = importlib.util.spec_from_file_location('frozen_current_entry', path)
    entry = importlib.util.module_from_spec(spec); spec.loader.exec_module(entry)
    got = entry.run()
    assert got == json.loads(path.with_name('local_current_entry_results.json').read_text('utf8'))
    return got


def curved_derivative_check():
    rng = np.random.default_rng(5871)
    x, y, vx, vy = rng.normal(size=(4, 5))*.3
    gx, gy = edge_gradients(x, y)
    expected = float(gx@vx+gy@vy)
    rows = []
    for step in (.02, .01, .005):
        got = float((original.distance_squared(x+step*vx, y+step*vy)-original.distance_squared(x-step*vx, y-step*vy))/(4*step))
        rows.append(dict(step=step, error=abs(got-expected)))
    assert rows[-1]['error'] < rows[0]['error']/12
    U = gauge.group_exp(rng.normal(size=3), 2); z = np.exp(.37j)
    def rotate(q):
        w = z**3*(U@(q[:2]+1j*q[2:4]))
        return np.r_[w.real, w.imag, q[4]]
    a, b = edge_gradients(rotate(x), rotate(y))
    error = abs(float(a@rotate(vx)+b@rotate(vy))-expected)
    assert error < 1e-14
    return dict(directional_derivative=expected, finite_difference=rows, gauge_covariance_error=error)


def full_magnetic_check():
    rng = np.random.default_rng(5872)
    links, velocities, generators = [], [], []
    for k in range(4):
        C, W, z = gauge.group_exp(.3*rng.normal(size=8), 3), gauge.group_exp(.4*rng.normal(size=3), 2), np.exp(.13j*rng.normal())
        A = np.einsum('a,aij->ij', .4*rng.normal(size=8), gauge.generators(3))
        B = np.einsum('a,aij->ij', .4*rng.normal(size=3), gauge.generators(2))
        a = float(.2*rng.normal())
        links.append((C, W, z)); generators.append((A, B, a))
        velocities.append((1j*A@C, 1j*B@W, 1j*a*z))
    direct = np.array(compose_face(links, velocities))
    def expm(A, t):
        eig, U = np.linalg.eigh(A)
        return (U*np.exp(1j*t*eig))@U.conj().T
    def value(q):
        C = q[0][0]@q[1][0]@dag(q[2][0])@dag(q[3][0])
        W = q[0][1]@q[1][1]@dag(q[2][1])@dag(q[3][1])
        z = q[0][2]*q[1][2]/q[2][2]/q[3][2]
        return gauge.potential(C, W, z, PAR)
    rows = []
    for step in (.002, .001, .0005):
        numeric = []
        for k in range(4):
            vals = []
            for sign in (1, -1):
                q = list(links); C, W, z = links[k]; A, B, a = generators[k]
                q[k] = (expm(A, sign*step)@C, expm(B, sign*step)@W, np.exp(1j*sign*step*a)*z)
                vals.append(value(q))
            numeric.append((vals[0]-vals[1])/(2*step))
        rows.append(dict(step=step, error=float(np.max(abs(numeric-direct)))))
    assert rows[-1]['error'] < 2e-5 and rows[-1]['error'] < rows[0]['error']/10
    # Independent representative changes for each edge of the quotient.
    lifted, lifted_v = [], []
    for k, (q, v) in enumerate(zip(links, velocities)):
        power = k+1; fac = (np.exp(2j*np.pi*power/3), (-1)**power, np.exp(1j*np.pi*power/3))
        lifted.append(tuple(a*b for a,b in zip(fac,q))); lifted_v.append(tuple(a*b for a,b in zip(fac,v)))
    lift_error = float(np.max(abs(np.array(compose_face(lifted, lifted_v))-direct)))
    # Actual node gauge transformations of both links and their tangents.
    nodes = [(gauge.group_exp(rng.normal(size=8),3), gauge.group_exp(rng.normal(size=3),2), np.exp(1j*rng.normal())) for _ in range(4)]
    ends = ((0,1),(1,2),(3,2),(0,3)); transformed, transformed_v = [], []
    for (i,j), q,v in zip(ends,links,velocities):
        transformed.append((nodes[i][0]@q[0]@dag(nodes[j][0]),nodes[i][1]@q[1]@dag(nodes[j][1]),nodes[i][2]*q[2]/nodes[j][2]))
        transformed_v.append((nodes[i][0]@v[0]@dag(nodes[j][0]),nodes[i][1]@v[1]@dag(nodes[j][1]),nodes[i][2]*v[2]/nodes[j][2]))
    gauge_error = float(np.max(abs(np.array(compose_face(transformed, transformed_v))-direct)))
    assert max(lift_error,gauge_error) < 1e-11 and np.linalg.norm(direct)>1e-3
    return dict(full_colour_weak_circle_derivatives=direct.tolist(), nonzero_face_potential=value(links),
                finite_difference=rows, independent_lift_error=lift_error, node_gauge_error=gauge_error)


def physical_source_rates(N):
    q = original.shared_source(N); continuum = original.geometry.make_source(N)
    eps = q['eps']; x,y,z = np.moveaxis(q['grid'],-1,0)
    psi = 1.1+.05*np.cos(x)+.02*np.sin(y)
    tests = [np.sin(x),np.cos(x+y),np.sin(x)*np.cos(y)]
    gradients = [np.stack((np.cos(x),0*x,0*x),axis=-1),
                 np.stack((-np.sin(x+y),-np.sin(x+y),0*x),axis=-1),
                 np.stack((np.cos(x)*np.cos(y),-np.sin(x)*np.sin(y),0*x),axis=-1)]
    phi = q['phi']; velocity = np.einsum('...ij,...j->...i',original.inverse(phi),q['P'])/(eps**3*psi[...,None]**6)
    W = original.lattice.su2(q['a'],eps); circle = np.exp(1j*eps*q['a0'])
    links, rates = [], []
    for mu in range(3):
        pe = (psi+np.roll(psi,-1,axis=mu))/2
        A = np.einsum('...a,aij->...ij',2*PAR['b'][1]/eps*pe[...,None]**-2*q['phat'][...,mu,:],original.lattice.T)
        omega = 2*PAR['b'][2]/eps*pe**-2*q['p0hat'][...,mu]
        links.append((np.broadcast_to(np.eye(3,dtype=complex),phi.shape[:-1]+(3,3)),W[...,mu,:,:],circle[...,mu]))
        rates.append((np.zeros(phi.shape[:-1]+(3,3),complex),1j*A@W[...,mu,:,:],1j*omega*circle[...,mu]))
    scalar = np.zeros(len(tests)); magnetic = np.zeros(len(tests))
    def shifted(q,mu):return tuple(np.roll(a,-1,axis=mu) for a in q)
    for mu in range(3):
        pe=(psi+np.roll(psi,-1,axis=mu))/2
        def transform(vector):
            X=vector[...,:2]+1j*vector[...,2:4]
            X=circle[...,mu,None]**3*np.einsum('...ij,...j->...i',W[...,mu,:,:],np.roll(X,-1,axis=mu))
            return original.real_phi(X,np.roll(vector[...,4],-1,axis=mu))
        gx,gy=edge_gradients(phi,transform(phi))
        dw_i=eps*pe**2*np.sum(gx*velocity,axis=-1)
        dw_j=eps*pe**2*np.sum(gy*transform(velocity),axis=-1)
        for k,f in enumerate(tests):scalar[k]+=float(np.sum((np.roll(f,-1,axis=mu)-f)*(dw_i-dw_j)/2))
        for nu in range(mu+1,3):
            pface=(psi+np.roll(psi,-1,axis=mu)+np.roll(psi,-1,axis=nu)+np.roll(np.roll(psi,-1,axis=mu),-1,axis=nu))/4
            face=[links[mu],shifted(links[nu],mu),shifted(links[mu],nu),links[nu]]
            fv=[rates[mu],shifted(rates[nu],mu),shifted(rates[mu],nu),rates[nu]]
            dw=[a/(eps*pface**2) for a in compose_face(face,fv)]
            for k,f in enumerate(tests):
                fm,fn=np.roll(f,-1,axis=mu),np.roll(f,-1,axis=nu)
                fmn=np.roll(fm,-1,axis=nu); ff=(f+fm+fn+fmn)/4
                fe=((f+fm)/2,(fm+fmn)/2,(fn+fmn)/2,(f+fn)/2)
                magnetic[k]+=sum(float(np.sum((ff-e)*d)) for e,d in zip(fe,dw))
    smom=np.einsum('...a,...ia->...i',continuum['p'],continuum['Dphi'])
    gmom=continuum['mom']-smom
    expected_s=np.array([-eps**3*np.sum(psi**-4*np.sum(df*smom,axis=-1)) for df in gradients])
    expected_g=np.array([-eps**3*np.sum(psi**-4*np.sum(df*gmom,axis=-1)) for df in gradients])
    gw,g0=original.lattice.residual(q)
    return dict(N=N,eps=eps,scalar_rate=scalar.tolist(),magnetic_rate=magnetic.tolist(),
                scalar_continuum=expected_s.tolist(),magnetic_continuum=expected_g.tolist(),
                scalar_error=float(np.max(abs(scalar-expected_s))),magnetic_error=float(np.max(abs(magnetic-expected_g))),
                total_error=float(np.max(abs(scalar+magnetic-expected_s-expected_g))),
                exact_Gauss_density_residual=float(max(np.max(abs(gw)),np.max(abs(g0)))/eps**3))


def continuum_check():
    rows=[physical_source_rates(n) for n in (8,12,20,32)]
    for key in ('scalar_error','magnetic_error','total_error'):
        assert rows[-1][key]<rows[0][key]/7, (key,rows)
    assert max(r['exact_Gauss_density_residual'] for r in rows)<1e-10
    assert max(abs(x) for x in rows[-1]['magnetic_continuum'])>.01
    return dict(rows=rows,nonzero_EW_magnetic_source=True,colour_in_source_zero=True,
                positive_fixed_smooth_geometry_input=True,classical_symbol_only=True,
                empirical_rate_not_uniform_quantum_limit=True)


def record_current_noise_check():
    # A one-coordinate differential diagnostic for the actual edge principal
    # vector. Full-operator product rule is proved in the report.
    x=np.array([.14,.53,-.17,.1,.42]);y=np.array([.23,.41,-.09,.14,.31])
    hbar=.7; weight=.8; coupling=.9
    def v(s):
        q=x.copy();q[-1]=s;g,_=edge_gradients(q,y)
        return float(-coupling/(2*weight)*(original.inverse(q)@g)[-1])
    def wave(s):return np.exp(-.2*s*s+.17j*s)
    expected=hbar**2*v(x[-1])**2*float(records.delta_a(x[-1]))
    rows=[]
    for step in (.004,.002,.001):
        def J(fun,s):
            der=(fun(s+step)-fun(s-step))/(2*step)
            divergence=(v(s+step)-v(s-step))/(2*step)
            return -1j*hbar*(v(s)*der+.5*divergence*fun(s))
        def funcs(s):
            a,_,b,_=records.instruments(np.array(s));return np.r_[a,b]
        first=[];second=[]
        for k in range(6):
            fun=lambda s,k=k:funcs(s)[k]*wave(s)
            first.append(funcs(x[-1])[k]*J(fun,x[-1]))
            second.append(funcs(x[-1])[k]*J(lambda s:J(fun,s),x[-1]))
        delta=(sum(second[:4])-sum(second[4:]))/wave(x[-1])
        first_error=max(abs(sum(first[:4])-J(wave,x[-1])),abs(sum(first[4:])-J(wave,x[-1])))
        rows.append(dict(step=step,first_moment_error=float(first_error),second_moment_gap_error=float(abs(delta-expected))))
    assert expected>1e-5 and rows[-1]['second_moment_gap_error']<2e-8
    assert rows[-1]['second_moment_gap_error']<rows[0]['second_moment_gap_error']/10
    return dict(pointwise_noise_gap_coefficient=expected,finite_difference=rows,
                fixed_neighbour_local_differential_diagnostic=True,full_Gauss_expectation_identity_is_analytic=True)


def run():
    evidence=dict(full_allocation=incidence_check(),curved_edge=curved_derivative_check(),
                  full_quotient_magnetic_derivative=full_magnetic_check(),
                  same_Gauss_source_continuum_current=continuum_check(),
                  same_record_current_fluctuations=record_current_noise_check())
    deps=('joint_curved_quantum_source.py','joint_record_source_compression.py','joint_quotient_gauge_completion.py',
          'joint_gauss_continuum_sampling.py','joint_gauss_einstein_initial_data.py',
          'round587_drafts/STATUS.md','round587_drafts/local_current_derivation.md',
          'round587_drafts/local_current_entry.py','round587_drafts/local_current_entry_results.json')
    return dict(round=587,tests_run=5,failures=0,errors=0,evidence=evidence,
                dependency_hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in deps},
                scope='full finite-model core current; fixed-graph semiclassical then smooth classical source limit; no joint fixed-hbar continuum or quantum HDA')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    elif TARGET.exists():assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
