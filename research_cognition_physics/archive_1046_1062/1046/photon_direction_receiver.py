"""1046: fixed finite-band dipole receiver, normal photon packets.

Default: read-only replay, with explicit floating-point comparison tolerance.
--write: exclusive creation of the initial scientific result only.
The continuous proofs live in proof.md; finite quadrature is a calibration.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGET = HERE / 'photon_direction_receiver_results.json'
HISTORY = (
    'archive_370_428/research_note_381.md',
    'archive_370_428/research_note_382.md',
    'archive_370_428/research_note_383.md',
    'archive_370_428/research_note_384.md',
    'archive_370_428/research_note_386.md',
    'archive_370_428/research_note_425.md',
    'archive_467_530/research_note_522.md',
    'archive_467_530/research_note_523.md',
    'archive_956_989/research_note_970.md',
    'archive_1009_1043/research_note_1038.md',
    'archive_1009_1043/research_note_1041.md',
    'archive_1044_/research_note_1044.md',
    'archive_1044_/research_note_1045.md',
)
I2 = np.eye(2, dtype=complex)
SIGMA = np.array([[[0,1],[1,0]], [[0,-1j],[1j,0]], [[1,0],[0,-1]]], complex)
G = .1
OMEGA = 1.5
TIME = 1.
A1 = F(71,960)
A2 = F(31,672)
B1 = F(3,160)
B2 = F(13,1680)
J2 = F(1,20)


def pauli(v):
    return np.einsum('i,ijk->jk', v, SIGMA)


def psqrt(a):
    eig, basis = np.linalg.eigh((a+a.conj().T)/2)
    assert eig.min() > -1e-12
    return (basis*np.sqrt(np.maximum(eig,0))) @ basis.conj().T


def frame(n):
    n = np.asarray(n, float); n /= np.linalg.norm(n)
    seed = np.eye(3)[np.argmin(abs(n))]
    u = seed-n*np.dot(seed,n); u /= np.linalg.norm(u)
    v = np.cross(n,u)
    return np.column_stack([u,v,n])


def angular(n, nc=12, nphi=16):
    """Direct physical transverse-vector quadrature, not assumed Bloch effects."""
    rot = frame(n); n=rot[:,2]; eps=(rot[:,0]+1j*rot[:,1])/np.sqrt(2)
    cs=[]; ws=[]
    x,w=np.polynomial.legendre.leggauss(nc)
    for lo,hi in [(-1.,.5),(.5,1.)]:
        cs.extend((lo+hi)/2+(hi-lo)*x/2)
        ws.extend((hi-lo)*w/2)
    norm=0.; tj=np.zeros((2,2),complex); self_energy=np.zeros((2,2),complex)
    mean=np.zeros(3); raw=np.zeros((3,3)); gram=np.zeros((3,3))
    for c,wc in zip(cs,ws):
        aperture=max(2*c-1,0)**2
        for phi in 2*np.pi*np.arange(nphi)/nphi:
            local=np.array([np.sqrt(1-c*c)*np.cos(phi),np.sqrt(1-c*c)*np.sin(phi),c])
            khat=rot@local; P=np.eye(3)-np.outer(khat,khat)
            weight=wc/(2*nphi)
            photon=aperture*(P@eps)/np.sqrt(float(A2))
            probability=float(np.vdot(photon,photon).real)
            norm+=weight*probability
            mean+=weight*probability*khat
            tj+=weight*pauli(photon)
            raw+=weight*aperture*P
            gram+=weight*aperture**2*P
            for i in range(3):
                for j in range(3): self_energy+=weight*P[i,j]*(SIGMA[i]@SIGMA[j])
    return {'n':n,'eps':eps,'norm':norm,'tj':tj,'C':tj/np.sqrt(2),
            'AAstar':self_energy,'mean_direction':mean,'B':raw,'S2':gram}


def radial(n):
    x,w=np.polynomial.legendre.leggauss(n)
    k=1.5+x/2
    weight=(w/2)*(8/3)*np.sin(np.pi*(k-1))**4
    norm=float(weight.sum());weight/=norm
    return k,weight,norm


def unitary(H,t):
    eig,basis=np.linalg.eigh(H)
    return (basis*np.exp(-1j*t*eig))@basis.conj().T


def dynamics(n):
    k,w,raw_norm=radial(n)
    star=np.diag(np.r_[OMEGA,k]).astype(complex)
    star[0,1:]=star[1:,0]=G*np.sqrt(2*w)
    initial=np.r_[0.,np.sqrt(w)]
    out=unitary(star,TIME)@initial
    dark=np.sqrt(w)*np.exp(-1j*k*TIME)
    # Direct full bright/dark quantum Hamiltonian, independent of input n.
    H=np.diag(np.r_[np.full(2,OMEGA),np.repeat(k,2),np.repeat(k,2)]).astype(complex)
    for j in range(n):
        H[:2,2+2*j:4+2*j]=G*np.sqrt(2*w[j])*I2
        H[2+2*j:4+2*j,:2]=G*np.sqrt(2*w[j])*I2
    return {'k':k,'w':w,'raw_norm':raw_norm,'star':star,'initial':initial,
            'out':out,'dark':dark,'H':H,'U':unitary(H,TIME)}


def embedding(C,scalar_norm,dyn,at_time=False):
    D=psqrt(scalar_norm*I2-C.conj().T@C)
    if at_time:
        a=dyn['out'][0];bright=dyn['out'][1:];dark=dyn['dark']
    else:
        a=0.;bright=dark=np.sqrt(dyn['w'])
    return np.vstack([a*C]+[z*C for z in bright]+[z*D for z in dark])


def compare_values(actual,saved,path=''):
    if isinstance(actual,dict):
        assert actual.keys()==saved.keys(),path
        for key in actual:compare_values(actual[key],saved[key],path+'/'+key)
    elif isinstance(actual,list):
        assert len(actual)==len(saved),path
        for i,(x,y) in enumerate(zip(actual,saved)):compare_values(x,y,path+'/'+str(i))
    elif isinstance(actual,float):
        assert np.isclose(actual,saved,rtol=2e-10,atol=2e-12),(path,actual,saved)
    else:assert actual==saved,(path,actual,saved)


def calculate():
    # Exact analytic constants: no floating quadrature supplies the theorem.
    alpha2=2*A1*A1/A2
    assert alpha2>F(12,25)**2
    dyson_upper=F(1,7)**3/(6*(1-F(1,7)**2/20))
    assert dyson_upper==F(980,2014782) and dyson_upper<F(1,2000)
    amplitude_lower=F(1,10)*F(12,25)*F(7,8)-F(1,2000)
    gap_lower=amplitude_lower**2
    assert amplitude_lower==F(83,2000) and gap_lower==F(6889,4000000)
    directions=[np.array(v,float) for v in
                [(0,0,1),(0,0,-1),(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(1,2,3),(-2,1,1)]]
    dyn=dynamics(24);dense=dynamics(48)
    Q=np.zeros_like(dyn['H']);Q[:2,:2]=I2
    H0=np.diag(np.diag(dyn['H']))
    commutator=np.linalg.norm(dyn['H']@Q-Q@dyn['H'],2)
    assert abs(commutator-G*np.sqrt(2))<1e-12
    dephased_H=Q@dyn['H']@Q+(np.eye(len(Q))-Q)@dyn['H']@(np.eye(len(Q))-Q)
    assert np.linalg.norm(dephased_H-H0,2)<1e-12
    eigen=np.linalg.eigvalsh(dyn['H'])
    assert eigen[0]>6/7 and eigen[-1]<15/7
    common_p=abs(dyn['out'][0])**2*float(A1*A1/A2)
    dense_p=abs(dense['out'][0])**2*float(A1*A1/A2)
    assert common_p>float(gap_lower)
    # Independent RK4 check of the scalar star dynamics (fixed, modest scale).
    y=dyn['initial'].astype(complex);dt=TIME/256
    def rhs(z):return -1j*(dyn['star']@z)
    for _ in range(256):
        q1=rhs(y);q2=rhs(y+dt*q1/2);q3=rhs(y+dt*q2/2);q4=rhs(y+dt*q3)
        y+=dt*(q1+2*q2+2*q3+q4)/6
    rk_error=float(np.linalg.norm(y-dyn['out']))
    assert rk_error<1e-10
    checks=[];unpolarized=[];max_angle=0.;max_intertwine=0.;max_reference=0.
    rng=np.random.default_rng(1046)
    rho_seed=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
    rhoQR=rho_seed@rho_seed.conj().T;rhoQR/=np.trace(rhoQR)
    rhoR=np.trace(rhoQR.reshape(2,2,2,2),axis1=0,axis2=2)
    for direction in directions:
        a=angular(direction);n=a['n'];C=a['C']
        expected_C=float(A1)/np.sqrt(float(2*A2))*pauli(a['eps'])
        angle_err=max(abs(a['norm']-1),np.linalg.norm(a['AAstar']-2*I2,2),
                      np.linalg.norm(C-expected_C,2),np.linalg.norm(np.cross(a['mean_direction'],n)))
        max_angle=max(max_angle,float(angle_err))
        assert angle_err<2e-12 and np.dot(a['mean_direction'],n)>.5
        J=embedding(C,1.,dyn)
        exact=dyn['U']@J
        reduced=embedding(C,1.,dyn,True)
        identity_error=float(np.linalg.norm(exact-reduced,2))
        max_intertwine=max(max_intertwine,identity_error)
        E=exact[:2].conj().T@exact[:2]
        expected=common_p*(I2-pauli(n))/2
        assert np.linalg.norm(E-expected,2)<2e-12
        assert np.linalg.norm(exact.conj().T@exact-I2,2)<2e-12
        # Full instrument output on a genuinely unknown qubit/reference input.
        We=np.zeros_like(exact);We[:2]=exact[:2]
        Wg=exact-We
        output=[]
        for K in [We,Wg]:
            KR=np.kron(K,I2);output.append(KR@rhoQR@KR.conj().T)
        ref=np.trace(sum(output).reshape(len(Q),2,len(Q),2),axis1=0,axis2=2)
        referr=float(np.linalg.norm(ref-rhoR,2));max_reference=max(max_reference,referr)
        p_actual=float(np.trace(output[0]).real)
        p_effect=float(np.trace(np.kron(E,I2)@rhoQR).real)
        assert abs(p_actual-p_effect)<2e-12 and referr<2e-12
        bright=psqrt((I2-pauli(n))/2)
        # The absorbed spin is opposite to the original bright spin.
        absorbed=exact[:2]@bright
        assert np.linalg.norm(pauli(n)@absorbed-absorbed,2)<2e-12
        checks.append({'direction':n.tolist(),'effect_eigenvalues':np.linalg.eigvalsh(E).tolist(),
                       'probability_unknown_reference':p_actual,'reference_residual':referr,
                       'full_unitary_vs_exact_angular_compression':identity_error})
        # Finite-rank, genuinely pointwise unpolarized normal photon source.
        Eu=np.zeros((2,2),complex);total=np.zeros((2,2),complex)
        for j in range(3):
            Cj=pauli(a['B'][:,j])/np.sqrt(2)
            Jj=embedding(Cj,float(a['S2'][j,j]),dyn)/np.sqrt(float(2*J2))
            Wj=dyn['U']@Jj
            Eu+=Wj[:2].conj().T@Wj[:2]
            total+=Wj.conj().T@Wj
        expected_q=abs(dyn['out'][0])**2*float(2*A1*A1+B1*B1)/float(4*J2)
        assert np.linalg.norm(Eu-expected_q*I2,2)<2e-12
        assert np.linalg.norm(total-I2,2)<2e-12
        unpolarized.append({'direction':n.tolist(),'effect_eigenvalues':np.linalg.eigvalsh(Eu).tolist(),
                            'anisotropic_residual':float(np.linalg.norm(Eu-np.trace(Eu)*I2/2,2))})
    # Terminal readout's nonzero possible energy cost, not donated for free.
    eta=np.zeros(len(Q),complex);eta[0]=1/np.sqrt(2)
    for j,w in enumerate(dyn['w']):eta[2+2*j]=np.sqrt(w/2)
    energy_delta=float(np.vdot(eta,(dephased_H-dyn['H'])@eta).real)
    assert abs(energy_delta+G*np.sqrt(2))<1e-12
    return {
        'round':1046,'date':'2026-10-08','new_science_groups':1,'new_cognitive_axioms':0,
        'roadmap_completed_here':False,'all_scientific_checks_passed':True,
        'parameters':{'dimension_adopted':3,'frequency_window':[1.,2.],'Omega':OMEGA,'g':G,'time':TIME,
                      'angle_support_cosine_lower':.5,'quadrature_radial_sizes':[24,48],
                      'angular_gauss_per_interval':12,'azimuth_nodes':16},
        'exact_certificates':{'a1':str(A1),'a2':str(A2),'b1':str(B1),'b2':str(B2),'j2':str(J2),
                              'alpha_squared':str(alpha2),'dyson_remainder_upper':str(dyson_upper),
                              'amplitude_lower':str(amplitude_lower),'antipodal_gap_lower':str(gap_lower),
                              'H_lower':'6/7','H_upper':'15/7'},
        'continuous_gap_rational_lower':float(gap_lower),
        'calibration_p_t_one':float(common_p),'radial_refined_p':float(dense_p),
        'radial_refinement_difference':float(abs(common_p-dense_p)),
        'raw_radial_normalization_errors':[abs(dyn['raw_norm']-1),abs(dense['raw_norm']-1)],
        'radial_mean_frequency':float(dyn['w']@dyn['k']),
        'radial_frequency_second_moment':float(dyn['w']@(dyn['k']**2)),
        'independent_RK4_state_error':rk_error,
        'max_physical_angular_integral_residual':max_angle,
        'max_full_unitary_vs_compressed_map_residual':max_intertwine,
        'max_reference_reduction_residual':max_reference,
        'finite_H_spectrum_bounds':[float(eigen[0]),float(eigen[-1])],
        'terminal_measurement_commutator_norm':float(commutator),
        'terminal_measurement_energy_cost_bound':float(G*np.sqrt(2)),
        'terminal_cost_saturating_control':energy_delta,
        'direction_instruments':checks,'unpolarized_counterexamples':unpolarized,
        'historical_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in HISTORY},
        'scope':[
            '3D Maxwell kinematics, doublets, band-limited isotropic RWA vertex, pinning and preparation adopted',
            'rank-one and all-direction statements analytically exact for this parent, not grid-derived',
            'qubit is receiver spin; photon momentum directions are physical input labels, not its polarization qubit',
            'all outcomes retain joint atom-field state and passive reference',
            'terminal projection is adopted and does not commute with H; no apparatus energy closure',
            'only direction-instrument part of old spatial premises; no displacement/neighborhood or dimensional generation',
            'quadrature and RK refinements are numerical checks, not a QED matching error',
            'normal unpolarized source is rank three, not a nonnormal exact-momentum mixture',
        ]}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=calculate()
    if args.write:
        with TARGET.open('x',encoding='utf8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    else:compare_values(result,json.loads(TARGET.read_text(encoding='utf8')))
    print(json.dumps({'round':1046,'all_checks_passed':True,'mode':'exclusive_write' if args.write else 'read_only_compare',
                      'p_t_one':result['calibration_p_t_one'],'certified_gap_lower':result['continuous_gap_rational_lower'],
                      'max_angular_residual':result['max_physical_angular_integral_residual']},ensure_ascii=False))
