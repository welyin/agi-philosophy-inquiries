"""598: finite CAR completion of the original confining Gauss model.

Checks quotient representations/Yukawa intertwiners, the original potential
bound, and a finite matrix check of the record/source identity. No continuum
chiral regulator, full Fock spectrum, or Einstein constraint is computed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_curved_quantum_source as original
import joint_quotient_gauge_completion as gauge

HERE = Path(__file__).resolve().parent
TARGET = HERE / 'joint_fermion_gauss_completion_results.json'
SLICES = dict(Q=slice(0,12), u=slice(12,18), d=slice(18,24),
              L=slice(24,28), e=slice(28,30), nu=slice(30,32))
Y = dict(u=1.2+.3j, d=.5-.1j, e=.2+.07j, nu=.4-.12j, s=.31+.09j)


def representation(C, W, z):
    spin = np.eye(2)
    blocks = dict(Q=np.kron(np.kron(C,W),spin)*z,
                  u=np.kron(C,spin)*z**4, d=np.kron(C,spin)*z**-2,
                  L=np.kron(W,spin)*z**-3, e=spin*z**-6, nu=spin)
    out = np.zeros((32,32),complex)
    for name, block in blocks.items(): out[SLICES[name],SLICES[name]] = block
    return out


def mass_matrices(phi):
    """Number-conserving Dirac block h and antisymmetric Majorana block Delta.

    Physical right-handed convention; X/sqrt(F) is chosen Einstein-frame
    mass dependence, with the usual Higgs sqrt(2) absorbed in diagnostic Y.
    """
    X = phi[:2]+1j*phi[2:4]
    tilde = np.array([X[1].conjugate(),-X[0].conjugate()])
    h = np.zeros((32,32),complex)
    for left, right, col, colors in [('Q','u',tilde,3),('Q','d',X,3),
                                     ('L','nu',tilde,1),('L','e',X,1)]:
        block = Y[right]*np.kron(np.eye(colors),np.kron(col[:,None],np.eye(2)))
        h[SLICES[left],SLICES[right]] = block/np.sqrt(original.F(phi))
    h = h+h.conj().T
    delta = np.zeros_like(h)
    delta[SLICES['nu'],SLICES['nu']] = (
        Y['s']*phi[4]/np.sqrt(original.F(phi))*np.array([[0.,1.],[-1.,0.]]))
    return h,delta


def representation_check():
    rng = np.random.default_rng(598)
    zeta = representation(np.eye(3)*np.exp(2j*np.pi/3),-np.eye(2),np.exp(1j*np.pi/3))
    center_error = float(np.max(abs(zeta-np.eye(32))))
    errors = []
    for _ in range(24):
        phi = rng.normal(size=5)*.4
        C=gauge.group_exp(rng.normal(size=8),3)
        W=gauge.group_exp(rng.normal(size=3),2)
        z=np.exp(1j*rng.normal()); R=representation(C,W,z)
        X=z**3*W@(phi[:2]+1j*phi[2:4])
        transformed=np.r_[X.real,X.imag,phi[4]]
        h,d=mass_matrices(phi); ht,dt=mass_matrices(transformed)
        errors.extend([np.max(abs(ht-R@h@R.conj().T)),
                       np.max(abs(dt-R@d@R.T)),np.max(abs(R.conj().T@R-np.eye(32))),
                       np.max(abs(d+d.T)),np.max(abs(h-h.conj().T))])
    # Nontrivial charged hopping with independently transformed endpoints.
    Cs=[gauge.group_exp(rng.normal(size=8),3) for _ in range(3)]
    Ws=[gauge.group_exp(rng.normal(size=3),2) for _ in range(3)]
    zs=np.exp(1j*rng.normal(size=3))
    transformed_link=representation(Cs[0]@Cs[2]@Cs[1].conj().T,
                                   Ws[0]@Ws[2]@Ws[1].conj().T,zs[0]*zs[2]/zs[1])
    hopping_error=float(np.max(abs(transformed_link-
        representation(Cs[0],Ws[0],zs[0])@representation(Cs[2],Ws[2],zs[2])@
        representation(Cs[1],Ws[1],zs[1]).conj().T)))
    assert max(center_error,max(errors),hopping_error)<2e-13
    return dict(one_particle_modes_per_generation_per_node=32,
                physical_Weyl_multiplicity_per_generation=16,center_error=center_error,
                maximal_covariance_error=float(max(errors)),hopping_covariance_error=hopping_error,
                sample_count=24,diagnostic_Y={k:[v.real,v.imag] for k,v in Y.items()},
                no_continuum_chirality_or_anomaly_claim=True)


def confinement_check():
    L,u,_=original.lattice.scalar.parameters()
    lam=float(np.linalg.det(L)/np.trace(L)); D=float(6*original.M-u.sum())
    A=32/(lam*D*D); B0=144/(D*D)
    rng=np.random.default_rng(5981)
    directions=rng.normal(size=(600,5));directions/=np.linalg.norm(directions,axis=1)[:,None]
    f=np.geomspace(2e-6,original.M*.999,600)
    phi=directions*np.sqrt(6*(original.M-f))[:,None]
    actual_f=original.F(phi);U=original.node_potential(phi)
    ratios=actual_f**-2/(A*U+B0)
    assert np.max(ratios)<=1+1e-11
    # Direct finite one-particle mass norm versus the inherited confinement.
    rows=[]
    direction=np.array([1.,2.,-.4,.6,.8]);direction/=np.linalg.norm(direction)
    for f0 in (1e-1,1e-2,1e-3,1e-4,1e-5):
        p=direction*np.sqrt(6*(original.M-f0));h,d=mass_matrices(p)
        norm=float(np.linalg.norm(h,2)+np.linalg.norm(d,2))
        potential=float(original.node_potential(p))
        rows.append(dict(F=float(original.F(p)),one_particle_coefficient_norm=norm,
                         original_potential=potential,coefficient_over_potential=norm/potential))
    assert all(rows[i+1]['coefficient_over_potential']<rows[i]['coefficient_over_potential']
               for i in range(len(rows)-1))
    young=[]
    for c,eps in ((.1,.3),(2.,.01),(11.,.2)):
        constant=3*c**(4/3)/(4**(4/3)*eps**(1/3))
        stationary=(c/(4*eps))**(4/3)
        gap=eps*stationary+constant-c*stationary**.25
        grid=stationary*np.geomspace(1e-8,1e8,401)
        assert np.min(eps*grid+constant-c*grid**.25)>-1e-10
        assert abs(gap)<1e-10
        young.append(dict(c=c,epsilon=eps,constant=constant,stationary_equality_error=float(abs(gap))))
    return dict(original_M=original.M,lambda_lower_bound=lam,D=D,A=A,B0=B0,
                sampled_max_bound_ratio=float(np.max(ratios)),boundary_rows=rows,
                young_checks=young,full_form_bound_is_analytic=True,
                no_cutoff_uniform_constant_claim=True)


def instrument_check():
    # Finite matrix fixture checks algebra, not the continuum spectrum or D_B.
    s=np.linspace(-.8,.9,6);m=len(s); n=4
    T=np.diag(np.full(m,2.))-np.diag(np.ones(m-1),1)-np.diag(np.ones(m-1),-1)
    Hb=np.kron(T+np.diag(.2+s*s),np.eye(n))
    B=np.zeros_like(Hb,dtype=complex);dB=np.zeros_like(B)
    for i,ss in enumerate(s):
        p=np.array([.4,.3,-.2,.1,ss]);h,d=mass_matrices(p)
        # Different charged components mix; matrix potential does not commute
        # with a generic internal observable but does commute with scalar L.
        ids=[24,25,30,31];block=h[np.ix_(ids,ids)]
        B[i*n:(i+1)*n,i*n:(i+1)*n]=block
        dB[i*n:(i+1)*n,i*n:(i+1)*n]=(1+ss)*block
    Ls=[np.kron(np.diag(np.sqrt(.5+sign*np.sin(s)/4)),np.eye(n)) for sign in (-1,1)]
    def channel(A):return sum(L@A@L for L in Ls)
    Db=channel(Hb)-Hb
    H=Hb+B;G=-6*Hb+dB
    injection_error=float(np.max(abs(channel(H)-H-Db)))
    source_error=float(np.max(abs(channel(G)-G+6*Db)))
    complete_error=float(np.max(abs(sum(L@L for L in Ls)-np.eye(m*n))))
    assert max(injection_error,source_error,complete_error)<2e-14
    assert np.linalg.norm(Db)>1e-3
    return dict(matrix_dimension=m*n,injection_identity_error=injection_error,
                source_identity_error=source_error,completeness_error=complete_error,
                fixture_only=True,continuum_identity_proved_in_note=True)


def run():
    rep=representation_check();conf=confinement_check();inst=instrument_check()
    deps=('research_note_531.md','research_note_553.md','research_note_569.md','research_note_574.md',
          'research_note_579.md','research_note_591.md','research_note_593.md',
          'joint_curved_quantum_source.py','joint_quotient_gauge_completion.py',
          'joint_condition_compression_update_597.md')
    return dict(round=598,tests_run=3,failures=0,errors=0,representation=rep,confinement=conf,
                instrument=inst,dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
                scope=dict(finite_graph_CAR_extension=True,original_bosonic_potential_retained=True,
                           common_semibounded_compact_resolvent_H_proved=True,
                           record_and_average_source_identity_retained=True,
                           added_fermion_statistics_and_dynamics_declared=True,
                           no_chiral_continuum_or_GR_or_unique_ground_claim=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as stream:
            stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(dict(round=598,tests=3,all_passed=True,
                         covariance=result['representation']['maximal_covariance_error'],
                         source_error=result['instrument']['source_identity_error'])))
