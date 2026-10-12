"""Independent qualification diagnostics; default is read-only.

Only --write-results writes this file's matching JSON.  The finite CAR test
checks a restriction of the normal state, not a replacement field model.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
FROZEN = {
    'receiver1137_preflight_author.md': '24c1903ddaedd651ca2ad137158c7c5511b92069c20a5e2bb448ed14b30b28a7',
    'receiver1137_preflight_author.py': 'b5c750404c2cfe793356811b244d9403859f7501ad3bb8a49bce895b5effa90d',
    'receiver1137_preflight_author.json': '71606653ed4dd237a76c538725d9a525ccb6880bc759e82edfefecc7b4d255d5',
    '../1133/proof.md': '476dc1d1258b819a72e2650684e1096cfa800d7c7c885c77cc1d41a07ebbde8d',
    '../1134/proof.md': '4b75693a9efeb7c2747152763da6fb3bed2f1e9b25418e9827cdd615f49c5666',
    '../1135/proof.md': 'f3577777e8983a02d02bf66ef5f869b207ac4f94df1610202c1ce6e995e17068',
    '../1136/proof.md': '360c450f85f40798c6fb9524b00124364b604434d814c770503e0e97e2196462',
}


def pack(values):
    return {'fractions': {k: str(v) for k, v in values.items()},
            'decimal': {k: float(v) for k, v in values.items()}}


def exact_budget():
    M = F(10**80)
    sig, a, T, A = F(1, 10**30), F(1, 10**9), F(4), F(10**34)
    lam, kap, eps, dc = F(11, 7), F(11, 14), F(1, 10**5), F(1, 10**6)
    C = F(33536000000000000000)
    v = 16*lam*C + F(558, 7)
    mass = 3*M-v
    CY = F(24780800000000000000000, 49)+16316
    p1, p2 = 268+94*A, 28000+188*268*A+12000*A*A
    A1, A2, B2 = 1/sig+p1, 1/sig**2+2*p1/sig+p2, p1/sig+p2
    ds = 3*T/(eps*A)*(lam*C*512+42*8/dc+(lam+kap)*512)
    dy = CY*(2*T*sig+T*T/(2*M*sig))
    dd = 28*T/(2*M*M*sig*sig)
    bt = T*A2/(2*mass)+T*(A2+2*B2)/(2*M)+CY*A1*T*T/(2*M)
    l1, l2, C1, C2 = 40/a, 2000/a**2, F(36000), F(10**7)
    drot = (T*l2+2*C1*T*T*l1+F(3, 2)*C2*T*T+C1*C1*T**3)/(2*M)
    orient, dm = T*C1*a, 40*lam*eps
    rest = 42*dc+ds+dy+dd+bt+drot+orient
    task = 2*dm+rest
    author = json.loads((HERE/'receiver1137_preflight_author.json').read_text(encoding='utf-8'))
    assert task == F(author['task_rows']['code']['fractions']['complete_task_error'])
    assert drot == F(author['angular']['fractions']['delta_rotor'])
    # Separate conservative probability bounds; no normalized branch division.
    gap = F(1, 2)-4*task
    N = F(5*10**10)
    y = 2*sig+T*A1/mass
    delta = 32*lam
    ep = delta*(y+(N+F(1, 2))/mass)
    ry = 2*y+(N+F(1, 2))/mass
    ell = F(1, 1000)
    joint = task+F(3, 2)*ep+F(3, 2)*ry/ell
    assert mass > 2*M and gap > F(49, 100)
    assert ep < F(11, 10**29) and ry < F(5, 10**30)
    assert joint < F(145, 10**5)
    assert 4*F(6144, 18) < 40**2
    assert 104**2*F(6144, 18) < 2000**2
    assert 464+3300+2*sig*(600+94*A)+l1 < N
    return pack(dict(m_min=mass,delta_rotor=drot,delta_orientation=orient,
        extra_mode_single_particle_error=dm,common_remaining_error=rest,
        full_task_error=task,direction_probability_gap_lower=gap,
        mass_conjugacy_delta=delta,canonical_spin_bound=N,rest_position=y,
        bare_pointer_branch_commutator=ep,boost_position_RMS=ry,
        position_tail=ry**2/ell**2,
        complete_boost_joint_instrument_error=joint))


def spinor_axis():
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex)
    sz=np.diag([1,-1]).astype(complex)
    pauli=[sx,sy,sz]
    beta=np.kron(sz,np.eye(2))
    alpha=[np.kron(sx,s) for s in pauli]
    Sigma=[np.kron(np.eye(2),s) for s in pauli]
    vals,vecs=np.linalg.eigh(alpha[2])
    plus=vecs[:,vals>0]
    rng=np.random.default_rng(1137)
    samples=[np.array([1,0],complex),np.array([0,1],complex)]
    samples += [rng.normal(size=2)+1j*rng.normal(size=2) for _ in range(12)]
    worst=0.0
    for c in samples:
        xi=plus@(c/np.linalg.norm(c))
        b=np.array([np.vdot(xi,beta@Sigma[j]@xi).real for j in [0,1]])
        n=np.array([-b[1],b[0]])/np.linalg.norm(b) if np.linalg.norm(b)>1e-14 else np.array([1.,0.])
        xir=-1j*(n[0]*Sigma[0]+n[1]*Sigma[1])@xi
        errors=[abs(np.vdot(xi,xir)),abs(np.vdot(xi,beta@xir))]
        for k in [np.array([.3,.7,1.1]),np.array([2.,-.4,.8])]:
            kin=sum(k[j]*alpha[j] for j in range(3))
            energy=math.sqrt(1+np.dot(k,k))
            sp=(np.eye(4)-(kin+beta)/energy)/2
            sm=(np.eye(4)-(-kin+beta)/energy)/2
            errors.append(abs(np.vdot(xi,(sp+sm)@xir)))
        worst=max(worst,*errors)
    # Fixed x is NOT generally valid: beta Sigma_x restricted to alpha_z+.
    w,v=np.linalg.eigh(plus.conj().T@beta@Sigma[0]@plus)
    xi=plus@v[:,np.argmax(w)]
    bad=abs(np.vdot(xi,beta@(-1j*Sigma[0]@xi)))
    assert worst<1e-13 and bad>.99
    return dict(samples=len(samples),maximum_axis_residual=worst,
                naive_fixed_x_beta_cross=bad)


def car_ann(j,n=8):
    out=np.zeros((2**n,2**n),complex)
    for b in range(2**n):
        if (b>>j)&1:
            out[b^(1<<j),b]=(-1)**((b&((1<<j)-1)).bit_count())
    return out


def unitary(H,t):
    val,vec=np.linalg.eigh(H)
    return (vec*np.exp(-1j*t*val))@vec.conj().T


def full_car():
    c=[car_ann(j) for j in range(8)]
    I=np.eye(256,dtype=complex)
    # The other six physical CAR modes have their actual normal-state
    # maximally mixed marginal, with an explicit 64-dimensional purifier.
    J=np.zeros((256,64,2),complex)
    for e in range(64):
        J[(e<<2)|1,e,0]=1/8
        J[(e<<2)|2,e,1]=1/8
    G=np.zeros_like(I)
    for aa,bb in [(0,2),(1,3),(4,6),(5,7)]:
        G+=c[aa].conj().T@c[bb]+c[bb].conj().T@c[aa]
    evolved=(unitary(G,math.pi/2)@J.reshape(256,-1)).reshape(256,64,2)
    xp=np.array([[0,1],[1,0]],complex)
    results=[]
    residual=0.0
    for orientation,(i,j) in [('correct',(2,3)),('wrong',(6,7))]:
        nz=c[i].conj().T@c[i]
        cp=(c[i]+c[j])/math.sqrt(2)
        nx=cp.conj().T@cp
        yy=1j*(c[j].conj().T@c[i]-c[i].conj().T@c[j])
        cont=unitary(yy,math.pi/4)
        for label,N in [('Z',nz),('X',nx)]:
            incoming=np.zeros((256,2,64,2),complex)
            incoming[:,0]=evolved
            full=(unitary(np.kron(N,xp),math.pi/2)@incoming.reshape(512,-1)).reshape(256,2,64,2)
            full=np.einsum('ab,bpeq->apeq',cont,full)
            ideal=np.stack([cont@(I-N)@evolved.reshape(256,-1),
                            -1j*cont@N@evolved.reshape(256,-1)],axis=1).reshape(256,2,64,2)
            residual=max(residual,float(np.linalg.norm((full-ideal).reshape(-1,2),2)))
            p1=full[:,1].reshape(-1,2).conj().T@full[:,1].reshape(-1,2)
            target=(np.diag([1.,0.]) if label=='Z' else np.ones((2,2))/2) if orientation=='correct' else np.eye(2)/2
            assert np.linalg.norm(p1-target)<1e-12
            gram=full.reshape(-1,2).conj().T@full.reshape(-1,2)
            assert np.linalg.norm(gram-np.eye(2))<1e-12
            results.append(dict(orientation=orientation,setting=label,
                  effect_real=p1.real.tolist(),effect_imaginary_norm=float(np.linalg.norm(p1.imag)),
                  complete_output_gram_residual=float(np.linalg.norm(gram-np.eye(2)))))
    assert residual<1e-12
    z=np.diag([1.,0.]);x=np.ones((2,2))/2
    return dict(physical_CAR_modes=8,physical_dimension=256,purifier_dimension=64,
        complete_record_and_continuation_residual=residual,
        logical_effect_commutator_norm=float(np.linalg.norm(z@x-x@z,2)),rows=results,
        finite_restriction_not_replacement_field=True)


def calculate():
    sha={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in FROZEN}
    assert sha==FROZEN
    return dict(status='PASS',scope='Three qualification gates only; no numbered round.',
                frozen_sources=sha,exact_budget=exact_budget(),
                selected_axis_diagnostic=spinor_axis(),complete_CAR_diagnostic=full_car())


def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-7,abs_tol=1e-12),(a,b)
    else:assert a==b,(a,b)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args(); data=calculate()
    output=HERE/'receiver1137_preflight_review.json'
    if args.write_results:output.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:compare(data,json.loads(output.read_text(encoding='utf-8')))
    print('PASS: seven frozen identities; exact normal-rotor/rest/boost budgets; selected spin axis; complete eight-mode CAR records and continuation. Qualification only.')
