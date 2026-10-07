"""768: full BRST block and relative-source algebra calibrations.

Analytic Green kernels, Hadamard singularities and the complete variational
source are proved/mapped in note768; finite matrices do not simulate them.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_covariant_gauge_complex as old
import joint_scalar_propagation_matching as scalar

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_brst_relative_source_results.json'


def mx(x):
    return float(np.max(np.abs(x)))


def herm(x):
    return (x+x.conj().T)/2


def block_matrix(blocks, sizes):
    rows=[]
    for i,row in enumerate(blocks):
        rows.append([np.zeros((sizes[i],sizes[j]),complex) if b is None else b
                     for j,b in enumerate(row)])
    return np.block(rows)


def setup(k):
    phi=np.array([.13,.64,-.11,.08,.37])
    weights=np.repeat([.73,1.17,.41],[8,3,1])
    H,W,Kr,Kar,Pr,k2=old.principal(k,phi,weights)
    # Fourier symbol and the W -> -W auxiliary convention.
    K=1j*Kr; Ks=-1j*Kar; P=-Pr.astype(complex)
    sizes=[63,16,16,16]; I=np.eye(16)
    D0=Ks@K
    L=block_matrix([[P,K,None,None],[Ks,-I,None,None],
                    [None,None,None,D0],[None,None,D0,None]],sizes)
    gamma=block_matrix([[None,None,K,None],[None,None,None,None],
                        [None,None,None,None],[None,I,None,None]],sizes)
    pairing=block_matrix([[H,None,None,None],[None,-W,None,None],
                          [None,None,-W,None],[None,None,None,-W]],sizes)
    ga=np.linalg.solve(pairing,gamma.conj().T@pairing)
    return H,W,K,Ks,P,D0,L,gamma,ga,pairing,k2


def green_blocks():
    H,W,K,Ks,P,D0,L,gamma,ga,pairing,k2=setup([1.3,.2,-.3,.4])
    G1=np.eye(63)/(-k2);G0=np.eye(16)/(-k2)
    G=block_matrix([[G1,K@G0,None,None],[Ks@G1,None,None,None],
                    [None,None,None,G0],[None,None,G0,None]],[63,16,16,16])
    errors=dict(auxiliary_adjoint=mx(Ks-np.linalg.solve(-W,K.conj().T@H)),
                nilpotency=mx(gamma@gamma), BRST_symmetry=mx(ga@L-L@gamma),
                formal_self_adjoint=mx(pairing@L-L.conj().T@pairing),
                inverse_left=mx(L@G-np.eye(111)),
                inverse_right=mx(G@L-np.eye(111)))
    assert max(errors.values())<1e-13
    return dict(errors=errors,
                scope='One noncharacteristic full 111-by-111 principal block. Inverse identities only; retarded support and variable-coefficient propagation are analytic.')


def covariance_lift():
    d=np.array([1.,2.,-3.]); d/=np.linalg.norm(d)
    H,W,K,Ks,P,D0,L,gamma,ga,pairing,k2=setup(np.r_[1.,d])
    G=np.linalg.inv(K.conj().T@K)
    M=np.linalg.inv(K.conj().T@H@H@K)
    R=herm(np.eye(63)-K@G@K.conj().T-H@K@M@K.conj().T@H)
    val,vec=np.linalg.eigh(R)
    E=vec[:,val>.5]
    noise=E@np.diag(np.linspace(.012,.043,E.shape[1]))@E.conj().T
    delta=np.linalg.solve(H,noise)
    n=.17
    def lift(a1,a0,ghost_sign):
        return block_matrix([[a1,K@a0,None,None],[a0@Ks,None,None,None],
                             [None,None,None,ghost_sign*a0],
                             [None,None,ghost_sign*a0,None]],[63,16,16,16])
    plus=lift((1+n)*np.eye(63)+delta,(1+n)*np.eye(16),1)
    minus=lift(n*np.eye(63)+delta,n*np.eye(16),-1)
    causal_weight=lift(np.eye(63),np.eye(16),1)
    parity=np.diag(np.r_[np.ones(79),-np.ones(32)])
    extension=block_matrix([[delta,None,None,None],[None,None,None,None],
                           [None,None,None,None],[None,None,None,None]],[63,16,16,16])
    errors=dict(exact_intertwining=mx(delta@K),
                subsidiary_bisolution=mx(Ks@delta),
                original_P_bisolution=mx(P@delta),
                two_point_equation_plus=mx(L@plus),
                two_point_equation_minus=mx(L@minus),
                graded_CCR=mx(plus-parity@minus-causal_weight),
                covariance_adjoint=mx(pairing@plus-plus.conj().T@pairing),
                smooth_extension_equation=mx(L@extension),
                positive_BRST_intertwining=mx(gamma@plus-plus@ga),
                negative_BRST_intertwining=mx(gamma@minus+minus@ga))
    physical_minima=[float(np.linalg.eigvalsh(herm(E.conj().T@H@a@E))[0])
                    for a in ((1+n)*np.eye(63)+delta,n*np.eye(63)+delta)]
    assert max(errors.values())<1e-13 and min(physical_minima)>0
    assert E.shape[1]==31
    return dict(errors=errors, physical_minima=physical_minima,
                auxiliary_state_difference_norm=mx(extension[63:,:]),
                scope='One null-frequency algebra calibration. Mode occupation is not a continuum Hadamard prescription. Physical covariance changes occupy only the original field block.')


def potential_gradient_jet(phi,v):
    # Original572 potential V/F^2; coefficients are exact Taylor jets through
    # t^2 of grad U(phi+t*v), not finite differences.
    mat,u,_=scalar.parameters()
    group=np.array([[1.,1.,1.,1.,0.],[0.,0.,0.,0.,1.]])
    delta=group@(phi*phi)-u
    d1=2*group@(phi*v);d2=group@(v*v)
    w0=group.T@(mat@delta);w1=group.T@(mat@d1);w2=group.T@(mat@d2)
    gv=[phi*w0,v*w0+phi*w1,v*w1+phi*w2]
    V=[float(delta@mat@delta/4),
       float(d1@mat@delta/2),
       float((d1@mat@d1+2*d2@mat@delta)/4)]
    F=2-phi@phi/6; a=-(phi@v)/(3*F);b=-(v@v)/(6*F)
    def inverse_power(p):
        return np.array([1.,-p*a,p*(p+1)*a*a/2-p*b])*F**(-p)
    f2,f3=inverse_power(2),inverse_power(3)
    second=[V[0]*phi,V[1]*phi+V[0]*v,V[2]*phi+V[1]*v]
    jet=[]
    for k in range(3):
        jet.append(sum(gv[j]*f2[k-j]+(2/3)*second[j]*f3[k-j]
                       for j in range(k+1)))
    return jet


def source_jet_identity():
    phi=np.array([.13,.64,-.11,.08,.37]);v=np.array([.27,-.14,.31,.18,-.09])
    t=old.scalar_generators()[8]
    R=t@phi
    e,hv,half_cvv=potential_gradient_jet(phi,v)
    _,hr,_=potential_gradient_jet(phi,R)
    cvv=2*half_cvv
    source_term=float(R@cvv)
    contact_term=float(2*hv@(t@v))
    errors=dict(potential_invariance=abs(float(e@R)),
                first_Noether=mx(hr+t.T@e),
                second_Noether=abs(source_term+contact_term))
    assert max(errors.values())<2e-14 and abs(source_term)>1e-5
    return dict(errors=errors, quadratic_source_term=source_term,
                necessary_contact_term=contact_term,
                scope='Original H5 potential jet only, at a declared off-shell local calibration point. It verifies the differentiated identity, not full quantum Ward or a conserved isolated potential source.')


def run():
    result=dict(round=768, tests_run=3, failures=0, errors=0,
                complete_BRST_blocks=green_blocks(),
                covariance_extension=covariance_lift(),
                original_potential_source_identity=source_jet_identity())
    deps=('research_note_734.md','research_note_735.md','research_note_754.md',
          'research_note_756.md','research_note_764.md','research_note_765.md',
          'research_note_766.md','research_note_767.md',
          'joint_covariant_gauge_complex.py','round768_drafts/joint_source_entry.md')
    result['dependency_hashes']={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in deps}
    result['scope']=('Same complete on-shell linear bosonic theory has a BRST two-point extension '
                     'with the original physical state. Common smooth physical covariance changes '
                     'have unchanged auxiliary blocks, give finite complete bosonic source differences '
                     'obeying the joint linear Noether identity, and enter the inherited formal '
                     'first-order constrained response. Absolute all-sector quantum Ward, interacting '
                     'QME, original Q continuum equivalence and autonomous records remain open.')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
    else:
        assert result==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:result[k] for k in ('round','tests_run','failures','errors')}))
