"""589 entry: positive full-metric coupling of the ORIGINAL quotient face potential.
Not a completed round; electric/matter operator domains and the full graph are pending.
"""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_quotient_gauge_completion as old
P=old.parameters()
TARGET=HERE/'full_metric_magnetic_entry_results.json'
PAIRS=((0,1),(0,2),(1,2))

def adj(U):
    T=old.generators(len(U))
    return 2*np.einsum('aij,jk,bkl,li->ab',T,U,T,U.conj().T).real

def channels(g):
    C,W,z=g;delta,coeff=old.coefficients(P)
    return [(delta*P['wq'],np.kron(C,W)*z),
            (delta*P['wq'],C.conj()*z**-4),(delta*P['wq'],C.conj()*z**2),
            (delta*P['wl'],W*z**-3),(delta*P['wl'],np.array([[z**6]])),
            (coeff[0],adj(C)),(coeff[1],adj(W)),(coeff[2],np.array([[z**6]]))]
    # The trivial neutrino representation contributes exactly zero.

def gram(gamma):
    gi=np.linalg.inv(gamma)
    return np.array([[gi[a,c]*gi[b,d]-gi[a,d]*gi[b,c] for c,d in PAIRS] for a,b in PAIRS])

def energies(faces,Gamma):
    byface=[channels(g) for g in faces]
    odd=np.zeros((3,3));even=np.zeros((3,3));single=np.zeros(3)
    for r in range(len(byface[0])):
        alpha=byface[0][r][0]
        O=[];E=[]
        for g in byface:
            R=g[r][1]
            O.append((R-R.conj().T)/(2j))
            E.append(np.eye(len(R))-(R+R.conj().T)/2)
        for f in range(3):
            for h in range(3):
                odd[f,h]+=.5*alpha*np.vdot(O[f],O[h]).real
                even[f,h]+=.5*alpha*np.vdot(E[f],E[h]).real
    energy=float(np.sum(Gamma*odd)+np.diag(Gamma)@np.diag(even))
    naive=float(np.sum(Gamma*(odd+even)))
    return energy,naive,np.diag(odd+even),float(np.trace(even))

def run():
    rng=np.random.default_rng(589)
    faces=[(old.group_exp(.5*rng.normal(size=8),3),old.group_exp(.4*rng.normal(size=3),2),
            np.exp(.17j*rng.normal())) for _ in range(3)]
    shape=np.array([[1.2,.22,-.1],[.22,.9,.16],[-.1,.16,1.1]])
    shape/=np.linalg.det(shape)**(1/3)
    psiface=np.array([1.04,1.11,.98])
    Gamma=gram(shape)/(psiface[:,None]*psiface[None,:])
    energy,naive,V,even=energies(faces,Gamma)
    oldV=np.array([old.potential(*g,P) for g in faces])
    scalar_error=float(np.max(abs(V-oldV)))
    conformal,_,_,_=energies(faces,np.diag(psiface**-2))
    conformal_error=abs(conformal-float(psiface**-2@oldV))
    eig=np.linalg.eigvalsh(Gamma);lower=eig[0]*sum(V);upper=eig[-1]*sum(V)
    assert scalar_error<1e-12 and conformal_error<1e-12 and lower<=energy<=upper
    # Independent quotient representative choices on each face.
    lifted=[]
    for j,(C,W,z) in enumerate(faces):
        n=j+1;lifted.append((np.exp(2j*np.pi*n/3)*C,(-1)**n*W,np.exp(1j*np.pi*n/3)*z))
    lift_error=abs(energies(lifted,Gamma)[0]-energy)
    # All faces have one common base node; a node gauge change is common conjugation.
    KC,KW=old.group_exp(rng.normal(size=8),3),old.group_exp(rng.normal(size=3),2)
    rotated=[(KC@C@KC.conj().T,KW@W@KW.conj().T,z) for C,W,z in faces]
    gauge_error=abs(energies(rotated,Gamma)[0]-energy)
    # Reverse one face orientation and the corresponding two-form coordinate.
    inverse=list(faces);C,W,z=faces[0];inverse[0]=(C.conj().T,W.conj().T,z.conjugate())
    S=np.diag([-1,1,1])
    changed,naive_changed,_,_=energies(inverse,S@Gamma@S)
    orientation_error=abs(changed-energy)
    naive_defect=abs(naive_changed-naive)
    assert max(lift_error,gauge_error,orientation_error)<1e-11 and naive_defect>1e-3
    # Small holonomy curvature with ALL original Lie-algebra factors.
    colour=rng.normal(size=(3,8));weak=rng.normal(size=(3,3));circle=rng.normal(size=3)
    expected=.5*(P['K'][0]*np.einsum('fg,fa,ga->',Gamma,colour,colour)
                  +P['K'][1]*np.einsum('fg,fa,ga->',Gamma,weak,weak)
                  +P['K'][2]*(circle@Gamma@circle))
    rows=[]
    for eta in (.04,.02,.01,.005):
        small=[(old.group_exp(eta*colour[j],3),old.group_exp(eta*weak[j],2),np.exp(1j*eta*circle[j])) for j in range(3)]
        value,_,_,ev=energies(small,Gamma)
        rows.append(dict(angle_scale=eta,energy_over_eta_squared=value/eta**2,
                         error=abs(value/eta**2-expected),even_channel_norm=ev))
    assert rows[-1]['error']<rows[0]['error']/40
    return dict(status='589 magnetic interface entry; not a completed round',checks_passed=True,
                exact_old_face_potential_error=scalar_error,exact_conformal_branch_error=conformal_error,
                eigenvalues_of_metric_Gram=eig.tolist(),positive_full_metric_energy=energy,
                comparison_bounds=[float(lower),float(upper)],independent_quotient_lift_error=lift_error,
                common_node_gauge_error=gauge_error,orientation_covariance_error=orientation_error,
                naive_full_embedding_orientation_defect=naive_defect,
                original_curvature_Hessian_prediction=float(expected),small_holonomy=rows,
                full_quantum_electric_scalar_and_domain_interface_pending=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
