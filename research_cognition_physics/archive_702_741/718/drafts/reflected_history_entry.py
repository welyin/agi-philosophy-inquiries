"""718 entry: original mass reflection and retention of the actual outcome."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_smooth_mode_contract as shared
car=shared.car;matter=car.matter
TARGET=HERE/'reflected_history_entry_results.json'


def exponential(H,t):
    d,v=np.linalg.eigh(H)
    return (v*np.exp(-1j*t*d))@v.conj().T


def quadratic_matrix(h,d):
    n=len(h);B=np.zeros((2**n,2**n),complex)
    for col in range(2**n):
        for row,z in car.quadratic({col:1.},h,d).items():B[row,col]=z
    return B


def run():
    # An exact invariant neutral sector of the original one-node mass matrix.
    # This conditional fixture is not the full dynamic bosonic Hamiltonian.
    phi=np.array([.4,0.,0.,0.,.35])
    h,d=matter.mass_matrices(phi)
    ids=np.array([26,27,30,31]);other=np.array([i for i in range(32) if i not in ids])
    assert np.linalg.norm(h[np.ix_(other,ids)])<1e-14
    assert np.linalg.norm(d[np.ix_(other,ids)])<1e-14
    H=quadratic_matrix(h[np.ix_(ids,ids)],d[np.ix_(ids,ids)])
    n=np.diag([float(bool(i&(1<<2))) for i in range(16)])
    m=np.diag([float(bool(i&(1<<3))) for i in range(16)])
    R=np.eye(16)-2*n
    Hr=R@H@R
    psi=np.zeros(16,complex);psi[0]=np.sqrt(.7);psi[12]=np.sqrt(.3)*np.exp(.7j)
    rho=np.outer(psi,psi.conj())
    rows=[]
    for t in (.3,1.2,2.7):
        V=m@exponential(H,t);Vr=m@exponential(Hr,t)
        ordinary=float(np.trace(V@rho@V.conj().T).real)
        reflected=float(np.trace(Vr@rho@Vr.conj().T).real)
        cross=np.trace(V@rho@Vr.conj().T@R)
        actual=[];matched=[]
        for sign in (1,-1):
            K=(np.eye(16)+sign*R)/2
            actual.append(float(np.trace(V@K@rho@K@V.conj().T).real))
            matched.append((ordinary+reflected+2*sign*cross.real)/4)
        error=float(max(abs(np.array(actual)-matched)))
        loss=float(max(abs(np.array(actual)-(ordinary+reflected)/4)))
        assert error<1e-13 and loss>1e-4
        rows.append(dict(time=t,actual_joint_probabilities=actual,
                         two_branch_with_interference=matched,
                         original_history=ordinary,reflected_history=reflected,
                         interference_real=float(cross.real),identity_error=error,
                         lost_actual_record_if_cross_omitted=loss))
    rng=np.random.default_rng(7182)
    points=np.array([[.4,-.3,.2,.1,.35],[-.2,.15,.1,-.25,.3],[.22,.34,-.12,.08,-.27]])
    links=[shared.group.sample(rng) for _ in range(3)]
    mass,dd,hop=shared.coefficients(points,links,np.array([.21,.27,.18]))
    hh=mass+hop
    f,_=shared.profile(np.array([.8,1.1,.6]),rng.normal(size=(3,2))+1j*rng.normal(size=(3,2)))
    U=np.eye(96)-2*np.outer(f,f.conj())
    hprime,dprime=U@hh@U,U@dd@U.T
    dh,delta=shared.difference_coeff(hh,dd,f)
    coeff=max(np.linalg.norm(hprime-hh-2*dh),np.linalg.norm(dprime-dd-2*delta))
    nambu=np.block([[hh,dd],[dd.conj().T,-hh.T]])
    transformed=np.block([[hprime,dprime],[dprime.conj().T,-hprime.T]])
    lifted=np.block([[U,np.zeros_like(U)],[np.zeros_like(U),U.conj()]])
    finite_time=float(np.linalg.norm(exponential(transformed,.83)-
                                    lifted@exponential(nambu,.83)@lifted))
    assert max(coeff,finite_time)<1e-12
    names=('research_note_623.md','research_note_624.md','research_note_634.md',
           'research_note_704.md','research_note_717.md','joint_smooth_mode_contract.py')
    return dict(entry_round=718,new_formal_round=False,
                actual_initial_projective_outcome_retained=True,
                conditional_neutral_sector_modes=4,conditional_Fock_dimension=16,
                conditional_history_rows=rows,full_coefficient_modes=96,
                reflection_difference_coefficient_error=float(coeff),
                conditional_Nambu_finite_time_transport_error=finite_time,
                full_boson_history_identity_proved_analytically=True,
                frozen_coefficient_diagnostics_not_full_physical_propagation=True,
                no_new_Hamiltonian_or_autonomous_device_claim=True,
                dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in names})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args()
    r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
