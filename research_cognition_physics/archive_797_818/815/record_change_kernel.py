"""815 working: full ordered kernel of record-change square.
Original neutral mass coefficients; finite Fock/boson calibration only.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'803'))
import quasifree_pairing_probe as old
TARGET=HERE/'record_change_kernel_results.json'

def run():
    car=old.car(4);xi=car+[a.conj().T for a in car];eye=np.eye(16)
    def Q(k):return sum((xi[i].conj().T@xi[j]*k[i,j]/2 for i in range(8) for j in range(8)),np.zeros((16,16),complex))
    ids=[24,25,30,31,56,57,62,63]
    take=lambda m:m[np.ix_(ids,ids)]
    phi=np.array([0.,.6987151138647895,0.,0.,.545863690226873])
    phi2=phi*np.array([1.,.8,1.,1.,1.1])
    H0=Q(take(old.vertex.mass(phi)));H=Q(take(old.vertex.mass(phi2)))
    ev,v=np.linalg.eigh(H0);state=v[:,0]
    ev,v=np.linalg.eigh(H);U=(v*np.exp(-.63j*ev))@v.conj().T
    state=U@state
    expect=lambda a:np.vdot(state,a@state)
    P=np.array([[expect(a@b.conj().T) for b in xi] for a in xi])
    f=(np.eye(8)[:,2]+1j*np.eye(8)[:,7])/np.sqrt(2)
    cf=np.concatenate([f[4:].conj(),f[:4].conj()])
    D=-np.outer(f,f.conj())+np.outer(cf,cf.conj())
    c=sum((f[i].conjugate()*xi[i] for i in range(8)),np.zeros((16,16),complex))
    p=eye-c.conj().T@c
    keys=[take(old.vertex.dmass(phi2,np.eye(5)[i])) for i in (1,4)]
    L=[1j*(D@k-k@D) for k in keys];physical=[1j*(p@Q(k)-Q(k)@p) for k in keys]
    m=np.array([expect(a) for a in physical])
    connected=np.array([[np.trace((np.eye(8)-P)@a@P@b)/2 for b in L] for a in L])
    kernel=connected+np.outer(m,m)
    direct_kernel=np.array([[expect(a@b) for b in physical] for a in physical])
    X=np.array([[0,1],[1,0]],complex);Y=np.array([[0,-1j],[1j,0]],complex)
    bosons=[X,Y];vac=np.array([1.,0.]);WB=np.array([[np.vdot(vac,a@b@vac) for b in bosons] for a in bosons])
    joint=np.kron(vac,state)
    change=sum((np.kron(b,a) for b,a in zip(bosons,physical)),np.zeros((32,32),complex))
    measured=float(np.vdot(joint,change@change@joint).real)
    predicted=float(np.sum(WB*kernel).real)
    mean_part=float(np.sum(WB*np.outer(m,m)).real)
    errors=dict(commutator_dictionary=max(float(np.max(abs(Q(l)-a))) for l,a in zip(L,physical)),
        ordered_kernel=float(np.max(abs(kernel-direct_kernel))),
        full_square=float(abs(predicted-measured)),selfdual_purity=float(np.max(abs(P@P-P))))
    assert max(errors.values())<1e-12 and measured>1e-5 and mean_part>1e-7,(errors,measured,mean_part,m)
    return dict(round=815,status='working_not_formal',all_checks_passed=True,
        errors=errors,change_square=measured,conditional_fermion_mean_square_term=mean_part,
        omitted_mean_term_defect=mean_part,
        fermion_Gram_min_eigenvalue=float(np.linalg.eigvalsh(kernel).min()),
        original_neutral_mass_subblock_used=True,
        diagnostic_record='Admissible Nambu mixture of the two sterile spin modes, not the original fixed occupation record.',
        original_continuum_W_or_retention_time_computed=False,
        protected_record_classified=False,
        scope='Finite nonstationary coefficient calibration of the ordered record-change square; no new state supplied to the original continuum model.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
