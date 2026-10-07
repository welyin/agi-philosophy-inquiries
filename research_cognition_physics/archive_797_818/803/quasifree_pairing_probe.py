"""803 working: ordered self-dual CAR contraction and a marginal selection rule.

Fock matrices independently check the covariance trace and both output changes.
The diagnostic state uses a constant four-mode subblock of the original mass;
it is NOT the fixed nonstationary continuum W or the coupled boson background.
"""
from pathlib import Path
import argparse, json, sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'802'))
import original_bff_vertex as vertex
TARGET=HERE/'quasifree_pairing_probe_results.json'

def car(n):
    result=[]
    for j in range(n):
        a=np.zeros((2**n,2**n),complex)
        for mask in range(2**n):
            if (mask>>j)&1:
                a[mask^(1<<j),mask]=(-1)**((mask&((1<<j)-1)).bit_count())
        result.append(a)
    return result

def bdg(h,d):
    return np.block([[h,d],[-d.conj(),-h.T]])

def run():
    rng=np.random.default_rng(803)
    n=4;ann=car(n);fields=ann+[a.conj().T for a in ann]
    ident=np.eye(2**n);one=np.eye(2*n)
    conjugation=np.block([[np.zeros((n,n)),np.eye(n)],[np.eye(n),np.zeros((n,n))]])
    def quadratic(k):
        return sum((fields[i].conj().T@fields[j]*k[i,j]/2
                    for i in range(2*n) for j in range(2*n)),np.zeros_like(ident,dtype=complex))
    selected=[24,25,30,31];indices=selected+[i+32 for i in selected]
    phi=np.array([0.,.55,0.,0.,.4])
    original_mass=vertex.mass(phi)[np.ix_(indices,indices)]
    energies,vectors=np.linalg.eigh(quadratic(original_mass))
    ground=vectors[:,0]
    assert energies[1]-energies[0]>1e-3
    def expect(a):return ground.conj()@a@ground
    pstate=np.array([[expect(a@b.conj().T) for b in fields] for a in fields])
    errors=dict(projector=float(np.max(abs(pstate@pstate-pstate))),
                self_duality=float(np.max(abs(conjugation@pstate.conj()@conjugation+ pstate-one))),
                ordered_covariance=0.,record_dictionary=0.,marginal_first_order=0.,joint_response=0.)
    values,vecs=np.linalg.eigh(pstate)
    blank=vecs[:,np.argmax(values)]
    native=np.eye(2*n,dtype=complex)[:,2]  # Original sterile annihilation mode.
    cases=[];wrong_order=[];wrong_weight=[]
    # Three oscillator levels suffice for these vacuum first-order contractions.
    a=np.diag(np.sqrt([1.,2.]),1);b=a+a.T
    phi_b=-1j*(a-a.T);rho_b=np.diag([1.,0.,0.])
    state=np.kron(rho_b,np.outer(ground,ground.conj()))
    def total_expect(x):return np.trace(state@x)
    for trial in range(10):
        def random_k():
            h=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));h=(h+h.conj().T)/10
            d=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));d=(d-d.T)/10
            return bdg(h,d)
        k1,k2=random_k(),random_k()
        k=k1+1j*k2
        assert np.max(abs(conjugation@k.T@conjugation+k))<1e-14
        for label,f in [('native_sterile',native),('P_adapted_blank',blank)]:
            cf=conjugation@f.conj()
            d=-np.outer(f,f.conj())+np.outer(cf,cf.conj())
            mode=sum((f[i].conjugate()*fields[i] for i in range(2*n)),np.zeros((16,16),complex))
            record=ident-mode.conj().T@mode
            errors['record_dictionary']=max(errors['record_dictionary'],float(np.max(abs(record-ident/2-quadratic(d)))))
            q=quadratic(k)
            covariance=expect(record@q)-expect(record)*expect(q)
            trace=np.trace((one-pstate)@d@pstate@k)/2
            errors['ordered_covariance']=max(errors['ordered_covariance'],float(abs(covariance-trace)))
            if label=='native_sterile':
                wrong_order.append(float(abs(covariance-np.trace((one-pstate)@k@pstate@d)/2)))
                wrong_weight.append(float(abs(covariance-2*trace)))
            interaction=np.kron(b,quadratic(k1))
            field=np.kron(phi_b,ident);pointer=np.kron(np.eye(3),record)
            dp=1j*(pointer@interaction-interaction@pointer)
            df=1j*(field@interaction-interaction@field)
            marginal=total_expect(dp)
            joint=total_expect(df@pointer+field@dp)-total_expect(df)*total_expect(pointer)-total_expect(field)*marginal
            # W_B(Phi,B)=-i, so K_f=-i K1, including the Nambu half in Q.
            formula=-np.trace((one-pstate)@d@pstate@(-1j*k1)).imag
            errors['marginal_first_order']=max(errors['marginal_first_order'],float(abs(marginal)))
            errors['joint_response']=max(errors['joint_response'],float(abs(joint-formula)))
            if trial==0:
                cases.append(dict(record=label,blank_probability=float(expect(record).real),
                    off_diagonal_record_norm=float(np.linalg.norm((one-pstate)@d@pstate)),
                    first_order_marginal=float(marginal.real),
                    first_order_joint=float(joint.real),ordered_complex_covariance=[float(covariance.real),float(covariance.imag)]))
            if label=='P_adapted_blank':assert abs(joint)<1e-12 and abs(covariance)<1e-12
    assert max(errors.values())<3e-12,errors
    assert max(wrong_order)>1e-3 and max(wrong_weight)>1e-3
    assert abs(cases[0]['first_order_joint'])>1e-3
    return dict(working_round=803,all_checks_passed=True,samples=20,maximum_residuals=errors,
        examples=cases,wrong_order_maximum_defect=max(wrong_order),
        missing_Nambu_half_maximum_defect=max(wrong_weight),
        original_mass_subblock_indices=selected,
        original_continuum_state_evaluated=False,
        calibration_vertices_are_original_continuum_vertices=False,
        universal_blank_preparation_not_assumed=True,
        formal_round_completed=False,new_numbered_test_groups=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite saved evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
