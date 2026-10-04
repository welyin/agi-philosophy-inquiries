"""687 entry: exact full-group one-site holonomy and complete S9 moment."""
import json
import math
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_subgroup_measure_source as old
import joint_nonflat_mass_measure as gauge
TARGET=HERE/'full_holonomy_probe_results.json'


def rotation(r):
    out=[]
    for t in old.T:
        target=r.conj()@t@r.conj().T
        coeff=np.array([np.trace(v.conj().T@target)/16 for v in old.T])
        assert np.max(abs(coeff.imag))<2e-13
        assert np.linalg.norm(target-sum(c.real*v for c,v in zip(coeff,old.T)),2)<2e-13
        out.append(coeff.real)
    out=np.array(out).T
    assert np.linalg.norm(out.T@out-np.eye(10),2)<2e-13
    return out


def frames(r):
    eye=np.eye(16);a=(eye-r.conj().T)/2;b=(eye+r.conj().T)/2
    u=1j*(np.kron(a,old.VP)+np.kron(b,old.spin.GAMMA[3]@old.VP))
    v=1j*(np.kron(a,old.VM)+np.kron(b,old.spin.GAMMA[3]@old.VM))
    c=(r+r.conj().T)/2;s=(r-r.conj().T)/(2j)
    X=np.kron(c,np.eye(4))-1j*np.kron(s,old.spin.GAMMA[3])
    H=np.kron(eye,old.spin.G5)@X
    D=(np.eye(64)+X)/2
    return u,v,D,H


def moment(O,degree=8):
    A=(2*np.eye(10)+O+O.T)/4
    eig=np.linalg.eigvalsh(A)
    assert min(eig)>-2e-13 and max(eig)<1+2e-13
    coeff=[1.]
    for k in range(1,degree+1):
        coeff.append(float(sum(np.sum(eig**j)*coeff[k-j] for j in range(1,k+1))/(2*k)))
    rising=math.prod(range(5,5+degree))
    return math.factorial(degree)/rising*coeff[-1]


def weight(r):
    M=moment(rotation(r))
    W=float(abs(np.linalg.det((np.eye(16)+r)/2))**2)
    return W*M


def run():
    rng=np.random.default_rng(687)
    rows=[];rs=[]
    for j in range(7):
        triple=gauge.group(68700+j, .6+.13*j)
        r=gauge.rep(*triple);rs.append(r)
        assert np.linalg.norm(r-old.old.exterior(old.old.carrier(*triple)),2)<2e-13
        O=rotation(r)
        u,v,D,H=frames(r);g5=np.kron(np.eye(16),old.spin.G5)
        err=max(np.linalg.norm(H@H-np.eye(64),2),
                np.linalg.norm(u@u.conj().T-(np.eye(64)-H)/2,2),
                np.linalg.norm(v@v.conj().T-(np.eye(64)+H)/2,2),
                np.linalg.norm(u.conj().T@v,2))
        assert err<2e-13
        vals=[]
        for _ in range(4):
            E=rng.normal(size=10);E/=np.linalg.norm(E)
            t=sum(e*x for e,x in zip(E,old.T))
            A=u.T@np.kron(t,old.B)@u
            shifted=(E+O@E)/2
            expected=-np.kron(sum(e*x for e,x in zip(shifted,old.T)),old.EPS)
            pf=old.pfaffian(A)
            scalar=float((shifted@shifted)**8)
            assert np.linalg.norm(A-expected,2)<3e-13
            assert abs(pf-scalar)<2e-13
            vals.append(float(abs(pf-scalar)))
        bar=np.kron(np.eye(16),old.VP.conj().T)
        det=np.linalg.det(bar@D@v)
        W=abs(np.linalg.det((np.eye(16)+r)/2))**2
        assert abs(det-W)<2e-13
        M=moment(O)
        assert M>=25.**-8
        rows.append(dict(index=j,frame_error=float(err),Pfaffian_error=max(vals),
            Weyl_determinant_error=float(abs(det-W)),full_S9_factor=M,
            full_physical_weight=float(W*M)))
    # Representation law, noncommuting pairs, and all-group Gram with exact S9 moments.
    assert np.linalg.norm(rs[0]@rs[1]-rs[1]@rs[0],2)>.1
    assert np.linalg.norm(rotation(rs[0]@rs[1])-rotation(rs[0])@rotation(rs[1]),2)<5e-13
    gram=np.array([[weight(a.conj().T@b) for b in rs] for a in rs])
    assert np.max(abs(gram-gram.T))<2e-13
    ev=np.linalg.eigvalsh((gram+gram.T)/2)
    assert min(ev)>-1e-12 and np.max(abs(np.diag(gram)-1))<2e-13
    olderrors=[]
    for theta in (0.,.17,.31,np.pi/3,np.pi/2):
        r=gauge.rep(np.eye(3),np.eye(2),np.exp(1j*theta))
        olderrors.append(abs(moment(rotation(r))-old.measure(theta)[0]))
    assert max(olderrors)<2e-13
    return dict(entry_round=687,latest_formal_round=686,not_formal_round=True,
        original_16_channel_group_used=True,full_S9_moment_exact_not_sampled=True,
        rows=rows,Gram_eigenvalues=ev.tolist(),old615_Beta_formula_error=max(olderrors),
        analytic_S9_uniform_lower_bound=25.**-8,
        analytic_full_group_Haar_weight_lower_bound=2.**-40,
        physical_time_RP_or_complete_673_Hb_average_not_proved=True)


if __name__=='__main__':
    result=run()
    with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(dict(entry_round=687,all_checks_passed=True,
        smallest_Gram_eigenvalue=result['Gram_eigenvalues'][0],
        max_pf_error=max(r['Pfaffian_error'] for r in result['rows']))))

