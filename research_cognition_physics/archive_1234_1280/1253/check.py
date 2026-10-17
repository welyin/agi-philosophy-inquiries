"""1253: original Material(1), constant swap coupling; read-only by default."""
from pathlib import Path
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
SOURCE = BASE / 'archive_1106_1119/1112/check.py'
EXPECTED = '6e49d71ec5535fd6999432481938663a5543e2704db993ec7a00b8f475b6cd9f'
RESULT = HERE / 'results.json'

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def close(x, label, tol=2e-8):
    v = float(np.linalg.norm(x))
    assert v < tol, (label, v)
    return v

def calculate():
    assert sha(SOURCE) == EXPECTED
    spec = importlib.util.spec_from_file_location('old1112', SOURCE)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    m = old.Material(1)
    n, shape = 864, (4,8,3,3,1,3)
    ident = np.eye(n, dtype=complex)
    basis = ident.reshape(shape+(n,))
    S = m.swap_external(basis).reshape(n,n)
    MD = m.h_action(basis).reshape(n,n)
    V = ident-S
    M = MD+ident+V
    mask = np.broadcast_to((np.arange(8)%2)[None,:,None,None,None,None], shape).reshape(n)
    Pi = np.diag(mask)
    algebra = {
        'swap_involution':close(S@S-ident, 'S2'),
        'swap_hermitian':close(S-S.conj().T, 'SH'),
        'commutator_mass':close(S@MD-MD@S, 'SM'),
        'actual_Mr_commutator':close(Pi@S-S@Pi, 'MrS'),
    }
    e, w = np.linalg.eigh(M)
    assert abs(e[0]-(m.baseline+1))<1e-10 and abs(e[-1]-(m.baseline+6))<1e-10
    algebra.update(mass_min=float(e[0]),mass_max=float(e[-1]),mass_width=float(e[-1]-e[0]))
    def evolve(s,B):
        return w @ (np.exp(-1j*s*e)[:,None]*(w.conj().T@B))
    delta, star = math.pi/4, math.pi/2
    J3 = m.loaded_joint(np.eye(3),delta).reshape(n,3)
    J = J3[:,1:3]
    target = np.exp(-1j*star)*m.blank_joint(np.eye(3),delta+star).reshape(n,3)
    task = {'all_qutrit_endpoint':close(evolve(star,J3)-target,'true R endpoint')}
    task['encoded_interaction_energy'] = close(J.conj().T@V@J-np.eye(2),'V code')
    task['encoded_mean_mass'] = close(J.conj().T@M@J-(20.5+2*math.sqrt(19))*np.eye(2),'M code')
    mean = 20.5+2*math.sqrt(19)
    task['encoded_mass_variance'] = close(J.conj().T@((M-mean*ident)@(M-mean*ident))@J-1.75*np.eye(2),'variance')
    weights = [1,3,4,4,3,1]
    spectral = []
    for j, count in enumerate(weights):
        eigmask = abs(e-(m.baseline+1+j))<1e-8
        C = w[:,eigmask].conj().T@J
        spectral.append(close(C.conj().T@C-(count/16)*np.eye(2),f'spectral {j}'))
    task['full_source_spectrum_residuals']=spectral
    ideal_bell = np.zeros(6,complex)
    ideal_bell[2]=ideal_bell[5]=1/math.sqrt(2)
    bell = np.outer(ideal_bell,ideal_bell.conj())
    erased = np.diag([.5,.5,0,0,0,0])
    rows=[]
    for s in [math.pi/4,star-4/75,star,star+4/75,7*math.pi/12]:
        B=evolve(s,J)
        clock=close(B.conj().T@Pi@B-math.sin((delta+s)/2)**6*np.eye(2),'clock')
        X=B.reshape(288,3,2)/math.sqrt(2)
        rr=np.einsum('ari,asj->risj',X,X.conj()).reshape(6,6)
        er=math.cos(s)**2
        err=close(rr-((1-er)*bell+er*erased),'R plus reference')
        # Compare complete isometries only after one common phase.
        centre=m.baseline+3.5
        full=float(np.linalg.norm(B-np.exp(-1j*centre*(s-star))*target[:,1:3],2))
        assert full <= 2.5*abs(s-star)+2e-9
        rows.append(dict(s=s,p=math.sin((delta+s)/2)**6,
                         clock_residual=clock,R_reference_residual=err,
                         erasure_probability=er,complete_isometry_distance=full))
    task['rows']=rows
    task['missing_source_expectation']=float(np.real(np.trace(J.conj().T@V@J)/2))
    assert abs(task['missing_source_expectation']-1)<1e-10
    bare= np.exp(-1j*star)*m.wait(m.loaded_joint(np.eye(2,3).T,delta),star).reshape(n,2)
    task['zero_coupling_R_nonvac_probability']=float(np.sum(abs(bare.reshape(288,3,2)[:,1:,:])**2)/2)
    assert task['zero_coupling_R_nonvac_probability']<1e-20

    eta, prob_error = F(1,100), F(1,50)
    time_error=F(8,3)*prob_error
    full_error=F(5,2)*time_error
    erase_error=time_error**2
    assert time_error==F(4,75) and full_error==F(2,15) and erase_error==F(16,5625)
    # A finite positive exponential sum proves 2 exp(-6)<1/200.
    exp_lower=sum((F(6)**j)/math.factorial(j) for j in range(19))
    assert exp_lower>400 and 2*30000*eta**2==6
    assert F(3,4)*F(11,12)**3>F(3,8)
    budget=dict(samples=30000,eta=str(eta),probability_error=str(prob_error),
                time_error=str(time_error),complete_output_bound=str(full_error),
                erasure_bound=str(erase_error),failure_bound='1/200',
                declared_composite_inventory_upper=32*30001,
                external_calibration_facilities_included=False)
    return dict(status='PASS',scientific_groups=3,source_sha256=EXPECTED,
                code_sha256=sha(Path(__file__)),algebra=algebra,task=task,budget=budget,
                scope='Given curved comoving composite, no dynamic geometry or external release.')

def compare(a,b):
    if isinstance(b,dict):
        assert a.keys()==b.keys()
        for k in b: compare(a[k],b[k])
    elif isinstance(b,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(b,float):
        assert math.isclose(a,b,rel_tol=2e-8,abs_tol=3e-9),(a,b)
    else:assert a==b,(a,b)

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--write-once',action='store_true')
    args=ap.parse_args()
    result=calculate()
    if args.write_once:
        with RESULT.open('x',encoding='utf-8') as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write('\n')
        print('FIRST_WRITE_PASS_1253')
    else:
        compare(result,json.loads(RESULT.read_text('utf-8-sig')))
        print('READ_ONLY_MATCH_1253')

