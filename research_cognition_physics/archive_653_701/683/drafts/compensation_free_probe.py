"""Exact free-background obstruction to the declared scalar reflection of K*K."""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_local_source_lift as body
base=body.base;TARGET=HERE/'compensation_free_probe_results.json'

def run():
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    _,_,_,h,_=base.kernel(links);n=len(h)
    g5=np.kron(np.eye(4),np.kron(base.internal.spin.G5,np.eye(16)))
    gamma0=np.kron(np.eye(4),np.kron(base.internal.spin.GAMMA[3],np.eye(16)))
    swap=np.array([[0.,1.],[1.,0.]])
    qt=np.kron(np.kron(np.array([[0.,1.],[-1.,0.]]),np.eye(2)),np.eye(64))
    W=np.eye(n)-np.kron(np.kron(np.eye(2),swap),np.eye(64))
    A=g5@gamma0@qt
    P=np.kron(np.kron(swap,np.eye(2)),np.eye(64));pm=np.kron(np.eye(2),P)
    assert body.err(A@A-np.eye(n))==0 and body.err(W@A-A@W)==0
    assert body.err(h-g5@W-A)==0
    rows=[]
    for a in (.2,.1):
        k,_,_=body.blocks(g5@h,g5,a,1);G=k.conj().T@k
        delta=pm.T@G@pm-G
        block=8*a*(np.eye(n)+a*W)@A
        predicted=np.block([[block,np.zeros((n,n))],[np.zeros((n,n)),-block]])
        assert body.err(delta-predicted)<1e-13
        defect=float(np.linalg.norm(delta,2));bound=8*a*(1+2*a)
        assert abs(defect-bound)<2e-13
        rows.append(dict(a=a,formula_error=body.err(delta-predicted),
            exact_defect_norm=bound,numerical_defect_norm=defect,
            minimum_Gram_eigenvalue=float(np.linalg.eigvalsh(G)[0])))
    return dict(original_678_L1_Gram=True,original_free_2_by_2_box=True,rows=rows,
        claim='Declared componentwise scalar reflection fails; original physical RP not refuted.')

if __name__=='__main__':
    data=run()
    if TARGET.exists():assert data==json.loads(TARGET.read_text('utf8'))
    else:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(data))
