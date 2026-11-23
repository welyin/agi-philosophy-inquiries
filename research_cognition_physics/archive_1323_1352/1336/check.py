"""Three-port constructive competitor for the finite direct-matching test."""
import argparse, json
from pathlib import Path
import numpy as np
I=np.eye(3)
def amplitude_basis(s12,s23,s13,delta):
    c12,c23,c13=np.sqrt(1-np.array([s12,s23,s13])**2)
    z=np.exp(1j*delta)
    return (np.array([[1,0,0],[0,c23,s23],[0,-s23,c23]]) @
            np.array([[c13,0,s13/z],[0,1,0],[-s13*z,0,c13]]) @
            np.array([[c12,s12,0],[-s12,c12,0],[0,0,1]]))

def model():
    s13=.00431
    c13=np.sqrt(1-s13**2)
    # Inputs at one scale; this is an existence witness.
    v=amplitude_basis(.2257/c13,.0416/c13,s13,1.2)
    targets=[np.array(x) for x in ([.00127,.619,171.7],[.00290,.055,2.89])]
    t=343.4
    M=30*t
    blocks=[]
    for masses,frame in zip(targets,[I,v]):
        singular=masses*np.sqrt(t*t+M*M-masses**2)/(t*M)
        blocks.append(np.diag(singular)@frame.conj().T)
    rest=I-sum(b.conj().T@b for b in blocks)
    lam,basis=np.linalg.eigh(rest)
    assert lam.min()>0
    third=(basis*np.sqrt(lam))@basis.conj().T
    w=np.vstack([*blocks,third])
    assert np.allclose(w.conj().T@w,I)
    compl=np.linalg.qr(w,mode="complete")[0][:,3:]
    left=[];right=[];reproduced=[]
    for j,target in enumerate(targets):
        p=np.eye(9)[:,3*j:3*j+3]
        kernel=np.column_stack([t*p,M*compl])
        ul,values,vr=np.linalg.svd(kernel)
        ids=np.argsort(values)[:3]
        assert np.allclose(values[ids],target,atol=1e-10,rtol=1e-8)
        reproduced.append(values[ids].tolist())
        left.append(ul[:,ids]);right.append(vr.conj().T[3:,ids])
    vl=left[0].conj().T@left[1]
    vr=right[0].conj().T@right[1]
    error=np.max(np.abs(np.abs(vl)-np.abs(v)))
    assert error<.0001 and np.abs(vl[1,2])>.041
    return dict(status="PASS", M_GeV=M, portal_mass_GeV=t,
         input_CKM_phase=1.2,masses_GeV=reproduced,left_current_abs=np.abs(vl).tolist(),
         right_current_norm=float(np.linalg.norm(vr,2)),max_entry_modulus_error=float(error),
         third_port_eigenvalues=lam.tolist(),third_port_unpolarized_probability=float(lam.mean()))

def compare(a,b):
    if isinstance(a,dict):
        assert a.keys()==b.keys()
        for k in a:compare(a[k],b[k])
    elif isinstance(a,str):assert a==b
    else:np.testing.assert_allclose(a,b,rtol=1e-8,atol=1e-10)

if __name__=="__main__":
    p=argparse.ArgumentParser(); group=p.add_mutually_exclusive_group(required=True)
    for opt in ('--check','--write-new'):group.add_argument(opt,action='store_true')
    args=p.parse_args(); result=model(); out=Path(__file__).with_name('results.json')
    if args.write_new:
        with out.open('x',encoding='utf8') as f:json.dump(result,f,indent=2)
    else:compare(result,json.loads(out.read_text(encoding='utf8')))
    print(json.dumps(result))