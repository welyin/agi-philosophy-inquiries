"""672 exploratory off-diagonal half-history kernel. No full Gauss claim."""
import json
import sys
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;sys.path.insert(0,str(BASE))
import joint_nonflat_mass_measure as prior
old=prior.old;internal=prior.internal;mass=prior.mass


def configuration(left,right,lam):
    links=np.tile(np.eye(16,dtype=complex),(2,4,1,1))
    links[0,:2]=left;links[0,2:]=right
    u,v,d,h=prior.kernel(links)
    count=4
    jm=np.kron(np.kron(np.eye(count),internal.VM),np.eye(16))
    jp=np.kron(np.kron(np.eye(count),internal.VP),np.eye(16))
    m=np.kron(np.eye(count),np.kron(internal.B,internal.T[0]))
    mb=np.kron(np.eye(count),np.kron(internal.B,internal.T[0].conj().T))
    pair,_=mass.mass_pairing(count,mass.car.PHI[0])
    q=prior.joint(u,v,d,jm,jp,m,mb,pair,lam=lam)
    reflect=np.eye(count)[[2,3,0,1]];tr=np.kron(reflect,np.eye(32))
    theta=np.block([[np.zeros_like(tr),tr],[tr,np.zeros_like(tr)]])
    idx=np.r_[np.arange(64,128),np.arange(192,256)]
    q['cross']=(theta@q['covariance'])[np.ix_(idx,idx)]
    q['gap']=float(min(abs(np.linalg.eigvalsh(h))))
    return q


def packed(z):return dict(real=z.real.tolist(),imag=z.imag.tolist())


def probe():
    scales=[0.,.10,.30]
    us=[prior.rep(*prior.group(67281,s)) for s in scales]
    rng=np.random.default_rng(67282)
    packets=rng.normal(size=(128,3))+1j*rng.normal(size=(128,3))
    packets/=np.linalg.norm(packets,axis=0)
    out=[]
    for lam in (0.,.37):
        scalar=np.zeros((3,3),complex);linear=np.zeros_like(scalar);gaps=[]
        for i in range(3):
            for j in range(3):
                q=configuration(us[i],us[j],lam)
                scalar[i,j]=q['weight'];linear[i,j]=q['weight']*np.vdot(packets[:,i],q['cross']@packets[:,j])
                gaps.append(q['gap'])
        scale=np.sqrt(abs(np.diag(scalar)))
        scalar/=scale[:,None]*scale[None,:];linear/=scale[:,None]*scale[None,:]
        out.append(dict(mass_scale=lam,minimum_Wilson_gap=min(gaps),
            scalar=packed(scalar),linear=packed(linear),
            scalar_hermiticity=old.err(scalar-scalar.conj().T),
            linear_hermiticity=old.err(linear-linear.conj().T),
            scalar_eigenvalues=np.linalg.eigvalsh((scalar+scalar.conj().T)/2).tolist(),
            linear_eigenvalues=np.linalg.eigvalsh((linear+linear.conj().T)/2).tolist()))
    return dict(round=672,exploratory=True,nt=2,nx=2,scales=scales,
        actual_internal_channels=16,shared_temporal_links_identity=True,
        bosonic_heat_kernel_and_Gauss_not_included=True,rows=out)


if __name__=='__main__':
    q=probe();target=HERE/'half_history_probe_results.json'
    if '--write-results' in sys.argv:
        with target.open('x',encoding='utf8') as f:json.dump(q,f,ensure_ascii=False,indent=2)
    else:assert q==json.loads(target.read_text('utf8'))
    print(json.dumps(q))
