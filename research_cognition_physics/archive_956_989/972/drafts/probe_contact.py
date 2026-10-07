from pathlib import Path
import importlib.util, json, math
import numpy as np
STAGE=Path(__file__).resolve().parents[2]
s=importlib.util.spec_from_file_location('old968',STAGE/'968/internal_relay.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
old,h,q,f,w,ident=m.material()
def entropy(r):
    p=np.linalg.eigvalsh((r+r.conj().T)/2);p=p[p>1e-14]
    return float(-p@np.log(p))
N=48;om=.5;g=.003;beta=2*math.log(2)
a=np.diag(np.sqrt(np.arange(1,N)),1);num=np.diag(np.arange(N))
pr=2.**(-np.arange(N));pr/=sum(pr);tau=np.diag(pr)
rf0=np.kron((np.outer(w[:,0],w[:,0])+np.outer(w[:,1],w[:,1]))/2,tau)
for eta in (0.,.05):
    hm=h+eta*f;HM=np.kron(hm,np.eye(N));HB=np.kron(np.eye(4),om*num)
    VI=g*np.kron(q,a+a.T)+g*g/om*np.kron(q@q,np.eye(N))
    H=HM+HB+VI;ev,V=np.linalg.eigh(H)
    for t in (100.,1000.,10000.):
        U=(V*np.exp(-1j*t*ev))@V.T
        reduced=[];full=[]
        for k in (0,1):
            rho=U@np.kron(np.outer(w[:,k],w[:,k]),tau)@U.conj().T
            rr=rho.reshape(4,N,4,N)
            reduced.append(np.einsum('anbn->ab',rr));full.append(rho)
        rho=(full[0]+full[1])/2;sm=(reduced[0]+reduced[1])/2
        rb=np.einsum('aman->mn',rho.reshape(4,N,4,N))
        info=entropy(sm)-(entropy(reduced[0])+entropy(reduced[1]))/2
        td=sum(abs(np.linalg.eigvalsh(reduced[0]-reduced[1])))/2
        Q=float(np.trace((rho-rf0)@HB).real)
        changes=[float(np.trace((rho-rf0)@A).real) for A in (HM,HB,VI)]
        print(json.dumps(dict(eta=eta,t=t,I_record=info,loss=math.log(2)-info,
          trace_distance=td,heat=Q,energy_changes=changes,S=entropy(sm)),ensure_ascii=False))
