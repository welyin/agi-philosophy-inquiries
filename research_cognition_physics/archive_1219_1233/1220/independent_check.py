from __future__ import annotations
import argparse, hashlib, json, math, os
from pathlib import Path
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["MKL_NUM_THREADS"]="1"
os.environ["OMP_NUM_THREADS"]="1"
import numpy as np

HERE=Path(__file__).resolve().parent
PBOUND=14
NCUT=6
MASS=25.
SOURCE_MASS=20.
KAPPA=2.
COUPLING=.15

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def clean(x):
    if isinstance(x,dict): return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)): return [clean(v) for v in x]
    if isinstance(x,np.ndarray): return clean(x.tolist())
    if isinstance(x,(float,np.floating)): return float(format(float(x),".13g"))
    if isinstance(x,(bool,np.bool_)): return bool(x)
    if isinstance(x,(int,np.integer)): return int(x)
    return x

def basis(hin,pbound=PBOUND,ncut=NCUT):
    return [(int((hin-n)%2),p,a,n-a)
        for n in range(ncut+1) for a in range(n+1)
        for p in range(-pbound,pbound+1)]

def edges(s,kappa,g):
    h,p,a,b=s
    yield (h,p-1,a,b),-kappa/2,"source"
    yield (h,p+1,a,b),-kappa/2,"source"
    yield (1-h,p-1,a+1,b),g*math.sqrt(a+1),"field"
    yield (1-h,p+1,a,b+1),g*math.sqrt(b+1),"field"
    if a: yield (1-h,p+1,a-1,b),g*math.sqrt(a),"field"
    if b: yield (1-h,p-1,a,b-1),g*math.sqrt(b),"field"

def model(hin,gap,kappa=KAPPA,g=COUPLING,pbound=PBOUND,ncut=NCUT):
    states=basis(hin,pbound,ncut); ix={s:i for i,s in enumerate(states)}
    n=len(states); H=np.zeros((n,n)); Vs=np.zeros_like(H); Vf=np.zeros_like(H)
    outside={}
    for col,s in enumerate(states):
        h,p,a,b=s; q=-p-a+b
        H[col,col]=math.hypot(MASS+gap*h,p)-MASS+math.hypot(SOURCE_MASS,q)-SOURCE_MASS+a+b+kappa
        for dest,w,kind in edges(s,kappa,g):
            if not w: continue
            hd,pd,ad,bd=dest; qd=-pd-ad+bd
            assert p+q+a-b == pd+qd+ad-bd == 0
            if kind=="field": assert qd==q and hd==1-h
            else: assert (ad,bd,hd)==(a,b,h) and qd==q-(pd-p)
            assert (h+a+b)%2 == (hd+ad+bd)%2 == hin
            if dest in ix:
                row=ix[dest]; H[row,col]+=w
                (Vs if kind=="source" else Vf)[row,col]+=w
            else:
                outside.setdefault(dest,{})
                outside[dest][col]=outside[dest].get(col,0.)+w
    outstates=sorted(outside)
    R=np.zeros((len(outstates),n))
    for i,s in enumerate(outstates):
        for j,w in outside[s].items(): R[i,j]=w
    assert np.max(np.abs(H-H.T)) == 0
    Z=np.zeros((n,3))
    for j,p in enumerate((-1,0,1)): Z[ix[(hin,p,0,0)],j]=1
    st=np.array(states)
    h,p,a,b=st.T; q=-p-a+b
    bins=(h*6+(np.sign(b-a).astype(int)+1)*2+(q>0).astype(int))
    cross=bins[:,None] != bins[None,:]
    Vscross=Vs*cross
    assert np.max(np.abs(Vf*(~cross))) == 0
    # Cross-source edges form disjoint q=0 / q=1 two-dimensional blocks.
    assert np.all(np.count_nonzero(Vscross,axis=1)<=1)
    assert np.max(np.sum(np.abs(Vscross),axis=1)) <= kappa/2
    return dict(states=states,H=H,Vs=Vs,Vf=Vf,R=R,Z=Z,h=h,p=p,a=a,b=b,q=q,bins=bins,Vscross=Vscross,
                leaked_rows=len(outstates),merged_rows=sum(len(v)>1 for v in outside.values()))

ALPHA=np.array([.5j,1/math.sqrt(2),-.5j])
WEIGHTS=np.array([.25,.5,.25])

def expect(A,V,coherent=True):
    if coherent:
        psi=V@ALPHA
        if A.ndim==1: return float(np.vdot(psi,A*psi).real)
        return float(np.vdot(psi,A@psi).real)
    if A.ndim==1: return float(np.sum(WEIGHTS*np.sum(np.conj(V)*(A[:,None]*V),axis=0).real))
    return float(np.sum(WEIGHTS*np.sum(np.conj(V)*(A@V),axis=0).real))

def observations(m,V,hin,coherent=True):
    h,p,a,b,q=(m[x] for x in ("h","p","a","b","q"))
    prob=np.abs(V@ALPHA)**2 if coherent else np.abs(V)**2@WEIGHTS
    triple=(h!=hin)&(b>a)&(q>0)
    table=np.bincount(m["bins"],weights=prob,minlength=12).reshape(2,3,2)
    dsource=-expect(m["Vscross"],V,coherent)
    dfield=-expect(m["Vf"],V,coherent)
    means={key:float(prob@m[key]) for key in ("h","p","q","a","b")}
    mean_N=means["a"]+means["b"]
    assert abs(means["p"]+means["q"]+means["a"]-means["b"])<1e-10
    return dict(joint_event=float(prob[triple].sum()),
      hflip_and_plus=float(prob[(h!=hin)&(a>b)].sum()),
      hflip_and_minus=float(prob[(h!=hin)&(b>a)].sum()),
      hflip=float(prob[h!=hin].sum()),table_h_direction_sourcebin=table,
      means=means,shifted_energy=expect(m["H"],V,coherent),
      reference_pinching_mean_work=dict(source=dsource,field=dfield,total=dsource+dfield),
      reference_pinching_state_bound=KAPPA/2+2*COUPLING*math.sqrt(2*max(0,mean_N)),
      norm=float(prob.sum()))

def run_case(hin,gap,T,kappa=KAPPA,g=COUPLING,midpoints=0):
    m=model(hin,gap,kappa,g)
    H,Z,R=m["H"],m["Z"],m["R"]
    # A scalar shift is only a numerical rotating phase; energies use unshifted H.
    shift=float(np.trace(Z.T@H@Z)/3)
    K=H-shift*np.eye(len(H))
    lam,V=np.linalg.eigh(K); B=V.T@Z
    def evo(t): return V@(np.exp(-1j*t*lam)[:,None]*B)
    end=evo(T)
    co=observations(m,end,hin,True); inc=observations(m,end,hin,False)
    ini_co=observations(m,Z,hin,True); ini_inc=observations(m,Z,hin,False)
    hp=H*m["p"][None,:]-m["p"][:,None]*H
    hq=H*m["q"][None,:]-m["q"][:,None]*H
    force_p=expect(1j*hp,Z,True); force_q=expect(1j*hq,Z,True)
    force_inc=expect(1j*hp,Z,False)
    max_energy_error=max(abs(co["shifted_energy"]-ini_co["shifted_energy"]),abs(inc["shifted_energy"]-ini_inc["shifted_energy"]))
    assert abs(ini_co["shifted_energy"]-ini_inc["shifted_energy"])<1e-12
    assert abs(force_p+kappa/math.sqrt(2))<1e-12
    assert abs(force_q-kappa/math.sqrt(2))<1e-12
    assert abs(force_inc)<1e-12
    assert max_energy_error<1e-10
    # Diagnostic spectral errors, not an interval arithmetic certificate.
    gram=V.T@V-np.eye(len(H))
    residual=K@V-V*lam[None,:]
    out=dict(hin=hin,gap=gap,T=T,kappa=kappa,g=g,dimension=len(H),
      leakage_rows=m["leaked_rows"],merged_leakage_rows=m["merged_rows"],phase_shift=shift,
      coherent=co,incoherent=inc,
      coherent_minus_incoherent_joint=co["joint_event"]-inc["joint_event"],
      initial=dict(energy_coherent=ini_co["shifted_energy"],energy_incoherent=ini_inc["shifted_energy"],
                   p_force_coherent=force_p,q_force_coherent=force_q,p_force_incoherent=force_inc),
      diagnostic=dict(eigen_residual_fro=float(np.linalg.norm(residual)),eigen_orthogonality_fro=float(np.linalg.norm(gram)),
                      energy_conservation_error=max_energy_error))
    if midpoints:
        step=T/midpoints
        RV=R@V
        sum0=0.;sum1=0.
        for first in range(0,midpoints,128):
            ts=(np.arange(first,min(first+128,midpoints))+.5)*step
            phase=np.exp(-1j*lam[:,None]*ts[None,:])
            coeff=(phase[:,:,None]*B[:,None,:]).reshape(len(H),-1)
            vals=(RV@coeff).reshape(len(R),len(ts),3)
            deriv=(RV@(lam[:,None]*coeff)).reshape(len(R),len(ts),3)
            sum0+=float(np.sqrt(np.sum(np.abs(vals)**2,axis=(0,2))).sum())
            sum1+=float(np.sqrt(np.sum(np.abs(deriv)**2,axis=(0,2))).sum())
        rnorm=math.sqrt(float(np.abs(R).sum(axis=0).max())*float(np.abs(R).sum(axis=1).max()))
        C2=rnorm*float(np.linalg.norm(K@(K@Z)))
        main=step*sum0; correction=step*step/4*sum1
        remainder=T*step*step*C2/24
        out["residual_vector_taylor_diagnostic"]=dict(midpoints=midpoints,midpoint_integral=main,
          first_derivative_term=correction,second_derivative_remainder=remainder,
          total=main+correction+remainder,C2=C2,R_norm_bound=rnorm,
          status="Exact-arithmetic analytic bound evaluated in binary64; no outward arithmetic certification claimed.")
    return out

def dependencies():
    root=HERE.parents[1]
    names=["archive_1106_1119/1118/proof.md","archive_1106_1119/1119/proof.md",
           "archive_1159_1160/1160/proof.md","archive_1161_1200/1169/proof.md",
           "archive_1219_/1219/candidate_proof.md","archive_1219_/1220/admission.md"]
    return {name:sha(root/name) for name in names}

def recompute():
    train=[]; hold=[]
    for hin in (0,1):
        train.append(run_case(hin,.5,3*math.pi/4,midpoints=1024))
        hold.append(run_case(hin,2/3,5*math.pi/8,midpoints=1024))
    no_source=run_case(0,.5,3*math.pi/4,kappa=0)
    no_field=run_case(0,.5,3*math.pi/4,g=0)
    for c in train+hold:
        assert c["coherent_minus_incoherent_joint"]>0
        inc=c["incoherent"]
        assert abs(inc["hflip_and_plus"]-inc["hflip_and_minus"])<1e-10
        assert c["diagnostic"]["eigen_residual_fro"]<1e-9
        assert c["diagnostic"]["eigen_orthogonality_fro"]<1e-9
    assert abs(no_source["coherent"]["hflip_and_plus"]-no_source["coherent"]["hflip_and_minus"])<1e-10
    assert no_field["coherent"]["joint_event"]<1e-24
    return clean(dict(schema="1220-independent-numerical-audit-v1",status="PASS",all_checks_passed=True,check_groups=6,
       independent_code_sha256=sha(Path(__file__)),source_sha256=dependencies(),
       scope="Independent finite projection and residual-formula diagnostics. Analytic mother-domain and rigorous outward arithmetic must be audited separately.",
       checks=["basis, total momentum, parity, merged exterior rows","coherence preparation energy and force",
               "training joint event and complete endpoint accounting","frozen holdout joint event",
               "no-source and no-field controls","vector Taylor leakage formula and terminal 12-block reference"],
       training=train,holdout=hold,no_source=no_source,no_field=no_field))

def compare(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a: compare(a[k],b[k],path+"/"+k)
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+"/"+str(i))
    elif isinstance(a,float):
        assert math.isclose(a,b,rel_tol=3e-10,abs_tol=5e-12),(path,a,b)
    else: assert a==b,(path,a,b)

if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--write",action="store_true");args=ap.parse_args()
    target=HERE/"independent_results.json"
    if args.write and target.exists(): raise SystemExit("--write refuses to overwrite an existing result")
    result=recompute()
    if args.write:
        target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:
        compare(result,json.loads(target.read_text(encoding="utf-8")))
    print(json.dumps(dict(status="PASS",all_checks_passed=True,results_match=True,check_groups=6),ensure_ascii=False))
