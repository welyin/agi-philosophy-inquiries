"""Unnumbered 572 entry probe: the actual 571 source and CMC momentum balance.

Necessary condition for a flat conformal periodic seed, not a no-go for GR.
"""
import json
import sys
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_gauss_continuum_sampling as source


def run():
    N=16;d=3;L=source.LBOX;eps=L/N;volume=L**d
    points=np.stack(np.meshgrid(*([np.arange(N)*eps]*d),indexing='ij'),axis=-1)
    f=source.fields(points);step=2e-5
    derivatives=[]
    for mu in range(d):
        shift=np.zeros(d);shift[mu]=step
        plus,minus=source.fields(points+shift),source.fields(points-shift)
        derivatives.append({k:(plus[k]-minus[k])/(2*step) for k in ('X','s','a','a0')})
    covariant=np.zeros(points.shape);canonical=np.zeros_like(covariant)
    gh=[];gs=[]
    for mu in range(d):
        df=derivatives[mu]
        DX=df['X']+1j*np.einsum('...a,aij,...j->...i',f['a'][...,mu,:],source.T,f['X'])
        DX+=3j*f['a0'][...,mu,None]*f['X']
        covariant[...,mu]=np.sum(f['PX'].conj()*DX,axis=-1).real+f['Ps']*df['s']
        canonical[...,mu]=np.sum(f['PX'].conj()*df['X'],axis=-1).real+f['Ps']*df['s']
        for nu in range(d):
            F=df['a'][...,nu,:]-derivatives[nu]['a'][...,mu,:]-np.cross(f['a'][...,mu,:],f['a'][...,nu,:])
            FF=df['a0'][...,nu]-derivatives[nu]['a0'][...,mu]
            covariant[...,mu]+=np.sum(f['E'][...,nu,:]*F,axis=-1)+f['E0'][...,nu]*FF
            canonical[...,mu]+=np.sum(f['E'][...,nu,:]*df['a'][...,nu,:],axis=-1)+f['E0'][...,nu]*df['a0'][...,nu]
        gh.append(df['X'][...,1].real);gs.append(df['s'])
    Pc=np.sum(covariant,axis=(0,1,2))*eps**d
    Pcan=np.sum(canonical,axis=(0,1,2))*eps**d
    h=np.sqrt(source.PAR['h2'])
    expected=volume*np.array([.001*h+.0028,.001*h,0.])
    assert np.max(abs(Pc-expected))<3e-8 and np.max(abs(Pcan-expected))<3e-8
    gh=np.stack(gh,axis=-1);gs=np.stack(gs,axis=-1)
    H=gh.reshape(-1,3);S=gs.reshape(-1,3)
    gram=eps**d*(H.T@H+S.T@S)
    lam=np.linalg.pinv(gram,rcond=1e-12)@Pcan
    dp_h=-np.einsum('...i,i->...',gh,lam);dp_s=-np.einsum('...i,i->...',gs,lam)
    newPX=f['PX'].copy();newPX[...,1]+=dp_h
    qw,q0=source.charges(f['X'],f['PX']);newqw,newq0=source.charges(f['X'],newPX)
    charge_change=max(np.max(abs(qw-newqw)),np.max(abs(q0-newq0)))
    after=Pcan+eps**d*np.sum(dp_h[...,None]*gh+dp_s[...,None]*gs,axis=(0,1,2))
    distance2=float(eps**d*np.sum(dp_h**2+dp_s**2))
    assert charge_change<1e-14 and np.max(abs(after))<1e-12
    assert abs(distance2-Pcan@np.linalg.pinv(gram)@Pcan)<1e-12
    return dict(status='entry_probe_not_completed_round',source_round=571,dimension=d,N=N,
        total_covariant_momentum=Pc.tolist(),total_canonical_momentum=Pcan.tolist(),
        analytic_total_momentum=expected.tolist(),Gauss_integration_by_parts_error=float(np.max(abs(Pc-Pcan))),
        radial_counterflow_gram=gram.tolist(),radial_counterflow_coefficients=lam.tolist(),
        momentum_after_source_change=after.tolist(),gauge_charge_change=float(charge_change),
        counterflow_flat_kinetic_norm=float(np.sqrt(distance2)),
        scope='Flat conformal seed plus CMC on a periodic torus requires zero mean densitized momentum. The 571 source fails this route. Counterflow changes the physical source; no Einstein Hamiltonian constraint or general nonexistence is proved.')


if __name__=='__main__':
    result=run()
    with (HERE/'momentum_entry_probe_results.json').open('x',encoding='utf8',newline='\n') as stream:
        stream.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False))
