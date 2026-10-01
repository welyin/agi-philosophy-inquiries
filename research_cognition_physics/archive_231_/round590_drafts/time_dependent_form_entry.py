"""590 entry: explicit geometric work coefficients for the SAME 589 source.

This evaluates a frozen classical Gauss phase point along a background path.
It does not propagate a quantum state or solve dynamic gravity.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import joint_full_spatial_metric as model
TARGET=HERE/'time_dependent_form_entry_results.json'


def weighted(A,dA,p,dp,inverse=False):
    if inverse:
        r=1/p;dr=-dp/p**2
    else:r,dr=p,dp
    C=r[..., :,None]*A*r[...,None,:]
    dC=dr[..., :,None]*A*r[...,None,:]+r[..., :,None]*dA*r[...,None,:]+r[..., :,None]*A*dr[...,None,:]
    return C,dC


def coefficients(q,t):
    x,y,z=np.moveaxis(q['grid'],-1,0)
    psi0=1.1+.05*np.cos(x)+.02*np.sin(y); a=.05+.025*np.cos(x)
    psi=psi0*np.exp(a*np.sin(1.3*t)); dp=psi*a*1.3*np.cos(1.3*t)
    initial=model.shape_field(q['grid']); eig,U=np.linalg.eigh(initial)
    S=(U*np.log(eig)[...,None,:])@np.swapaxes(U,-1,-2)
    b=.3*np.sin(.8*t);db=.24*np.cos(.8*t)
    A=model.shape_exp(S,b);dA=db*S@A; inv=np.linalg.inv(A);dinv=-inv@dA@inv
    pe=[];dpe=[];pf=[];dpf=[]
    for mu in range(3):
        pe.append((psi+np.roll(psi,-1,axis=mu))/2)
        dpe.append((dp+np.roll(dp,-1,axis=mu))/2)
    for mu,nu in model.PAIRS:
        def face(v):return (v+np.roll(v,-1,axis=mu)+np.roll(v,-1,axis=nu)+np.roll(np.roll(v,-1,axis=mu),-1,axis=nu))/4
        pf.append(face(psi));dpf.append(face(dp))
    pe,dpe,pf,dpf=[np.stack(v,axis=-1) for v in (pe,dpe,pf,dpf)]
    B=model.gram(A)
    dB=np.stack([np.stack([dinv[...,a,c]*inv[...,b,d]+inv[...,a,c]*dinv[...,b,d]-dinv[...,a,d]*inv[...,b,c]-inv[...,a,d]*dinv[...,b,c]
                          for c,d in model.PAIRS],axis=-1) for a,b in model.PAIRS],axis=-2)
    cg,dg=weighted(inv,dinv,pe,dpe)
    ce,de=weighted(A,dA,pe,dpe,True)
    cm,dm=weighted(B,dB,pf,dpf,True)
    w=q['eps']**3*psi**6; dw=6*w*dp/psi
    return dict(w=w,gradient=cg,electric=ce,magnetic=cm),dict(w=dw,gradient=dg,electric=de,magnetic=dm)


def generalized(A,B):
    eig,U=np.linalg.eigh(A)
    inv=(U*eig[...,None,:]**-.5)@np.swapaxes(U,-1,-2)
    return np.linalg.eigvalsh(inv@B@inv)


def work_bound(C,dC):
    bounds=[float(np.max(abs(dC['w']/C['w'])))]
    bounds.extend(float(np.max(abs(generalized(C[k],dC[k])))) for k in ('gradient','electric','magnetic'))
    return max(bounds)


def source_value(q,raw,C,dC=None):
    phi=q['phi']; P=q['P']; eps=q['eps']
    node=.5*np.einsum('...a,...ab,...b->...',P,model.original.inverse(phi),P)
    U=model.original.node_potential(phi)
    if dC is None:
        scalar=float(np.sum(node/C['w'])); onsite=float(np.sum(C['w']*U)); coef=C
    else:
        scalar=float(np.sum(-node*dC['w']/C['w']**2));onsite=float(np.sum(dC['w']*U));coef=dC
    grad=float(eps/2*np.sum(coef['gradient']*raw['grad']))
    electric=float(np.sum(coef['electric']*raw['electric'])/eps)
    mag=float((np.sum(coef['magnetic']*raw['odd'])+np.sum(np.diagonal(coef['magnetic'],axis1=-2,axis2=-1)*np.diagonal(raw['even'],axis1=-2,axis2=-1)))/eps)
    return dict(scalar_kinetic=scalar,onsite=onsite,gradient=grad,electric=electric,magnetic=mag,total=scalar+onsite+grad+electric+mag)


def run():
    q=model.original.shared_source(8);raw=model.graph_data(q,np.ones(q['phi'].shape[:-1]))
    base,_=coefficients(q,0.);baseE=source_value(q,raw,base)['total'];rows=[]
    for t in (0.,.23,.71,1.25):
        C,dC=coefficients(q,t); E=source_value(q,raw,C);dE=source_value(q,raw,C,dC);lam=work_bound(C,dC)
        spectra=[C['w']/base['w'],base['w']/C['w']]
        spectra.extend(generalized(base[k],C[k]) for k in ('gradient','electric','magnetic'))
        lo=min(float(v.min()) for v in spectra);hi=max(float(v.max()) for v in spectra)
        assert lo*baseE-1e-12<=E['total']<=hi*baseE+1e-12
        assert abs(dE['total'])<=lam*E['total']+1e-12
        errors=[]
        for h in (.004,.002,.001):
            cp,_=coefficients(q,t+h);cm,_=coefficients(q,t-h)
            plus,minus=source_value(q,raw,cp),source_value(q,raw,cm)
            errors.append(max(abs((plus[k]-minus[k])/(2*h)-dE[k]) for k in E))
        assert errors[-1]<errors[0]/12 and errors[-1]<1e-4
        rows.append(dict(time=t,energy=E,geometric_work_rate=dE,relative_form_rate_bound=lam,
                         base_form_comparison=[lo,hi],centered_difference_errors=errors))
    T=1.25;c0,_=coefficients(q,0.);c1,_=coefficients(q,T)
    exact=source_value(q,raw,c1)['total']-source_value(q,raw,c0)['total'];quad=[]
    for n in (8,16):
        x,w=np.polynomial.legendre.leggauss(n);v=0.
        for u,z in zip(x,w):
            C,dC=coefficients(q,T*(u+1)/2);v+=z*T/2*source_value(q,raw,C,dC)['total']
        quad.append(dict(nodes=n,integrated_fixed_source_work=float(v),endpoint_error=abs(v-exact)))
    assert quad[-1]['endpoint_error']<1e-10 and abs(exact)>1e-3
    return dict(status='590 time-dependent form entry; complete round pending',checks_passed=True,
                explicit_path_rows=rows,frozen_source_endpoint_energy_change=exact,work_quadrature=quad,
                propagated_quantum_state=False,Einstein_geometry_not_solved=True,
                sources={n:hashlib.sha256((HERE.parent/n).read_bytes()).hexdigest() for n in ('joint_full_spatial_metric.py','research_note_589.md')})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write-results',action='store_true');args=parser.parse_args()
    result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False,indent=2))
