"""Identify radial Yukawa response contrasts from physical diagonal pseudoscalars."""
import argparse
import cmath
import json
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
EDGES=[(0,1),(1,0),(1,1),(2,1),(2,2),(1,2),(1,2)]


def assemble(values):
    out=np.zeros((3,3),complex)
    for (i,j),v in zip(EDGES,values,strict=True):
        out[i,j]+=v
    return out


def responses(edges,slopes):
    slopes=np.asarray(slopes)
    Y,excess=assemble(edges),assemble(slopes*edges)
    H=Y@Y.conj().T
    ul,m,vh=np.linalg.svd(Y)
    ix=np.argsort(m)
    ul,m,vr=ul[:,ix],m[ix],vh.conj().T[:,ix]
    g=ul.conj().T@(Y+excess)@vr
    p=(np.diag(g)/m).imag
    l=m*m
    ring=(H[0,1]*H[1,2]*H[2,0]).imag
    a,b,c,e,f,r,d=edges
    base1=r.real*d.imag
    base2=ring/H[0,0].real
    T=float(H[1,1].real+abs(f)**2)
    assert np.isclose(base2,c.real*e.real*f.real*d.imag,rtol=1e-10)
    moments=np.array([l@p,(l*(l-T))@p])
    contrasts=np.array([slopes[6]-slopes[5],slopes[6]+slopes[3]-slopes[2]-slopes[4]])
    recovered=moments/np.array([base1,base2])
    np.testing.assert_allclose(recovered,contrasts,rtol=1e-9,atol=1e-10)
    assert abs(p.sum())<1e-11
    # Compare the derivative with a finite variation of the actual named paths.
    t=1e-5
    def mass_at(u):
        return assemble(edges*np.exp((1+slopes)*u))
    num=(mass_at(t)-mass_at(-t))/(2*t)
    np.testing.assert_allclose(num,Y+excess,atol=2e-9)
    # The local readout is independent of chiral SVD phases once the masses are real positive.
    ratios=abs(g)**2
    return {'diagonal_p':p.tolist(),'diagonal_scalar':(np.diag(g)/m).real.tolist(),
            'two_contrasts':contrasts.tolist(),'measured_contrasts':recovered.tolist(),
            'bs_plus_sb_over_bb_massless_tree':(ratios[2,1]+ratios[1,2])/ratios[2,2],
            'bb_partial_width_ratio_massless_tree':ratios[2,2]/l[2],'y_squared':l.tolist(),
            'two_moments':moments.tolist(),'ring':ring}


def calculate():
    cfg=json.loads((HERE/'inputs.json').read_bytes())
    old=json.loads((HERE/cfg['quark_texture']).read_bytes())['reconstruction']
    a,b,c,e,f=old['normalized_real_edges_A_B_C_E_F']
    r,d=old['two_real_coefficients']
    phase=np.radians(old['eta_example_deg'])
    edges=np.array([a,b,c,e,f,r,d*cmath.exp(1j*phase)])
    np.testing.assert_allclose(edges[-2]+edges[-1],complex(*old['Z_real_imag']),atol=1e-14)
    P=np.zeros((7,6))
    for j,(i,k) in enumerate(EDGES):
        P[j,[i,k+3]]=1.
    C=np.array([[0,0,0,0,0,-1,1],[0,0,-1,1,-1,0,1]])
    assert np.linalg.matrix_rank(P)==5 and np.linalg.matrix_rank(C)==2 and not np.any(C@P)
    out={}
    for label,s in cfg['example_slopes'].items():
        data=responses(edges,s)
        other=responses(edges.conj(),s)
        np.testing.assert_allclose(data['diagonal_p'],-np.array(other['diagonal_p']),atol=1e-12)
        shifted=responses(edges,np.asarray(s)+P@np.array([.1,-.8,.45,.24,.6,-.21]))
        np.testing.assert_allclose(data['diagonal_p'],shifted['diagonal_p'],atol=1e-12)
        assert not np.allclose(data['diagonal_scalar'],shifted['diagonal_scalar'])
        out[label]=data
    return {'response_quotient_dimension':2,'invisible_dimension_for_pseudoscalars_only':5,'sample':out}


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,list):
        return [clean(v) for v in x]
    if isinstance(x,(float,np.floating)):
        return 0. if abs(x)<1e-14 else float(f'{x:.11g}')
    return x


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    data=clean(calculate())
    path=HERE/'results.json'
    if args.check:
        assert data==json.loads(path.read_bytes())
        print('PASS: radial derivatives, chiral mass basis, both cycle responses, kernel and independent rates')
    else:
        with path.open('x',encoding='utf-8') as out:
            out.write(json.dumps(data,indent=2)+'\n')
