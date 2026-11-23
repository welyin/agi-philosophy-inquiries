"""Reachability and observable row bounds for a bipartite three-port family."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def band(row, widths):
    center, *errors=row
    return center-widths*sum(errors), center+widths*sum(errors)


def bond(i,j):
    a=np.zeros((3,3))
    a[i,j]=a[j,i]=1
    return a


def gate(a,angle):
    d,W=np.linalg.eigh(a)
    return (W*np.exp(-1j*angle*d))@W.conj().T


def period(angles,generators):
    result=np.eye(3,dtype=complex)
    for n,angle in enumerate(angles):
        result=gate(generators[n%2], angle)@result
    return result


def lie_dimension(generators):
    basis=[]
    def add(matrix):
        v=np.concatenate([matrix.real.flatten(),matrix.imag.flatten()])
        for old,_ in basis:
            v-=old*(old@v)
        norm=np.linalg.norm(v)
        if norm>1e-9:
            v/=norm
            basis.append((v,v[:9].reshape(3,3)+1j*v[9:].reshape(3,3)))
    for a in generators:
        add(-1j*a)
    n=0
    while n<len(basis):
        current=basis[n][1]
        for _,other in basis[:n]:
            add(current@other-other@current)
        n+=1
    return len(basis)


def row_certificate(cfg):
    data=read(cfg['first_row'])
    cb=read(cfg['semileptonic'])['BelleII_Vcb']
    k=cfg['widths']
    ud=band(data['ud_beta'],k)
    product,form=band(data['K_product'],k),band(data['K_form'],k)
    us=(product[0]/form[1],product[1]/form[0])
    ratio,cb=band(data['ub_over_cb'],k),band(cb,k)
    ub=(ratio[0]*cb[0],ratio[1]*cb[1])
    box=[ud,us,ub]
    equal_gap=min(max(0,a[0]-b[1],b[0]-a[1]) for a,b in itertools.combinations(box,2))
    sum_gap=min(max(0,1-a[1]-b[1],a[0]+b[0]-1) for a,b in itertools.combinations(box,2))
    # Nonnegative entries and upper box endpoints overestimate any row fidelity.
    tops=np.array(box)[:,1]
    f_axis,f_hel=0.,0.
    for k in range(3):
        a=tops[k]
        b,c=np.delete(tops,k)
        f_axis=max(f_axis,math.hypot(a,(b+c)/math.sqrt(2)))
        f_hel=max(f_hel,(b+c)/2+math.hypot((b-c)/2,a/math.sqrt(2)))
    f=max(f_axis,f_hel)
    assert f<1 and sum_gap>.02 and equal_gap>.21
    return {'direct_marginal_moduli':box,'smallest_moduli_equality_gap':equal_gap,
            'smallest_sum_to_one_gap':sum_gap,'axis_overlap_upper':f_axis,
            'helicity_overlap_upper':f_hel,'minimum_optimally_rephased_unit_row_distance':math.sqrt(2*(1-f)),
            'minimum_projective_angle_deg':math.degrees(math.acos(f))}


def calculate():
    cfg=read('inputs.json')
    bonds=[bond(0,1),bond(1,2)]
    G=np.diag([1,1j,1])
    periods=[period(seq,bonds) for seq in cfg['pulses_rad']]
    bases=[]
    ev=[]
    axes=[]
    for U in periods:
        R=G.conj().T@U@G
        assert np.max(abs(R.imag))<1e-14 and abs(np.linalg.det(R)-1)<1e-14
        eig,V=np.linalg.eig(U)
        order=np.argsort(np.angle(eig))
        assert min(np.diff(np.angle(eig[order])))>.1
        ev.append(np.angle(eig[order]).tolist())
        V=V[:,order]
        np.testing.assert_allclose(V.conj().T@V,np.eye(3),atol=1e-13)
        bases.append(V)
        axes.append(V[:,1])
    K=bases[0].conj().T@bases[1]
    x=abs(axes[0].conj()@axes[1])
    y=math.sqrt((1-x*x)/2)
    low,high=(1-x)/2,(1+x)/2
    expected=np.array([[high,y,low],[y,x,y],[low,y,high]])
    assert min(np.max(abs(np.abs(K)-expected)),np.max(abs(np.abs(K)-expected[::-1])))<1e-12
    quartet=(K[0,0]*K[1,1]*K[0,1].conjugate()*K[1,0].conjugate()).imag
    assert abs(quartet)<1e-12
    a,b=cfg['oneway_test_rad']
    test=period([a,b],bonds)
    forward,back=abs(test[2,0])**2,abs(test[0,2])**2
    assert abs(forward-(math.sin(a)*math.sin(b))**2)<1e-14 and back<1e-25
    dims=[lie_dimension(bonds),lie_dimension([*bonds,np.diag([2,-1,-1])/3])]
    assert dims==[3,8]
    return {'row':row_certificate(cfg),'ordered_single_cycle_transition':[forward,back],
            'sample_quasiphases':ev,'sample_mixing_moduli':np.abs(K).tolist(),
            'sample_J_abs_below':1e-12,'control_lie_dimensions':dims}


def clean(x):
    if isinstance(x,dict):
        return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [clean(v) for v in x]
    return float(f'{x:.11g}') if isinstance(x,(float,np.floating)) else x


if __name__=='__main__':
    cli=argparse.ArgumentParser()
    cli.add_argument('--check',action='store_true')
    args=cli.parse_args()
    data=clean(calculate())
    path=HERE/'results.json'
    if args.check:
        assert data==json.loads(path.read_bytes())
        print('PASS: ordered exchange, spin-one overlap, full direct row box, Lie closure')
    else:
        with path.open('x',encoding='utf-8') as out:
            out.write(json.dumps(data,indent=2)+'\n')
