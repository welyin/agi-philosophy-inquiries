"""Enclose the comparator on the archived SK normal-ordering atmospheric grid."""
import argparse
import hashlib
import itertools
import json
import runpy
import tarfile
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
FLAVOUR=runpy.run_path(str(HERE.parent/'1669/independent_check.py'))


def read(name):
    return json.loads((HERE/name).read_bytes())


def verify_public_archive(path, inputs):
    with open(path,'rb') as stream:
        assert hashlib.sha256(stream.read()).hexdigest()==inputs['archive_sha256']
    with tarfile.open(path,'r:gz') as archive:
        data=archive.extractfile(inputs['member']).read()
    assert hashlib.sha256(data).hexdigest()==inputs['member_sha256']
    rows=np.loadtxt(data.splitlines())
    assert rows[:,4].min()==0
    table={tuple(r[:4]):r[4] for r in rows}
    for section in inputs['sections'].values():
        grid=np.array(section['q'])
        for idx,point in zip(np.ndindex(grid.shape),itertools.product(*section['axes'])):
            assert table[point]==grid[idx]
    return len(rows)


def value(axes,grid,points):
    # Repeated affine contractions of the sixteen public cell vertices.
    idx=[]
    alpha=[]
    for axis,point in zip(axes,points):
        assert np.all(point>=axis[0]-1e-12) and np.all(point<=axis[-1]+1e-12)
        k=np.clip(axis.searchsorted(point,side='right')-1,0,len(axis)-2)
        idx.append(k)
        alpha.append((point-axis[k])/(axis[k+1]-axis[k]))
    total=0.
    for bits in itertools.product((0,1),repeat=4):
        weight=1.
        for bit,u in zip(bits,alpha):
            weight*=u if bit else 1-u
        total+=weight*grid[tuple(k+bit for k,bit in zip(idx,bits))]
    return total


def masslimits(rs, accelerator):
    # d32/d21=R^2-1; both candidate expressions for each endpoint
    # increase with R. Intersection with the optional beam cut is monotone.
    low=np.maximum(7.232e-5*(rs**2-1),.002459*(1-1/rs**2))
    high=np.minimum(7.544e-5*(rs**2-1),.002563*(1-1/rs**2))
    if accelerator:
        low,high=np.maximum(low,.00237),np.minimum(high,.00251)
    return low,high


def audit(name,prior,section,beam):
    axes=[np.array(x) for x in section['axes']]
    grid=np.array(section['q'])
    np.testing.assert_array_equal(np.take(grid,0,axis=2),np.take(grid,-1,axis=2))
    # Every partial within each cell is a positive average of edge slopes.
    lips=[]
    for dim,axis in enumerate(axes):
        shape=[1]*4
        shape[dim]=-1
        lips.append(np.max(np.abs(np.diff(grid,axis=dim)/np.diff(axis).reshape(shape))))

    a,b=prior['with_JUNO_two_widths']
    # Pad published R endpoints for their printed precision, cap lower edge
    # with the analytically admissible mass endpoints by the mask below.
    cuts=np.linspace(a-2e-9,b+2e-9,4097)
    rs=(cuts[:-1]+cuts[1:])/2
    h=np.diff(cuts).max()/2
    (tl,th),(dl,dh)=FLAVOUR['coordinates'](name,rs,h)
    (t1,t2),(d1,d2)=FLAVOUR['coordinates'](name,rs)
    t=(t1+t2)/2
    delta=(d1+d2)/2
    _,p3=FLAVOUR['projectors'](name,rs)
    z=p3[:,0,0].real
    dt=np.maximum(th-t,t-tl)
    dd=np.maximum(dh-delta,delta-dl)
    eps=lips[1]*dt + lips[3]*(2*h+3e-11)+1e-8

    # Enclose the full parameter box in each R cell, and independently
    # evaluate a true feasible point (cell centre) for comparison.
    lo_o=masslimits(cuts[:-1],beam)[0]
    hi_o=masslimits(cuts[1:],beam)[1]
    lo_i,hi_i=masslimits(rs,beam)
    def extrema(lo,hi,phase_free):
        # For a fixed angle tuple every possible extrema along d32 lies at
        # a boundary or a grid knot. The same applies independently to a
        # free periodic phase.
        mask=(lo<=hi)
        phases=axes[2] if phase_free else [np.mod(delta,2*np.pi)]
        low=np.full(rs.shape,np.inf)
        high=np.full(rs.shape,-np.inf)
        for phase in phases:
            for knot in (lo,hi,*axes[0][1:-1]):
                m=np.clip(knot,lo,hi)
                q=value(axes,grid,[m,t,phase,z])
                low=np.minimum(low,q)
                high=np.maximum(high,q)
        return low[mask],high[mask],mask

    output={}
    for mode,free in [('frozen',False),('unrestricted_phase',True)]:
        li,ui,_=extrema(lo_i,hi_i,free)
        lo,up,mask=extrema(lo_o,hi_o,free)
        error=eps+(0 if free else lips[2]*dd)
        # Bounds include table rounding and numerical noise for both sides.
        margin=0.000051
        minimum=[float(min(lo-error[mask])-margin),float(min(li)+margin)]
        maximum=[float(max(ui)-margin),float(max(up+error[mask])+margin)]
        assert minimum[0]<minimum[1] and maximum[0]<maximum[1]
        output[mode]={'minimum_bracket':np.round(minimum,7).tolist(),
                      'maximum_bracket':np.round(maximum,7).tolist()}
    return output


def solve():
    data=read('inputs.json')
    old=read('../1668/results.json')['branches']
    result={}
    for beam,key in [(False,'two_reactor_gap_windows'),(True,'with_accelerator_mass_cut')]:
        output={name:audit(name,old[name],section,beam) for name,section in data['sections'].items()}
        third,second=output['third_leg_3_to_2'],output['second_leg_counterpart']
        margin=second['frozen']['minimum_bracket'][0]-third['frozen']['maximum_bracket'][1]
        assert margin>0
        output['minimum_score_separation_frozen']=round(margin,6)
        result[key]=output
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--archive',type=Path)
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    if args.archive:
        print('Public SHA-256 and all excerpted nodes verified. Full rows:',verify_public_archive(args.archive,read('inputs.json')))
    result=solve()
    if args.check:
        assert result==read('independent_results.json')
        print('1671 independent PASS: full-domain four-axis interpolation enclosures')
    else:
        with (HERE/'independent_results.json').open('x',encoding='utf8') as out:
            out.write(json.dumps(result,indent=2)+'\n')
