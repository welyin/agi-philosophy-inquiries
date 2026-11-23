"""Full-polytope bounds on the heavy mixed and electromagnetic dimension-eight response."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads((HERE/path).read_bytes())


def strengths(ld, lt):
    d=sum(math.exp(4*l) for l in ld)
    k=sum(math.exp(4*l) for l in lt)
    return (k/18)/(d+k/27)


def compare(a, b):
    assert math.isclose(a,b,rel_tol=5e-12,abs_tol=1e-80)


def boundaries(v, mz, mass_d, mass_t, api):
    t,x,k=api['coefficients'](v)
    dmax=t-math.log(mass_d/mz)
    tmax=t-math.log(mass_t/mz)
    c=3*math.pi*(x-7*k)/8
    delta=3*math.pi*k
    # The entire triangular/quadrilateral domains in the report apply at these masses.
    smax=8*(tmax-c)/7
    bmax=2*dmax-delta
    assert 0<smax<dmax and 0<bmax<tmax and dmax<delta<2*dmax
    assert x > 2*tmax/(3*math.pi)
    # Exact extrema incl. their masses. Using A_lo means the lower mixed/photon ratio.
    cases={'A_lo':([smax,0],[tmax]),'A_hi':([smax/2]*2,[tmax]),
           'B_lo':([dmax,delta-dmax],[0,0]),'B_hi':([dmax,dmax],[bmax,0])}
    result={name:strengths(*z) for name,z in cases.items()}
    # A numerator correction S_max is convex and the triangle has a third vertex.
    assert result['A_lo']<=strengths([0,0],[c])<=result['A_hi']
    compare(result['A_lo'],1.5/(1+27*math.exp(-4*c)*(math.exp(smax/2)+math.exp(-7*smax/2))))
    compare(result['A_hi'],1.5/(1+54*math.exp(-4*c-1.5*smax)))
    smin=2*math.exp(-4*c-1.5*smax)
    result['A_deficit_from_1p5_at_A_hi']=40.5*smin/(1+27*smin)
    for rows in cases.values():
        example=api['example'](v,mz,*rows)
        assert example['common_trace_coefficient']>0
    U=mz*math.exp(t)
    result.update(A_log_min_D_over_T=(8*c-tmax)/7,B_log_min_T_over_D=delta-dmax,
                  B_lightest_charged_D_upper_GeV=U*math.exp(-delta/2),
                  B_T_lower_GeV=mass_d*math.exp(delta-dmax),
                  B_fourphoton_weight_lower_TeV_m4=2*(1000/U)**4*math.exp(2*delta),
                  A_T_upper_GeV=U*math.exp(-c))
    return result


def run():
    cfg=read('inputs.json')
    ew=read(cfg['gauge_inputs'])
    additional=read(cfg['strong_inputs'])
    api=runpy.run_path(str(HERE/cfg['gauge_function']))
    saved=read(cfg['published_points'])
    mz=ew['mu_ref_GeV']
    gs=[ew[name] for name in ('g1_GUT','g2','g3')]
    alpha,s_alpha=additional['strong_lattice_5f']['ALPHA_2026']
    mt,s_mt=additional['top_MS_matching_GeV']
    sources={}
    for key in ['common_MZ','independent_lattice']:
        low=[]
        for ns in itertools.product([-2,2],repeat=4):
            a=np.array([4*math.pi/(g+n*s)**2 for (g,s),n in zip(gs,ns)])
            if key=='independent_lattice':
                a[2]=1/(alpha+ns[2]*s_alpha)+math.log((mt+ns[3]*s_mt)/mz)/(3*math.pi)
            low.append(a)
        source={}
        for name,d,trip in [('original_window',mz,mz),
                            ('conditional_search_floors',cfg['conditional_D_mass_floor_GeV'],
                             cfg['conditional_T_mass_floor_GeV'])]:
            rows=[boundaries(v,mz,d,trip,api) for v in low]
            source[name]={k:[min(row[k] for row in rows),max(row[k] for row in rows)] for k in rows[0]}
        for kind in ['two_doublets_one_triplet','two_split_complete_pairs']:
            m=saved[key][kind]
            d=np.array(m['M_D_GeV'])/1000
            tr=np.array(m['M_T_GeV'])/1000
            D=sum(d**-4);T=sum(tr**-4)
            source[kind]=dict(P_TeV_m4=D+T/27,M_TeV_m4=T/18,M_over_P=T/18/(D+T/27),
                              deficit_from_1p5=40.5*(D/T)/(1+27*D/T))
        sources[key]=source
    for name in ['original_window','conditional_search_floors']:
        A=min(s[name]['A_lo'][0] for s in sources.values())
        B=max(s[name]['B_hi'][1] for s in sources.values())
        assert A>1.49 and B<1e-27
    return sources


def stable(x):
    if isinstance(x,dict):
        return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):
        return [stable(v) for v in x]
    if isinstance(x,float):
        return float(format(x,'.12g'))
    return x


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args=parser.parse_args()
    raw=json.dumps(stable(run()),indent=2,allow_nan=False)+'\n'
    out=HERE/'results.json'
    if args.check:
        assert json.loads(raw)==json.loads(out.read_bytes())
    else:
        out.write_bytes(raw.encode())
    print('1809', hashlib.sha256(raw.encode()).hexdigest()[:12])
