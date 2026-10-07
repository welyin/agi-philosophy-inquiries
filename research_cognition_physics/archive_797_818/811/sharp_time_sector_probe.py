"""811 working: sharp-time TT--fermion principal-sector domain diagnostic.

Reuses the original 64-dimensional Gamma matrices. This is a homogeneous
principal-symbol calculation, NOT the constrained nonstationary PDE, not an
autonomous apparatus and not a proof about the entire original representation.
"""
from pathlib import Path
from itertools import product
import argparse, json, sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as old
TARGET=HERE/'sharp_time_sector_probe_results.json'
IDX=[30,31]
G=np.array([a[np.ix_(IDX,IDX)] for a in old.GAMMA])
I=np.eye(2)
def gamma(x):return np.einsum('i,ijk->jk',x,G)
def kernel(p,r):
    """Pair creation cross block; k=r-p. No energy delta at sharp time."""
    k=r-p;kn=np.linalg.norm(k);n=k/kn
    a=p+r-n*np.dot(n,p+r);a/=np.linalg.norm(a)
    b=np.cross(n,a)
    e=(np.outer(a,b)+np.outer(b,a))/np.sqrt(2)
    pp=(I+gamma(p/np.linalg.norm(p)))/2
    pm=(I-gamma(r/np.linalg.norm(r)))/2
    A=pp@gamma(e@(p+r))@pm/4
    square=float(np.vdot(A,A).real)
    energy=np.linalg.norm(p)+np.linalg.norm(r)+kn
    defects=max(np.linalg.norm(e@n),abs(np.trace(e)),abs(np.linalg.norm(e)-1))
    return square/(2*kn),square,energy,float(defects)

def integrate(scale,order,tau=0.):
    # A compact six-dimensional sector about p=e_x and r=e_y.
    # Quadrature certifies the code, not a positive continuum lower bound.
    centers=np.array([1.,0.,0.,0.,1.,0.]);halfwidth=.04
    x,w=np.polynomial.legendre.leggauss(order)
    result=0.;minimum=1e100;defect=0.
    for ix in product(range(order),repeat=6):
        z=scale*(centers+halfwidth*x[list(ix)])
        value,square,energy,d=kernel(z[:3],z[3:])
        weight=(scale*halfwidth)**6*np.prod(w[list(ix)])
        result+=weight*value*np.exp(-(tau*energy)**2)
        minimum=min(minimum,square/scale**2);defect=max(defect,d)
    return float(result),float(minimum),float(defect)

def run():
    point=kernel(np.array([1.,0,0]),np.array([0.,1,0]))
    assert abs(point[1]-1/32)<1e-14
    assert abs(point[0]-1/(64*np.sqrt(2)))<1e-14
    clifford=max(float(np.max(abs(G[i]@G[j]+G[j]@G[i]-2*(i==j)*I)))
                  for i in range(3) for j in range(3))
    scales=[1.,2.,4.,8.,16.,32.]
    rows=[]
    for scale in scales:
        bare,mn,d=integrate(scale,3)
        smooth,_,_=integrate(scale,3,.25)
        rows.append(dict(scale=scale,sharp_time_sector_integral=bare,
            normalized_by_scale_power_seven=bare/scale**7,
            Gaussian_time_average_sector_integral=smooth,
            sampled_minimum_pair_weight=mn,TT_residual=d))
    base=rows[0]['sharp_time_sector_integral']
    scaleerr=max(abs(r['normalized_by_scale_power_seven']/base-1) for r in rows)
    fine,_,_=integrate(1.,4)
    qerr=abs(fine-base)/fine
    assert clifford<1e-14 and scaleerr<1e-12 and qerr<1e-6
    assert max(r['TT_residual'] for r in rows)<1e-13
    assert rows[-1]['Gaussian_time_average_sector_integral']<1e-200
    # A lower-order mass vertex does not eliminate this leading symbol.
    return dict(working_round=811,all_checks_passed=True,
        original_Nambu_dimension=64,particle_sterile_indices=IDX,
        clifford_residual=clifford,point_pair_amplitude_squared=point[1],
        point_boson_weighted_kernel=point[0],sector_rows=rows,
        homogeneous_scaling_relative_residual=scaleerr,
        order_three_vs_four_quadrature_relative_difference=qerr,
        sharp_time_principal_sector_power=7,time_average_width=.25,
        interpretation='Positive TT--sterile pair-creation principal sector has a Lambda^7 sharp-time norm-square scaling. Fixed smooth temporal smearing suppresses its high-frequency shells.',
        actual_nonstationary_constrained_PDE_analyzed=False,
        full_original_Hamiltonian_nonexistence_proven=False,
        finite_graph_continuum_identification_proven=False,
        autonomous_internal_record_implemented=False,
        formal_round_completed=False,new_numbered_test_groups=0)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();result=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite existing evidence.'
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert result==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(result,ensure_ascii=False,indent=2))

