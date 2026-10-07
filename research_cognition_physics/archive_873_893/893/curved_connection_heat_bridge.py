"""893: keep the actual nonflat horizontal connection in a joint heat process.
Full-species theorem uses finite-rank band occupation sectors and heat domination.
Numerics test original electron neutral-pair sector with TWO quantum coordinates;
the declared three-coordinate full-species theorem is not replaced by this grid.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'892'))
import multisource_transport_curvature as prior
old=prior.old
TARGET=HERE/'curved_connection_heat_bridge_results.json'
KAPPA=.001
HALFWIDTHS=(.025,.02,.02)
CENTER=-1/12

def band_trace_audit(m,q):
    rows=[]
    for k,a in (([0.,0.,0.],[0.,0.,0.]),([1.,2.,-1.],[.02,-.01,.015]),([0.,1.,0.],[-1/12,.012,-.009])):
        _,P,ds,_=prior.jets(m,q,np.array(a),np.array(k))
        traces=[abs(np.trace(P@prior.comm(ds[i],ds[j]))) for i,j in ((0,1),(0,2),(1,2))]
        # Sum is not the proof of individual cancellation; charged4 Clifford
        # blocks have zero trace separately, checked on every charge subspace.
        errors=[]
        for charge in sorted(set(q.tolist())):
            ids=np.flatnonzero(q==charge)
            Pc=P[np.ix_(ids,ids)]
            for i,j in ((0,1),(0,2),(1,2)):
                Fc=prior.comm(ds[i],ds[j])[np.ix_(ids,ids)]
                errors.append(abs(np.trace(Pc@Fc)))
        rows.append(dict(momentum=k,alpha=a,full_negative_band_trace=max(traces),
                         largest_charge_block_trace=max(errors)))
    assert max(r['largest_charge_block_trace'] for r in rows)<1e-10
    return rows

def log_heat(e,beta):
    e=np.asarray(e);e0=e.min()
    return float(-beta*e0+np.log(np.exp(-beta*(e-e0)).sum()))

def pair_grid(m,q,N,points):
    widths=HALFWIDTHS[:2]
    axes=[np.linspace(-d+2*d/(points+1),d-2*d/(points+1),points) for d in widths]
    coords=[np.array([CENTER+x,y,0.]) for x in axes[0] for y in axes[1]]
    weights=[KAPPA/(2*d/(points+1))**2 for d in widths]
    sites=len(coords);T=np.diag(np.full(sites,2*sum(weights)))
    HP=np.zeros((4*sites,4*sites),complex)
    scalar=np.zeros(sites);frames=[];hs=[]
    ids=[26,27,28,29]
    mass=np.exp(old.XI)*m.MASS[np.ix_(ids,ids)]
    _,V0=np.linalg.eigh(mass)
    pairbits=[b for b in range(16) if (b&3).bit_count()==1 and ((b>>2)&3).bit_count()==1]
    assert len(pairbits)==4
    for j,a in enumerate(coords):
        # Only first direction sees a cut; second/third shifts remain exact.
        p1=N//2-6*a[0]
        b=old.buffered(N,2*np.pi/N*p1)
        h=b*m.GAMMA[0][np.ix_(ids,ids)]-6*a[1]*m.GAMMA[1][np.ix_(ids,ids)]+mass
        energy=np.sqrt(np.trace(h@h).real/4)
        WF=prior.prior.exterior(prior.segment(h,mass)@V0)
        frames.append(WF);hs.append(h);scalar[j]=2*energy
        HP[4*j:4*j+4,4*j:4*j+4]=(2*sum(weights)+2*energy)*np.eye(4)
    leakage=0.;vacerr=0.;holonomies={}
    for i in range(points):
        for j in range(points):
            site=i*points+j
            for axis,(di,dj) in enumerate(((1,0),(0,1))):
                ni,nj=i+di,j+dj
                if ni>=points or nj>=points:continue
                other=ni*points+nj
                link=prior.prior.exterior(prior.segment(hs[site],hs[other]))
                C=frames[site].conj().T@link@frames[other]
                block=C[np.ix_(pairbits,pairbits)]
                remaining=[b for b in range(16) if b not in pairbits]
                leakage=max(leakage,float(np.max(abs(C[np.ix_(remaining,pairbits)]))))
                vacerr=max(vacerr,float(abs(C[3,3]-1)))
                ss=slice(4*site,4*site+4);oo=slice(4*other,4*other+4)
                HP[ss,oo]=-weights[axis]*block;HP[oo,ss]=-weights[axis]*block.conj().T
                T[site,other]=T[other,site]=-weights[axis]
                holonomies[site,other]=block;holonomies[other,site]=block.conj().T
    # Check a genuinely curved plaquette in the retained pair connection.
    maximum_loop=0.
    for i in range(points-1):
        for j in range(points-1):
            a=i*points+j;b=a+points;c=b+1;d=a+1
            loop=holonomies[a,d]@holonomies[d,c]@holonomies[c,b]@holonomies[b,a]
            maximum_loop=max(maximum_loop,float(np.linalg.norm(loop-np.eye(4))))
    HS=T+np.diag(scalar);ep=np.linalg.eigvalsh(HP);es=np.linalg.eigvalsh(HS)
    logratio=log_heat(ep,old.BETA)-np.log(4)-log_heat(es,old.BETA)
    assert leakage<1e-11 and vacerr<1e-11 and maximum_loop>1e-3
    assert ep[0]>=es[0]-1e-10 and logratio<=1e-10
    # The rank factor counts existing spin multiplicity, not extra coordinates.
    return dict(N=N,points_per_coordinate=points,quantum_coordinate_dimensions=2,
        neutral_pair_matrix_dimension=4*sites,neutral_pair_sector_leakage=leakage,
        filled_vacuum_connection_link_error=vacerr,
        maximum_neutral_pair_plaquette_holonomy_norm=maximum_loop,
        curved_pair_ground=float(ep[0]),scalar_comparison_ground=float(es[0]),
        log_heat_ratio_to_rank_times_scalar=logratio,
        heat_ratio_to_rank_times_scalar=float(np.exp(logratio)))

def run():
    m=old.load();_,matter=old.charges(m)
    q=prior.prior.prior.prior.last.em_charges(m,matter)
    traces=band_trace_audit(m,q)
    rows=[pair_grid(m,q,17,n) for n in (3,5,7)]
    return dict(round=893,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3678,
        argument_scope='Retaining892 nonflat horizontal connections, the declared three-coordinate transported/vacuum-subtracted free branch has invariant finite-rank band-occupation sectors. Vanishing filled-vacuum determinant curvature, Kato norm inequality and unitary heat-kernel domination yield full-species cut thermal suppression, a continuous trace-norm joint Gibbs limit, and compatible finite physical-CAR histories without a flat completion. This is not original bare Gauss or interacting Q/E equivalence.',
        declared_three_coordinate_box_halfwidths=list(HALFWIDTHS),kappa=KAPPA,
        original_full64_vacuum_curvature_checks=traces,
        original_two_coordinate_neutral_pair_heat_checks=rows,
        canonical_nonflat_connection_retained=True,
        extra_flat_completion_introduced=False,
        filled_vacuum_line_curvature_zero_on_declared_patch=True,
        full_original_species_matrix_sector_heat_bound=True,
        full_joint_continuous_trace_norm_bridge=True,
        compatible_finite_physical_CAR_histories=True,
        original_bare_Gauss_kinetic_recovered=False,
        original_vacuum_sources_matched=False,
        interacting_Q_E_or_gravity_bridge_complete=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args()
    r=run()
    if a.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n','utf-8')
    print(json.dumps(r,ensure_ascii=False,indent=2))
