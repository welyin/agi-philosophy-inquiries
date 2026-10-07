"""928: conditional natural source allocation and a finite quantum boundary witness.
The stochastic matrix acts on an operator list, not on the quantum state.
"""
from pathlib import Path
from fractions import Fraction as F
import itertools,argparse,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent

def compositions(total):
    yield (total,)
    for first in range(1,total):
        for rest in compositions(total-first):yield (first,)+rest

def surjections(n):
    for m in range(1,n+1):
        for labels in itertools.product(range(m),repeat=n):
            if set(labels)==set(range(m)):yield m,labels

def sparse_rank(rows):
    piv={}
    for source in rows:
        row={j:F(v) for j,v in source.items() if v}
        while row:
            j=min(row)
            if j not in piv:
                scale=row[j];piv[j]={k:v/scale for k,v in row.items()};break
            scale=row[j]
            for k,v in piv[j].items():
                x=row.get(k,F(0))-scale*v
                if x:row[k]=x
                elif k in row:del row[k]
    return len(piv)

def constraints(weighted):
    objects=list(compositions(4)) if weighted else list(range(1,5))
    dims={w:len(w) if weighted else w for w in objects};offset={};size=0
    for w in objects:offset[w]=size;size+=dims[w]**2
    def ix(w,i,j):return offset[w]+i*dims[w]+j
    rows=[];morphisms=0
    for w in objects:
        n=dims[w]
        for j in range(n):rows.append({ix(w,i,j):1 for i in range(n)})
        for m,labels in surjections(n):
            v=tuple(sum(w[j] for j in range(n) if labels[j]==a) for a in range(m)) if weighted else m
            assert v in offset;morphisms+=1
            for a in range(m):
                for j in range(n):
                    row={ix(w,i,j):1 for i in range(n) if labels[i]==a}
                    k=ix(v,a,labels[j]);row[k]=row.get(k,0)-1
                    rows.append({k:x for k,x in row.items() if x})
    rank=sparse_rank(rows)
    out=dict(objects=len(objects),variables=size,morphisms=morphisms,exact_rank=rank,affine_dimension=size-rank)
    if weighted:
        delta=[F(0)]*size
        for w in objects:
            for i in range(len(w)):
                for j in range(len(w)):delta[ix(w,i,j)]=F(w[i],4)-int(i==j)
        residual=max(abs(sum(F(v)*delta[j] for j,v in row.items())) for row in rows)
        assert residual==0 and size-rank==1
        isolated={ix((1,1,1,1),0,3):1}
        out.update(candidate_null_vector_residual=float(residual),rank_with_one_isolated_pair=sparse_rank(rows+[isolated]))
        assert out['rank_with_one_isolated_pair']==size
    else:assert rank==size
    return out

def natural(w,b):
    w=np.asarray(w,float);w=w/w.sum()
    return (1-b)*np.eye(len(w))+b*np.outer(w,np.ones(len(w)))
def matrix(m,labels):
    c=np.zeros((m,len(labels)))
    for j,a in enumerate(labels):c[a,j]=1
    return c

def run():
    weighted=constraints(True);unweighted=constraints(False)
    w=np.array([.1,.2,.3,.4]);b=.4;maxerr=0.;count=0
    for m,lab in surjections(4):
        C=matrix(m,lab);err=np.linalg.norm(C@natural(w,b)-natural(C@w,b)@C,ord=1)
        maxerr=max(maxerr,float(err));count+=1
    assert maxerr<1e-14
    # Equal cell counts are not retained by unequal regrouping unless recorded as weights.
    C3=np.array([[1,0,0],[0,1,1]],float)
    forget=C3@natural([1,1,1],.5)-natural([1,1],.5)@C3
    assert abs(np.max(np.abs(forget))-1/12)<1e-14
    # The finite four-unit grid alone does not prove the universal b<=1 bound.
    finite_b=1.2
    finite_min=min(float(natural(w,finite_b).min()) for w in compositions(4))
    assert finite_min>=-1e-14
    smaller_min=float(natural([.01,.99],finite_b).min());assert smaller_min<0
    # One adjacent finite source allocation; same coarse totals do not determine boundary response.
    alpha=.25
    L=np.zeros((4,4))
    for i in range(3):L[i,i]+=1;L[i+1,i+1]+=1;L[i,i+1]-=1;L[i+1,i]-=1
    T=np.eye(4)-alpha*L
    C=np.array([[1,1,1,0],[0,0,0,1]],float)
    assert T.min()>=0 and np.max(abs(T.sum(axis=0)-1))<1e-14
    e1=np.eye(4)[:,0];e3=np.eye(4)[:,2]
    coarse_error=float(np.linalg.norm(C@e1-C@e3))
    response_gap=C@T@(e1-e3)
    assert coarse_error==0 and np.max(abs(response_gap-np.array([.25,-.25])))<1e-14
    # Keep the distinct boundary contribution. Positive 3-entry interface suffices for this source task.
    D=np.array([[1,1,0,0],[0,0,1,0],[0,0,0,1]],float)
    B=np.array([[1,.75,.25],[0,.25,.75]],float)
    boundary_error=float(np.linalg.norm(C@T-B@D));assert boundary_error==0 and B.min()>=0
    # Homogeneous disjoint pairs are a weaker allowed-grouping contract, not a universal one.
    pair=np.array([[.75,.25,0,0],[.25,.75,0,0],[0,0,.75,.25],[0,0,.25,.75]])
    Cp=np.array([[1,1,0,0],[0,0,1,1]],float)
    assert np.linalg.norm(Cp@pair-Cp)==0
    # Quantum witness: four excitation modes plus an internal two-level probe, all finite.
    energies=np.array([[(j>>(3-i))&1 for j in range(16)] for i in range(4)],float)
    H0=energies.sum(axis=0);SB=T[3]@energies
    nc=np.array([0.,1.]);omega=1.;eta=1.;time=4*np.pi
    Htot=np.kron(H0,np.ones(2))+omega*np.tile(nc,16)+eta*np.kron(SB,nc)
    phase=np.exp(-1j*time*Htot);plus=np.array([1,1],complex)/np.sqrt(2)
    X=np.array([[0,1],[1,0]],complex);effect=np.kron(np.eye(16),(np.eye(2)+X)/2)
    probs=[];energy_changes=[];initial_energies=[];coarse_totals=[]
    for j in (8,2):
        state=np.kron(np.eye(16)[:,j],plus);final=phase*state
        probs.append(float(np.vdot(final,effect@final).real))
        e=float(np.sum(abs(state)**2*Htot));ef=float(np.sum(abs(final)**2*Htot))
        initial_energies.append(e);energy_changes.append(ef-e)
        coarse_totals.append((C@energies[:,j]).tolist())
    assert max(abs(probs[0]-1),abs(probs[1]))<1e-14
    assert coarse_totals[0]==coarse_totals[1]==[1.,0.]
    assert abs(initial_energies[0]-1.5)<1e-14 and abs(initial_energies[1]-1.625)<1e-14
    assert max(abs(x) for x in energy_changes)<1e-14
    # Unknown coherent data remain in the full system; the probe is entangled, not a classical replacement.
    data=(np.eye(16)[:,8]+np.eye(16)[:,2])/np.sqrt(2)
    state=phase*np.kron(data,plus);amplitude=state.reshape(16,2)
    probe=amplitude.conj().T@amplitude
    purity=float(np.trace(probe@probe).real);assert abs(purity-.5)<1e-14
    return dict(round=928,date='2026-10-06',all_scientific_checks_passed=True,
        weighted_family_exact_calibration=weighted,unweighted_family_exact_calibration=unweighted,
        weighted_naturality_checks=count,weighted_naturality_max_error=maxerr,
        forgotten_weight_max_entry_defect=float(np.max(np.abs(forget))),
        finite_grid_b_above_one_positive_minimum=finite_min,smaller_weight_b_above_one_minimum=smaller_min,
        local_counterexample=dict(coarse_total_difference=coarse_error,coarse_source_difference=response_gap.tolist(),
            boundary_interface_error=boundary_error,restricted_pair_grouping_error=float(np.linalg.norm(Cp@pair-Cp))),
        autonomous_quantum_source_probe=dict(dimension=32,minimum_total_energy=float(Htot.min()),
            duration=time,coarse_initial_energies=coarse_totals,full_initial_energies=initial_energies,
            final_plus_probabilities=probs,total_energy_changes=energy_changes,coherent_probe_purity=purity),
        six_protocol_gap_map_covers_C01_to_C27=True,all_regrouping_summary_is_extra_assumption=True,
        primitive_operator_assignment_and_weights_remain_inputs=True,
        source_allocation_not_quantum_state_channel=True,all_six_protocols_jointly_realized=False,
        tamper_resistance_and_repeated_reset_constructed=False,locality_derived_from_protocols=False,
        unique_geometric_stress_derived=False,spacetime_or_GR_generated=False,
        physical_minimum_scale_or_exact_continuum_assumed=False,whole_stage_completed=False,full_goal_completed=False,
        source_hashes={str(Path('928')/Path(__file__).name):hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run()
    path=HERE/'source_grouping_selection_results.json'
    if a.write:
        assert not path.exists();path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in out.items() if k!='source_hashes'},ensure_ascii=False,indent=2))
