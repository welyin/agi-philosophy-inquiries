"""850: exact cubic projection for arbitrary ten-CAR quadratic responses.
The candidate response matrices here calibrate the algorithm, not the original
geometric response image or an original continuous spacetime solution.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse,itertools,json,random,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
STAGE=ROOT/'research_cognition_physics/archive_764_'
TARGET=HERE/'geometry_cubic_code_probe_results.json'
sys.path.insert(0,str(STAGE/'829'))
import majorana_code_source_bridge as code

def run():
    gamma,_,_,_,comp,_,_=code.code_data()
    logical=code.mul((0,0,1),code.product(gamma[i] for i in (0,1,4,5,12,14)))
    support=json.loads((STAGE/'845/record_source_hbar_filtration_results.json').read_text('utf-8'))['chosen_logical_axis_sixth_moment_support']
    assert len(support)==80
    pairs=list(itertools.combinations(range(20),2))
    ops={ij:code.mul((0,0,1),code.mul(gamma[ij[0]],gamma[ij[1]])) for ij in pairs}
    def direct(a,b,c):
        out=F(0)
        for (ij,x),(kl,y),(mn,z) in itertools.product(a.items(),b.items(),c.items()):
            op=code.mul(ops[ij],code.mul(ops[kl],ops[mn]))
            kind,phase=comp(code.mul(logical,op))
            if kind=='scalar':
                assert phase in (0,2)
                out+=x*y*z*(1 if phase==0 else -1)
        return out
    def pf(a,ids):
        if not ids:return F(1)
        return sum(((-1)**(j+1))*a.get((ids[0],ids[j]),F(0))*pf(a,ids[1:j]+ids[j+1:])
                   for j in range(1,len(ids)))
    def cubic(a):return -6*sum(row['sign']*pf(a,tuple(row['indices'])) for row in support)
    def add(*aa):
        out={}
        for a in aa:
            for ij,v in a.items():out[ij]=out.get(ij,F(0))+v
        return {ij:v for ij,v in out.items() if v}
    def polar(a,b,c):
        return (cubic(add(a,b,c))-cubic(add(a,b))-cubic(add(a,c))-cubic(add(b,c))
                +cubic(a)+cubic(b)+cubic(c))/6
    selected=[{(0,1):F(1)},{(4,5):F(1)},{(12,14):F(-1)}]
    assert direct(*selected)==polar(*selected)==1
    assert all(cubic(a)==0 for a in selected)
    assert cubic(add(*selected))==6
    cancelling=add(*selected,{(13,15):F(1)})
    assert cubic(cancelling)==direct(cancelling,cancelling,cancelling)==0
    # A 28-dimensional two-block response space is blind to all 80 directions.
    two_block_pairs=list(itertools.combinations(range(8),2))
    assert len(two_block_pairs)==28
    assert not any(set(row['indices'])<=set(range(8)) for row in support)
    blind=[{(0,1):F(1)},{(0,2):F(1)},{(0,3):F(1)}]
    assert all(polar(*[blind[i] for i in ids])==0 for ids in itertools.combinations_with_replacement(range(3),3))
    rng=random.Random(850);rows=[]
    for trial in range(4):
        mats=[]
        for index in range(3):
            a=dict(selected[index])
            # Include general cross-block matrix entries, not just local Pauli axes.
            for pair in rng.sample(pairs,12):a[pair]=F(rng.choice((-3,-2,-1,1,2,3)),rng.choice((1,2,3)))
            mats.append(a)
        d=direct(*mats);p=polar(*mats)
        assert d==p
        assert all(direct(*[mats[j] for j in perm])==d for perm in itertools.permutations(range(3)))
        combined=add(*mats)
        assert cubic(combined)==direct(combined,combined,combined)
        rows.append(dict(trial=trial,trilinear_coefficient=str(d),diagonal_cubic_coefficient=str(cubic(combined)),
                         exact_direct_and_Pfaffian_agreement=True))
    return dict(round=850,all_checks_passed=True,fresh_test_groups=1,
        exact_rational_and_original_Pauli_arithmetic=True,Majoranas=20,quadratic_coordinates=190,
        original_logical_six_word_support=80,Pfaffian_size=6,
        canonical_three_distinct_response_coefficient='1',canonical_combined_cubic_coefficient='6',
        each_canonical_direction_separately_has_zero_cubic=True,
        nonzero_entries_can_cancel_the_logical_cubic=True,
        blind_two_spatial_block_response_space_dimension=28,
        equally_three_dimensional_blind_and_detecting_response_spaces_exist=True,
        response_order_permutation_independence_verified=True,rows=rows,
        original_geometric_response_image_computed=False,
        original_pure_metric_cubic_nonzero_proven=False,
        original_mean_geometry_nonzero_proven=False,
        physical_quadratic_kernel_to_code_map_is_analytic=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
