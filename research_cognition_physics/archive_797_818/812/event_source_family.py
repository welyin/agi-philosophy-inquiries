"""812: common event-profile differentiation and retained source contacts.

Exact ten-dimensional BV diagnostic of the source identities. Not an evaluation
of continuum counterterms, an internal apparatus, or the Einstein PDE.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import factorial
import argparse,json,sys
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'811'))
import smooth_slice_connection as alg
p=alg.p;C=alg.C
TARGET=HERE/'event_source_family_results.json'
# Powers: epsilon, event u, old-source r, new-source z.
def zero_source(s,n):
    return {k:v for k,v in s.items() if k[n]==0}
def dexp(w,dw):
    out={};powers=[alg.ONE]
    for n in range(1,4):powers.append(alg.mul(powers[-1],w))
    for n in range(1,4):
        for k in range(n):
            term=alg.mul(alg.mul(powers[k],dw),powers[n-1-k])
            out=alg.add(out,alg.sc(term,Q(1,factorial(n))))
    return out
def run():
    X=alg.old.representative(p.base.X,Q(1,5))
    Y=alg.old.representative(p.base.Y,Q(1,7))
    Z=alg.old.representative(p.base.Z,Q(1,11))
    original={(1,0,0,0):p.scale(X,Q(1,3)),
              (1,0,1,0):p.scale(Z,Q(1,5))}
    event={(1,0,0,1):Y,(1,1,0,1):Z}
    contacts={(2,0,1,1):p.scale(X,Q(1,7)),(2,1,1,1):p.scale(Y,Q(1,13)),
              (2,0,0,2):p.scale(Z,Q(1,17)),(2,1,0,2):p.scale(X,Q(1,19))}
    w=alg.sc(alg.add(original,event,contacts),C(0,1))
    wold=alg.sc(original,C(0,1))
    sold=alg.exp(wold);s=alg.exp(w)
    relative=alg.mul(alg.adj(sold),s)
    assert alg.equal(zero_source(relative,3),alg.ONE)
    assert alg.equal(alg.mul(alg.adj(relative),relative),alg.ONE)
    assert all(not any(p.delta(v).flat) for v in relative.values())
    # Independent ordered Frechet sum, with all contact derivatives included.
    direct=alg.mul(alg.adj(sold),dexp(w,alg.diff(w,1)))
    assert alg.equal(direct,alg.diff(relative,1))
    obs=alg.sc(zero_source(alg.diff(relative,3),3),C(0,-1))
    assert alg.equal(obs,alg.adj(obs))
    assert alg.equal(alg.diff(alg.diff(obs,1),2),alg.diff(alg.diff(obs,2),1))
    # Incorrectly transporting only the classical event profile loses contacts.
    rawderivative=alg.diff(alg.sc(event,C(0,1)),1)
    wrongdu=alg.mul(alg.adj(sold),dexp(w,rawderivative))
    delta=alg.add(direct,alg.sc(wrongdu,-1))
    mixed=zero_source(zero_source(alg.diff(alg.diff(delta,3),2),3),2)
    twoevents=zero_source(zero_source(alg.diff(alg.diff(delta,3),3),3),2)
    assert mixed and twoevents
    # Source-independent moving slice frame; actual event derivative survives.
    u=alg.exp({(1,1,0,0):p.scale(X,C(0,1))})
    g=alg.sc(alg.mul(alg.diff(u,1),alg.adj(u)),C(0,-1))
    carried=alg.ad(u,obs)
    covariant=alg.add(alg.diff(carried,1),alg.sc(alg.comm(g,carried),C(0,-1)))
    assert alg.equal(covariant,alg.ad(u,alg.diff(obs,1)))
    assert alg.equal(alg.diff(carried,2),alg.ad(u,alg.diff(obs,2)))
    derivative=zero_source(alg.diff(obs,1),2)
    assert derivative
    return dict(round=812,all_checks_passed=True,epsilon_order=3,
        source_order=2,BV_dimension=10,old_source_family_exactly_preserved_at_new_source_zero=True,
        same_relative_source_unitarity_and_Ward=True,
        ordered_Frechet_derivative_matches_profile_derivative=True,
        profile_and_old_source_derivatives_commute=True,
        source_contact_derivatives_required=True,
        missed_mixed_old_new_source_contact=alg.first(mixed),
        missed_two_event_source_contact=alg.first(twoevents),
        moving_description_connection_removes_only_description_change=True,
        retained_actual_event_derivative=alg.first(derivative),
        original_continuous_quantum_counterterms_evaluated=False,
        original_profile_completed_in_declared_local_formal_class_by_argument=True,
        internal_material_record_carrier_or_finite_coupling_proven=False,
        new_numbered_test_groups=1)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true')
    args=parser.parse_args();r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

