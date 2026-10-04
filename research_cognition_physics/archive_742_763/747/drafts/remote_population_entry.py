"""Entry747: original hopping sends unknown encoded occupation to another node.

The bare occupation response is not yet an actual remote T readout.
"""
from pathlib import Path
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
sys.path.insert(0,str(ROOT))
import joint_record_mass_feedback as old
TARGET=HERE/'remote_population_entry_results.json'
def run():
    rng=np.random.default_rng(747);rows=[]
    for _ in range(8):
        points=rng.normal(size=(2,5))*.7
        h=np.zeros((64,64),complex);d=np.zeros_like(h)
        for v in range(2):
            hs,ds=old.mass_x(points[v]);sl=slice(32*v,32*v+32);h[sl,sl]=hs;d[sl,sl]=ds
        # Full original species; SU(3),SU(2),U(1) gauge links are nontrivial.
        C=old.matter.gauge.group_exp(rng.normal(size=8),3)
        W=old.matter.gauge.group_exp(rng.normal(size=3),2)
        R=old.matter.representation(C,W,np.exp(1j*rng.normal()))
        spin=np.array([[.4,.2+.13j],[-.17+.11j,.31]],complex)
        transfer=R@np.kron(np.eye(16),spin)
        h[32:,:32]=transfer;h[:32,32:]=transfer.conj().T
        source=[{0:1.},{(1<<30)|(1<<31):1.}]
        outputs=[old.car.quadratic(s,h,d) for s in source]
        nB=lambda state:((state>>62)&1)+((state>>63)&1)
        coefficients=[sum(nB(k)*abs(z)**2 for k,z in out.items()) for out in outputs]
        difference=coefficients[1]-coefficients[0]
        tau=transfer[30:32,30:32];expected=float(np.trace(tau.conj().T@tau).real)
        assert abs(difference-expected)<2e-14
        common=2*abs(old.matter.Y['s'])**2*points[1,4]**2
        assert abs(coefficients[0]-common)<2e-14
        rows.append(dict(full_modes=64,empty_t2=coefficients[0],pair_t2=coefficients[1],
            original_hop_t2_difference=difference,sterile_transfer_HS_squared=expected,
            local_Majorana_creation_common_to_both=common))
    deps=('research_note_577.md','research_note_598.md','research_note_746.md',
          'joint_vertex_shared_evolution.py','joint_fermion_gauss_completion.py')
    return dict(entry_round=747,latest_completed_round=746,formal_test_count_unchanged=3430,
        rows=rows,full_graph_first_step_relation_is_analytic=True,
        no_remote_actual_T_instrument_claim=True,
        dependencies={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in deps},
        next='Determine the actual remote existing T readout response, keeping both scalar-edge and fermion-hop channels and no intermediate bare-occupation measurement.',
        all_checks_passed=True)
if __name__=='__main__':
    r=run()
    with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    print(json.dumps(dict(entry_round=747,rows=len(r['rows']),signal=r['rows'][0]['original_hop_t2_difference'])))
