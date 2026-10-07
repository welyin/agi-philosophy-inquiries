"""844: distinguish a Wick chart change from a changed physical state.

Finite paired CAR diagnostics test exact reset, nonGaussian higher records,
and the same complete instrument. They do not certify finite-coupling QFT.
"""
from pathlib import Path
from functools import lru_cache
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'physical_state_wick_instrument_results.json'
sys.path.insert(0,str(HERE.parent/'839'))
from gaussian_reference_entropy_split import kron_all,state

def run():
    ident=np.eye(16);x=np.array([[0,1],[1,0]],complex);y=np.array([[0,-1j],[1j,0]]);z=np.diag([1.,-1.])
    gs=[kron_all([z]*j+[p]+[np.eye(2)]*(3-j)) for j in range(4) for p in (x,y)]
    parity=kron_all([z]*4);logical=kron_all([z]*3+[np.eye(2)])
    rng=np.random.default_rng(844);raw=rng.normal(size=(8,8));raw=(raw-raw.T)*.17
    k=sum(1j*raw[i,j]*gs[i]@gs[j] for i in range(8) for j in range(i+1,8))
    ref0=state(.7*k);ref1=state(1.3*k)
    nu=np.einsum('aiaj->ij',ref1.reshape(8,2,8,2));g=np.kron(np.eye(8)/8,nu)
    covariance=lambda rho:np.array([[np.trace(rho@a@b) for b in gs] for a in gs])
    def factory(cov):
        @lru_cache(None)
        def normal(word):
            if not word:return {():1.+0j}
            first,rest=word[0],word[1:]
            out={(first,)+w:c for w,c in normal(rest).items()}
            for j,label in enumerate(rest):
                for w,c in normal(rest[:j]+rest[j+1:]).items():
                    out[w]=out.get(w,0)-(-1)**j*cov[first,label]*c
            return {w:c for w,c in out.items() if abs(c)>1e-15}
        return normal
    n0=factory(covariance(ref0));n1=factory(covariance(ref1))
    def matrix(poly):
        out=np.zeros((16,16),complex)
        for word,c in poly.items():
            a=ident.copy()
            for j in word:a=a@gs[j]
            out+=c*a
        return out
    # Re-expand one old-chart quartic Wick observable in the new chart.
    word=(0,1,4,5);physical_poly=n0(word);remaining=dict(physical_poly);new_coeff={}
    for w in sorted(physical_poly,key=lambda a:(-len(a),a)):
        c=remaining.get(w,0)
        if abs(c)<1e-14:continue
        new_coeff[w]=c
        for ww,cc in n1(w).items():remaining[ww]=remaining.get(ww,0)-c*cc
    assert max(abs(v) for v in remaining.values())<1e-13
    old_observable=matrix(physical_poly)
    transported=sum((c*matrix(n1(w)) for w,c in new_coeff.items()),np.zeros((16,16),complex))
    chart_error=float(np.linalg.norm(old_observable-transported));assert chart_error<1e-12
    wrongly_same_coeff=matrix(n1(word));omitted_terms=float(np.linalg.norm(old_observable-wrongly_same_coeff));assert omitted_terms>1e-3
    physical_state_change=float(np.trace((ref1-ref0)@old_observable).real)
    assert abs(physical_state_change)>1e-5
    # Full finite instrument fixed once for all states; no retuning of A.
    a=ident+.3j*gs[0]@gs[1]+.4*logical+.2j*gs[6]@gs[7]
    eta=.6;k0=np.linalg.inv(ident+1j*eta*a);k1=ident-k0
    completeness=float(np.linalg.norm(k0.conj().T@k0+k1.conj().T@k1-ident));assert completeness<1e-12
    source=.5j*gs[6]@gs[7]+.2j*gs[0]@gs[6]
    rows=[]
    for bias in (0.,.6):
        sigma=(np.eye(8)+bias*kron_all([z]*3))/8;rho=np.kron(sigma,nu)
        out=np.zeros_like(ref1);sumk=np.zeros_like(ref1);homogeneous=0.
        for aa in range(8):
            for bb in range(8):
                kab=np.zeros((8,8));kab[aa,bb]=np.sqrt(sigma[aa,aa].real)
                kab=np.kron(kab,np.eye(2));out+=kab@ref1@kab.conj().T;sumk+=kab.conj().T@kab
                sign=(-1)**((aa.bit_count()+bb.bit_count())%2)
                homogeneous=max(homogeneous,float(np.linalg.norm(parity@kab@parity-sign*kab)))
        reset_error=float(np.linalg.norm(out-rho));assert reset_error<1e-12
        assert np.linalg.norm(sumk-ident)<1e-12 and homogeneous<1e-12
        two=float(np.linalg.norm(covariance(rho)-covariance(g)));assert two<1e-12
        four=max(abs(np.trace((rho-g)@matrix({ids:1}))) for ids in itertools.combinations(range(8),4))
        assert four<1e-12
        branches=[kk@rho@kk.conj().T for kk in (k0,k1)]
        probs=[float(np.trace(br).real) for br in branches]
        assert min(probs)>0 and abs(sum(probs)-1)<1e-12
        chart_expectation_error=abs(np.trace(rho@(old_observable-transported)))
        rows.append(dict(record_bias=bias,reset_CP_residual=reset_error,
            parity_homogeneity_residual=homogeneous,two_point_gaussianization_error=two,
            all_four_distinct_Majorana_difference=float(four),six_Majorana_record=float(np.trace(rho@logical).real),
            same_instrument_probabilities=probs,
            same_physical_source_before=float(np.trace(rho@source).real),
            same_physical_source_after=float(np.trace(sum(branches)@source).real),
            chart_expectation_residual=float(chart_expectation_error)))
    assert abs(rows[0]['same_instrument_probabilities'][1]-rows[1]['same_instrument_probabilities'][1])>1e-3
    assert abs(rows[0]['same_physical_source_before']-rows[1]['same_physical_source_before'])<1e-12
    return dict(round=844,all_checks_passed=True,fresh_test_groups=1,
        finite_CAR_modes=4,full_pairing_and_cross_block_reference=True,
        Wick_chart_operator_residual=chart_error,omitting_Wick_lower_terms_error=omitted_terms,
        true_state_change_at_fixed_observable=physical_state_change,
        instrument_completeness_residual=completeness,record_and_state_rows=rows,
        gaussianization_preserves_entire_state=False,normal_order_change_is_heating=False,
        diagnostic_is_original_continuum=False,interacting_finite_coupling_entropy_computed=False,
        autonomous_preparation_proven=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();out=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
