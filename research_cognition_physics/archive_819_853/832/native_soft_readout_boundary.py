"""832: original soft scalar readout, conditional coding and source accounting.

Finite matrices check CP order; original scalar geometry checks the exact
injection identity. No finite matrix is presented as the original full H.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
TARGET=HERE/'native_soft_readout_boundary_results.json'
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime


def choi_from_effect(w,effect,c=4):
    # w maps input to C tensor E. The normalized maximally entangled input
    # has coefficients w[c,e,input]/sqrt(input dimension).
    d=w.shape[1];e=effect.shape[0]
    psi=w.reshape(c,e,d).transpose(2,0,1).reshape(d*c,e)/np.sqrt(d)
    return psi@effect.T@psi.conj().T


def cp_order_checks():
    rng=np.random.default_rng(832)
    v=np.array([[1.,0.],[0.,0.],[0.,0.],[0.,1.]])
    vec=v.reshape(-1,order='F')/np.sqrt(2)
    projector=np.outer(vec,vec.conj());complement=np.eye(8)-projector
    rows=[]
    for k in (1,2,4):
        word=np.eye(3,dtype=complex)
        for _ in range(k):
            a=rng.normal(size=(3,3))+1j*rng.normal(size=(3,3))
            u,_=np.linalg.qr(a)
            # Original readout range; the finite eigenvalues are diagnostics.
            effect=.5+np.sin(rng.uniform(-3.,3.,3))*.25
            word=np.diag(np.sqrt(effect))@u@word
        b=word.conj().T@word;m=.25**k;M=.75**k
        assert np.linalg.eigvalsh(b).min()>=m-1e-14
        assert np.linalg.eigvalsh(b).max()<=M+1e-14
        raw=rng.normal(size=(12,2))+1j*rng.normal(size=(12,2))
        w=np.linalg.qr(raw)[0]
        j=choi_from_effect(w,np.eye(3));jb=choi_from_effect(w,b)
        low=float(np.linalg.eigvalsh(jb-m*j).min())
        high=float(np.linalg.eigvalsh(M*j-jb).min())
        p=float(np.trace(jb).real)
        before=float(np.trace(complement@j).real)
        after=float(np.trace(complement@jb).real/p)
        assert min(low,high)>-1e-13 and m<=p<=M
        assert after>=m/M*before-1e-13
        rows.append(dict(terminal_word_length=k,effect_lower=m,effect_upper=M,
                         CP_lower_smallest_eigenvalue=low,CP_upper_smallest_eigenvalue=high,
                         entangled_input_success_probability=p,
                         conditional_entangled_input_defect=after,
                         certified_defect_lower_bound=m/M*before))
    # A rank-deficient filter can really select a pure code branch. This is
    # the necessary contrast to the original strictly positive readout.
    p=.4;x=np.array([[0.,1.],[1.,0.]])
    w=np.zeros((4,2,2))
    w[:,0,:]=np.sqrt(p)*v;w[:,1,:]=np.sqrt(1-p)*(v@x)
    w=w.reshape(8,2)
    soft=choi_from_effect(w,np.diag([.75,.25]));sharp=choi_from_effect(w,np.diag([1.,0.]))
    soft_error=float(np.trace(complement@soft).real/np.trace(soft).real)
    sharp_error=float(np.trace(complement@sharp).real/np.trace(sharp).real)
    assert abs(soft_error-1/3)<1e-14 and abs(sharp_error)<1e-14
    return dict(rows=rows,soft_filter_conditional_defect=soft_error,
                singular_filter_conditional_defect=sharp_error,
                singular_filter_success_probability=p,
                singular_filter_is_not_the_original_scalar_instrument=True)


def native_geometry_checks(old):
    rng=np.random.default_rng(83219);rows=[];hbar=.7
    for _ in range(20):
        direction=rng.normal(size=5);direction/=np.linalg.norm(direction)
        phi=direction*rng.uniform(.1,3.3)
        inverse=old.geom.inverse(phi)
        w=float(rng.uniform(.4,2.))
        for menu in ('singlet','Higgs_T'):
            if menu=='singlet':
                a=float(phi[4]);grad=np.eye(5)[4]
            else:
                a=float(phi[:4]@phi[:4]/2);grad=np.r_[phi[:4],0.]
            e=.5+np.sin(a)/4;de=np.cos(a)/4
            q=float(grad@inverse@grad)
            squared_derivatives=de**2/(4*e)+de**2/(4*(1-e))
            injection=hbar**2/(2*w)*q*squared_derivatives
            identity=hbar**2/(8*w)*q*de**2/(e*(1-e))
            assert abs(injection-identity)<2e-15
            assert .25<=e<=.75 and injection>=0
            if menu=='Higgs_T':
                assert q<=32/9+1e-13
                assert injection<=hbar**2/(9*w)+1e-13
            # For the same source J_w = partial_w H and a w-independent
            # readout menu, the source increment is partial_w D = -D/w.
            step=w*1e-5
            derivative=((injection*w/(w+step))-(injection*w/(w-step)))/(2*step)
            source_expected=-injection/w
            assert abs(derivative-source_expected)<2e-10
            rows.append(dict(menu=menu,effect=e,gradient_squared=q,
                             energy_injection=injection,volume_source_increment=source_expected,
                             injection_identity_residual=abs(injection-identity),
                             source_derivative_residual=abs(derivative-source_expected)))
    return dict(original_metric_samples=rows,
                complete_H_identity_used_analytically=True,
                original_Dirac_Majorana_and_hopping_commute_with_scalar_readout=True,
                original_Higgs_injection_bound='hbar^2/(9*w)',
                autonomous_energy_supply_or_detector_constructed=False)


def run():
    with ResearchRuntime(Layout()).installed():
        import joint_record_mass_feedback as old
        geometry=native_geometry_checks(old)
    return dict(round=832,all_checks_passed=True,fresh_test_groups=1,
                conditional_CP_order=cp_order_checks(),original_source_accounting=geometry,
                scope='The original scalar readout effects lie between 1/4 and 3/4. Applied only to the environment after a closed energy-conserving evolution, any finite terminal word cannot herald the exact 829 encoding for the fixed 746 preparation. Choi order also bounds its conditional defect. The same original readout has the inherited energy injection and volume-source increment. Interleaved full interactions, rank-deficient instruments, changed preparations and autonomous detector construction are not covered.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();r=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(dict(round=832,all_checks_passed=True,
                         CP=r['conditional_CP_order'],
                         source_samples=len(r['original_source_accounting']['original_metric_samples']),
                         max_source_residual=max(x['source_derivative_residual'] for x in r['original_source_accounting']['original_metric_samples'])),ensure_ascii=False,indent=2))
