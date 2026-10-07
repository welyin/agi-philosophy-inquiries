"""810 working: the 809 nonselective readout and its covariance/source changes.

Original auxiliary-past 64-matrix block calibrates the CAR part. Bosonic data
are a canonical diagnostic, not the original E_Pi j or spacetime source.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'readout_source_probe_results.json'
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as old

def wick4(P,indices):
    W=P@old.CHARGE
    a,b,c,d=indices
    return W[a,b]*W[c,d]-W[a,c]*W[b,d]+W[a,d]*W[b,c]
def run():
    N=old.occupied(np.zeros(3));P=old.I-N;H=old.hamiltonian(np.zeros(3))
    f=np.zeros(64,complex);f[30]=1/np.sqrt(2);f[31]=1j/np.sqrt(2)
    cf=old.CHARGE@f.conj();F=np.outer(f,f.conj())+np.outer(cf,cf.conj())
    R=old.I-2*F
    assert np.max(abs(R@R-old.I))<1e-13
    reflected=R@P@R
    post=(P+reflected)/2;dN=-(post-P)
    assert np.linalg.eigvalsh(post).min()>-1e-12
    assert np.linalg.eigvalsh(old.I-post).min()>-1e-12
    preserved=float(abs(np.vdot(f,post@f)-np.vdot(f,P@f)))
    assert preserved<1e-13
    reality=float(np.max(abs(old.CHARGE@post.conj()@old.CHARGE+post-old.I)))
    assert reality<1e-13
    best=(0.,None)
    for i in range(32):
        for j in range(i+1,32):
            indices=[32+i,i,32+j,j]
            actual=(wick4(P,indices)+wick4(reflected,indices))/2
            Gaussian=wick4(post,indices)
            err=float(abs(actual-Gaussian))
            if err>best[0]:best=(err,[i,j])
    assert best[0]>1e-5
    source=[]
    for d in np.eye(5):
        dm=old.vertex.dmass(old.PHI,d)
        source.append(float((np.trace(dm@dN)/2).real))
    energy=float((np.trace(H@dN)/2).real);assert energy>1e-5
    # Nonselective sin readout: equal mixture of shifts +/- alpha E j/2.
    alpha=.7;V=np.array([[1.,.2],[.2,.8]]);h=np.array([0.,alpha/2])
    delta=np.outer(h,h);after=V+delta
    before_v=V[1,1];post_v=after[1,1]
    actual4=3*before_v**2+6*before_v*h[1]**2+h[1]**4
    excess=float(actual4-3*post_v**2);assert abs(excess+2*h[1]**4)<1e-13
    assert after[0,0]==V[0,0] and np.array_equal(delta[0],np.zeros(2))
    T=np.array([[2.,.3],[.3,1.5]])
    source_direct=(h@T@h)/2;source_covariance=np.trace(T@delta)/2
    assert abs(source_direct-source_covariance)<1e-14
    return dict(working_round=810,all_checks_passed=True,
        original_Nambu_dimension=64,original_auxiliary_momentum=[0,0,0],
        reflection_involution_residual=float(np.max(abs(R@R-old.I))),
        CAR_reality_residual=reality,measured_occupation_change=preserved,
        fermion_post_covariance_min=float(np.linalg.eigvalsh(post).min()),
        fermion_Gaussian_replacement_fourpoint_defect=best[0],defect_CAR_indices=best[1],
        original_five_mass_source_changes=source,
        diagnostic_auxiliary_free_energy_change=energy,
        boson_read_coordinate_variance_change=float(after[0,0]-V[0,0]),
        boson_positive_covariance_increment_eigenvalues=np.linalg.eigvalsh(delta).tolist(),
        boson_fourth_cumulant=excess,quadratic_source_change=float(source_direct),
        actual_original_coupled_boson_solution_used=False,
        full_spacetime_sources_or_autonomous_instrument_computed=False,
        post_state_generally_Gaussian=False,formal_round_completed=False,new_numbered_test_groups=0)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args()
    r=run()
    if args.write:
        assert not TARGET.exists(),'Do not overwrite evidence.'
        TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
