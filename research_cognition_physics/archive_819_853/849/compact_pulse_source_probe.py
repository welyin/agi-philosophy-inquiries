"""849 working probe: a compact zero-moment pulse retains a native source.

This is the original constant auxiliary Dirac symbol. It does not construct a
pure-metric observable on the original varying, constrained background.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'compact_pulse_source_probe_results.json'
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original

def run():
    width=.2
    sx=np.array([[0,1],[1,0]],complex)
    sy=np.array([[0,-1j],[1j,0]],complex);sz=np.diag([1.,-1.])
    rotation=(np.eye(2)+1j*sx)/np.sqrt(2)
    assert np.max(abs(rotation.conj().T@sy@rotation-sz))<1e-14
    def quad(n):
        s,w=np.polynomial.legendre.leggauss(n);u=1-s*s;eta=np.exp(-1/u)
        norm=np.dot(w,eta)
        first=-2*s/u**2*eta
        second=(4*s*s/u**4-2/u**2-8*s*s/u**3)*eta
        return s,w*eta/norm,w*first/norm,w*second/norm
    def average(k,n):
        s,base,first,weights=quad(n)
        eig,v=np.linalg.eigh(original.hamiltonian([0,0,k]))
        multiplier=np.einsum('t,tij->ij',weights,
            np.cos(width*s[:,None,None]*(eig[:,None]-eig[None,:])[None,:,:]/k))
        full=v@((v.conj().T@original.GAMMA[1]@v)*multiplier)@v.conj().T
        return full[30:32,30:32]
    s,base,first,second=quad(256)
    c=float(base@np.cos(2*width*s));expected=-4*width**2*c
    coefficient=float(second@np.cos(2*width*s))
    moments=[float(second@(s**j)) for j in (0,1)]
    antiderivative_integral=float(first.sum())
    assert max(abs(v) for v in moments+[antiderivative_integral])<1e-11
    assert abs(coefficient-expected)<1e-11 and abs(coefficient)>.1
    # Zero total area cancels a commuting Z source, but not the rotating Y.
    assert abs(second.sum())<1e-11
    rows=[];quad_error=0.
    for k in (2.,8.,32.,128.):
        a=average(k,256);rot=rotation.conj().T@a@rotation
        z=float(np.trace(rot@sz).real/2);y=float(np.trace(a@sy).real/2)
        triple=.54*z*z*y
        rows.append(dict(momentum=k,rotated_Z_coefficient=z,unrotated_Y_coefficient=y,
                         native_triple_source_content_difference=triple,
                         native_compression_error_from_principal=float(np.max(abs(a-expected*sy)))))
        if k==8.:quad_error=float(np.max(abs(a-average(k,192))))
    assert quad_error<1e-10
    assert abs(rows[-1]['native_triple_source_content_difference']-.54*expected**3)<1e-7
    return dict(round=849,status='working_not_formal',all_probe_checks_passed=True,
        original_full_Nambu_dimension=64,original_full_mass_matrix_retained=True,
        pulse_is_second_derivative_of_compact_bump=True,pulse_zero_moments=moments,
        compact_first_antiderivative_integral=antiderivative_integral,
        time_average_native_transverse_coefficient=coefficient,
        integration_by_parts_prediction=expected,
        first_two_block_mode_basis_rotates_Y_to_Z=True,
        principal_triple_source_content_difference=.54*expected**3,
        independent_quadrature_error=quad_error,rows=rows,
        pure_metric_source_test_on_original_background_constructed=False,
        original_gauge_and_matter_constraints_solved=False,
        original_inhomogeneous_propagation_bound_proven=False,
        fresh_formal_test_groups_counted=0)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');a=p.parse_args();r=run()
    if a.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
