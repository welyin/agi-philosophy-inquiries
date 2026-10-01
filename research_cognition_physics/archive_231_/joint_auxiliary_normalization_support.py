"""658: exact zero-sector support of the original full sphere average.

Analytic null modes, not a spectral threshold defining a new physical model.
Finite-box strict normalization does not imply a uniform continuum limit.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import joint_auxiliary_reflection_gluing as old

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_auxiliary_normalization_support_results.json'


def zero_frame(nx,length):
    ev,v=np.linalg.eigh(old.old.old.spin.GAMMA[3]);wp=v[:,ev>0]
    wm=old.old.old.spin.G5@wp
    z=np.zeros((4*nx*length,2*(length-1)),complex)
    for t in range(length-1):
        for x in range(nx):
            z[4*(t*nx+x):4*(t*nx+x+1),2*t:2*t+2]=wp/np.sqrt(2*nx)
            z[4*((t+1)*nx+x):4*((t+1)*nx+x+1),2*t:2*t+2]=-wm/np.sqrt(2*nx)
    return z


def support_check():
    rows=[]
    for nx,length in ((1,2),(2,2),(3,2),(3,3),(4,3)):
        f=old.old.frame(1,nx,2*length);d=old.reflection(f)['D'];z=zero_frame(nx,length)
        residual=float(np.max(abs(d@z)));orth=float(np.max(abs(z.conj().T@z-np.eye(z.shape[1]))))
        assert residual<3e-14 and orth<2e-14
        lifted=d+z@z.conj().T;eig=np.linalg.eigvalsh(lifted)
        assert min(eig)>0
        logdet=float(np.linalg.slogdet(lifted)[1])
        if nx==1:assert abs(logdet+4*np.log(2))<3e-13
        rows.append(dict(spatial_sites=nx,time_sites=2*length,analytic_zero_spin_modes=z.shape[1],
            null_residual=residual,orthogonality_residual=orth,lifted_smallest_eigenvalue=float(min(eig)),
            spin_pseudodeterminant_log=logdet,projection_cutoff_used=False))
    return dict(rows=rows,nullity_formula='2*(half_time_length-1); only spatial momentum zero',
                arbitrary_other_momentum_strictness_is_analytic_spectral_integral=True)


def restricted_pairing_check():
    rng=np.random.default_rng(65801);rows=[]
    for nx,length in ((1,2),(3,2),(3,3),(4,3)):
        z=zero_frame(nx,length);zz=np.kron(z,np.eye(16));baseline=np.zeros((nx*length,10));baseline[:,0]=1
        a0=zz.T@old.field_pair(baseline)@zz;p0=old.old.old.pfaffian(a0)
        assert abs(abs(p0)-1)<1e-12
        for width in (.12,.5,1.4):
            e=np.eye(10)[0]+width*rng.normal(size=(length,nx,10));e/=np.linalg.norm(e,axis=2)[:,:,None]
            actual=old.old.old.pfaffian(zz.T@old.field_pair(e.reshape(-1,10))@zz)/p0
            means=e.mean(axis=1)
            predicted=float(np.prod([np.linalg.norm((means[t]+means[t+1])/2)**16 for t in range(length-1)]))
            error=float(abs(actual-predicted));assert error<3e-13
            assert actual.real>=-1e-13
            rows.append(dict(spatial_sites=nx,half_time_length=length,width=width,
                actual_zero_sector_Pfaffian=[float(actual.real),float(actual.imag)],
                positive_open_chain_product=predicted,error=error))
    return dict(rows=rows,spatial_means_are_not_assumed_unit_vectors=True,
                signed_Pfaffian_not_replaced_by_absolute_value=True)


def rising(a,n):
    out=Fraction(1)
    for k in range(n):out*=a+k
    return out


def moment_and_bound_check():
    lam=rising(Fraction(9,2),8)/rising(Fraction(9),8)
    spectrum=[]
    for ell in range(9):
        value=lam*Fraction(math.factorial(8),math.factorial(8-ell))/rising(Fraction(17),ell)
        dim=math.comb(ell+9,9)-(math.comb(ell+7,9) if ell>=2 else 0)
        spectrum.append((dim,value))
    rows=[]
    for length in (1,2,3,4):
        amplitude=lam**(length-1)
        lower=Fraction(1,2**32)*amplitude**2
        exact_partition=sum(dim*value**(2*length) for dim,value in spectrum)
        assert exact_partition>=lower>0
        rows.append(dict(spatial_sites=1,time_sites=2*length,
            exact_zero_sector_mean=str(amplitude),exact_partition=str(exact_partition),
            support_lower_bound=str(lower),lower_bound_ratio=float(lower/exact_partition)))
    cap_lower=32/(35*np.pi)*(7/16)**3.5
    assert cap_lower>0
    # Uses original nontrivial spatial matrix; no integrations discarded.
    other=[]
    for length in (2,3):
        nx=3;d=old.reflection(old.old.frame(1,nx,2*length))['D'];z=zero_frame(nx,length)
        logdet=np.linalg.slogdet(d+z@z.conj().T)[1]
        log_lower=8*logdet+2*nx*length*np.log(cap_lower)-32*(length-1)*np.log(2)
        other.append(dict(spatial_sites=nx,time_sites=2*length,
            explicit_full_S9_log_partition_lower_bound=float(log_lower)))
    return dict(original_time_chain_lambda0=str(lam),open_chain_rows=rows,
        S9_cap_probability_lower_bound=float(cap_lower),nontrivial_spatial_bounds=other,
        cap_is_proof_subset_not_physical_cutoff=True,
        finite_box_lower_bound_not_asserted_uniform_in_volume_or_time=True)


def run():
    deps=('joint_auxiliary_reflection_gluing.py','joint_spatial_auxiliary_geometry.py',
          'joint_subgroup_measure_source.py','research_note_653.md','research_note_657.md')
    return dict(date='2026-10-02',round=658,tests_run=3,failures=0,errors=0,
        null_support=support_check(),actual_restricted_pairing=restricted_pairing_check(),
        sphere_mean_and_bounds=moment_and_bound_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope='Strictly positive original auxiliary S9 normalization on every finite free periodic spatial box and even AP time, using analytic zero modes and original average-vector support. No common one-step semigroup, physical CAR identification, gauge interaction, uniform continuum bound or quantum GR.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');args=p.parse_args();r=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8') as stream:json.dump(r,stream,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps({k:r[k] for k in ('round','tests_run','failures','errors')}))
