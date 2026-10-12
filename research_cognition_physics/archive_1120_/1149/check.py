"""Round 1149: finite matrix elements and static cross response; no full parent simulation."""
from pathlib import Path
import argparse, hashlib, importlib.util, json, math
import numpy as np
from numpy.polynomial.legendre import leggauss

BASE=Path(__file__).resolve().parent
QUAL=BASE.parent/'_shared/after1148_local_density_check.py'
spec=importlib.util.spec_from_file_location('density_qualification_1149',QUAL)
q=importlib.util.module_from_spec(spec)
spec.loader.exec_module(q)

def y(r):
    return np.exp(-r)/(4*np.pi*r)

def yp(r):
    return -np.exp(-r)/(4*np.pi)*(1/r+1/r**2)

def radial_probability(n):
    z,w=leggauss(n)
    r=(z+1)/2
    w=w/2
    v=q.radial(n)
    f=q.j0(v['q']*r)
    g=v['b']*q.j1(v['q']*r)
    rho=v['N_squared']*(f*f-g*g)/v['scalar_integral']
    return r,w,rho,v

def cross_group(v):
    r,w,rho,_=radial_probability(128)
    z,wz=leggauss(128)
    qr=v['normalized_density_yukawa_moment']
    qg=v['old_g_yukawa_moment']
    ratio=(qr/qg)**2
    match=qg/qr
    rows=[]
    for d in [2.05,3.,5.,8.]:
        dist=np.sqrt(d*d+r[:,None]**2-2*d*r[:,None]*z)
        # Integrate the second actual density against the first exterior potential.
        direct=2*np.pi*qr*np.sum(w*r*r*rho*(y(dist)@wz))
        fact=qr*qr*float(y(d))
        old=qg*qg*float(y(d))
        err=abs(direct-fact)/fact
        matched=direct*match**2
        assert err<1e-11
        assert abs(matched/old-1)<1e-11
        rows.append(dict(distance=d,direct_density_cross=float(direct),
                         factorized_cross=fact,old_cross=old,
                         matched_cross=float(matched),relative_residual=float(err)))
    delta=-np.pi*(ratio-1)
    # Fixed code, explicit phase-aligned Y-quadrature diagnostic.
    B=np.array([1,-np.exp(-1j*delta)],complex)/np.sqrt(2)
    B0=np.array([1,-1],complex)/np.sqrt(2)
    sigma_y=np.array([[0,-1j],[1j,0]],complex)
    effect=(np.eye(2)+sigma_y)/2
    p=float(np.vdot(B,effect@B).real)
    p0=float(np.vdot(B0,effect@B0).real)
    gap=abs(p-p0)
    assert abs(gap-abs(np.sin(delta))/2)<1e-14
    return dict(kind='all_separations_cross_matching_and_long_phase_diagnostic',passed=True,
                rows=rows,naive_relative_coefficient_difference=float(ratio-1),
                matching_a_over_kappa=float(match),
                conditional_energy_phase_difference=float(delta),
                aligned_code_operator_distance=float(2*abs(np.sin(delta/2))),
                same_Y_effect_probability_difference=gap,
                numerical_intervals_certified=False,full_parent_remainder_proved=False)

def force_group(v):
    amplitude=v['old_g_yukawa_moment'] # divide fields by lambda*j1
    rows=[]
    for d in [2.2,3.,5.]:
        for radius in [1.02,(d/2+1)/2]:
            assert 1<radius<d-1
            for n in [64,128]:
                z,w=leggauss(n)
                normal=np.column_stack([np.sqrt(1-z*z),np.zeros(n),z])
                pos=radius*normal
                off=pos-np.array([0,0,d])
                dist=np.linalg.norm(off,axis=1)
                phiA=-amplitude*float(y(radius))
                phiB=-amplitude*y(dist)
                gradA=-amplitude*float(yp(radius))*normal
                gradB=-amplitude*yp(dist)[:,None]*off/dist[:,None]
                dotBn=np.sum(gradB*normal,axis=1)
                dotAn=np.sum(gradA*normal,axis=1)
                dotAB=np.sum(gradA*gradB,axis=1)
                Tn=gradA*dotBn[:,None]+gradB*dotAn[:,None]-normal*(dotAB+phiA*phiB)[:,None]
                force_z=-2*np.pi*radius**2*np.sum(w*Tn[:,2])
                expected=-amplitude**2*float(yp(d))
                relative=abs(force_z-expected)/expected
                assert relative<2e-11
                rows.append(dict(distance=d,sphere_radius=radius,order=n,
                                 normalized_flux_force=float(force_z),
                                 normalized_potential_force=expected,
                                 relative_residual=float(relative)))
    return dict(kind='same_static_cross_field_surface_flux_and_force',passed=True,
                normalization='divide by lambda^2*j1^2',rows=rows,
                full_quantum_stress_claimed=False,autonomous_wall_claimed=False)

def compute():
    manifest=json.loads((BASE/'source_hashes.json').read_text(encoding='utf-8'))
    for path,sha in manifest.items():
        assert hashlib.sha256((BASE/path).read_bytes()).hexdigest()==sha, path
    radial=[q.radial(n) for n in [64,128,256,512]]
    spread={key:max(v[key] for v in radial)-min(v[key] for v in radial)
            for key in ['scalar_integral','normalized_density_yukawa_moment','old_g_yukawa_moment']}
    assert max(spread.values())<1e-12
    material=dict(kind='actual_ground_doublet_and_radial_density',passed=True,
                  reused_qualification_not_previously_counted=True,
                  spin=q.spin_group(),radial=radial,spread=spread)
    groups=[material,cross_group(radial[1]),force_group(radial[1])]
    return dict(round=1149,all_passed=True,numerical_check_groups=len(groups),
                groups=groups,frozen_source_count=len(manifest),
                scientific_calibration_increment=0,new_adopted_axioms=0,
                full_parent_construction=False,uniform_long_time_remainder=False)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-results',action='store_true')
    args=parser.parse_args()
    result=compute()
    path=BASE/'results.json'
    output=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.write_results:
        if path.exists():
            assert path.read_text(encoding='utf-8')==output,'refuse different result overwrite'
        else: path.write_text(output,encoding='utf-8')
    else:
        assert json.loads(path.read_text(encoding='utf-8'))==result
    print(json.dumps(dict(status='PASS',groups=3,
                          sources=result['frozen_source_count'],
                          static_flux_max_relative_residual=max(x['relative_residual'] for x in result['groups'][2]['rows']),
                          coefficient_mismatch=result['groups'][1]['naive_relative_coefficient_difference'],
                          goal_complete=False)))
if __name__=='__main__':
    main()
