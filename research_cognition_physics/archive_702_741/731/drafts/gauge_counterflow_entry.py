"""731 entry: a gauge-invariant momentum shift missing from the old radial menu."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
HERE=Path(__file__).resolve().parent;ARCHIVE=HERE.parent
sys.path.insert(0,str(ARCHIVE))
import joint_gravity_material_coordinates as old
TARGET=HERE/'gauge_counterflow_entry_results.json'

def run():
    N=16;f=old.fields(N);q=f['q'];a=q['f']['a'];E=q['f']['E'];dx=q['dx']
    x,y,z=np.moveaxis(q['grid'],-1,0)
    derivative=old.geo.derivative
    def D(v,i):return derivative(v,i)-np.cross(a[...,i,:],v)
    F=derivative(a[...,2,:],0)-derivative(a[...,0,:],2)-np.cross(a[...,0,:],a[...,2,:])
    w=np.sin(2*z);weighted=w[...,None]*F
    dE=np.zeros_like(E);dE[...,0,:]=D(weighted,2);dE[...,2,:]=-D(weighted,0)
    gauss=sum(D(dE[...,i,:],i) for i in range(3))
    momentum=np.zeros(3);covariant=np.zeros(3)
    for j in range(3):
        for i in range(3):
            da=derivative(a[...,i,:],j)
            momentum[j]+=dx**3*float(np.sum(dE[...,i,:]*da))
            fij=da-derivative(a[...,j,:],i)-np.cross(a[...,j,:],a[...,i,:])
            covariant[j]+=dx**3*float(np.sum(dE[...,i,:]*fij))
    exact=-2*np.pi**3*.12**2*.41**2*(1+.38**2/2)
    assert np.max(abs(gauss))<2e-13
    assert max(abs(momentum[:2]))<1e-13 and abs(momentum[2]-exact)<2e-13
    assert np.max(abs(momentum-covariant))<2e-13
    # Old radial Gram has only x/y range; the new direction closes the missing z.
    G=q['gram'];columns=np.column_stack((G[:,:2],momentum))
    determinant=float(np.linalg.det(columns))
    assert np.linalg.matrix_rank(G)==2 and abs(determinant)>1e-6
    target=np.array([.13,-.07,.11]) # declared source-space test, not a computed quantum source
    coeff=np.linalg.solve(columns,-target);residual=float(np.linalg.norm(columns@coeff+target))
    bw=old.geo.old.PAR['b'][1];weight=f['psi']**-2
    linear=2*bw*dx**3*float(np.sum(weight*np.sum(E*dE,axis=(-1,-2))))
    quadratic=bw*dx**3*float(np.sum(weight*np.sum(dE*dE,axis=(-1,-2))))
    assert quadratic>0 and residual<1e-13
    deps=('research_note_572.md','research_note_574.md','research_note_634.md',
          'research_note_730.md','joint_gravity_material_coordinates.py')
    return dict(entry_round=731,new_formal_round=False,N=N,
        covariant_Gauss_residual=float(np.max(abs(gauss))),
        canonical_momentum_shift=momentum.tolist(),covariant_momentum_shift=covariant.tolist(),
        analytic_z_shift=float(exact),original_radial_rank=int(np.linalg.matrix_rank(G)),
        augmented_determinant=determinant,declared_test_target=target.tolist(),
        coefficients=coeff.tolist(),total_momentum_cancellation_error=residual,
        old_geometry_electric_energy_linear=linear,old_geometry_electric_energy_quadratic=quadratic,
        no_new_configuration_or_scalar_charge=True,
        quantum_source_and_new_Hamiltonian_constraint_not_yet_solved=True,
        dependencies={n:hashlib.sha256((ARCHIVE/n).read_bytes()).hexdigest() for n in deps},
        scope='Original weak curvature supports an exact smooth gauge-invariant canonical momentum shift spanning the old missing z direction. Fixed configuration Gauss and global momentum checked. Source-space test is declared, not a transplanted634 flat-source measurement. Full state-dependent quantum Gauss, constraint backreaction and relative source definition remain open.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true');a=p.parse_args();r=run()
    if a.write_results:
        with TARGET.open('x',encoding='utf8') as f:json.dump(r,f,ensure_ascii=False,indent=2)
    else:assert r==json.loads(TARGET.read_text('utf8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))

