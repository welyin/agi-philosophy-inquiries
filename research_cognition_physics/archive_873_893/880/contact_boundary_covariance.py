"""880: boundary pinching of the original730/805 full past CAR covariance.
Odd periodic Fourier collocation regulator; all64 Nambu components retained.
Not a simulation of the actual future curved PDE; the continuum obstruction
is established analytically from its original Dirac principal symbol.
"""
from pathlib import Path
import argparse,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from research_layout import Layout,ResearchRuntime
with ResearchRuntime(Layout()).installed():
    sys.path.insert(0,str(HERE.parent/'805'))
    import original_past_covariance_action as original
TARGET=HERE/'contact_boundary_covariance_results.json'
def fro2(a):return float(np.vdot(a,a).real)
def groups():
    ev,u=np.linalg.eigh(original.MASS@original.MASS)
    sets=[]
    for i,e in enumerate(ev):
        if not sets or abs(e-sets[-1][0])>1e-9:sets.append([float(e),[i]])
        else:sets[-1][1].append(i)
    return [(e,u[:,idx]@u[:,idx].conj().T) for e,idx in sets]
def calculate(n):
    # Wave numbers on a torus of side2*pi; odd N retains k<->-k exactly.
    freq=np.fft.fftfreq(n,d=1/n)
    ks=np.stack(np.meshgrid(freq,freq,freq,indexing='ij'),axis=-1)
    radius=np.sum(ks**2,axis=-1)
    matrices=[];coefs=[]
    for ev,proj in groups():
        inv=1/np.sqrt(radius+ev)
        matrices.append(original.MASS@proj)
        coefs.append(.5*inv)
        for a in range(3):
            matrices.append(original.GAMMA[a]@proj)
            coefs.append(.5*ks[...,a]*inv)
    matrices=np.array(matrices);coefs=np.array(coefs)
    kernels=np.fft.ifftn(coefs,axes=(1,2,3)).reshape(len(matrices),-1)
    flat=matrices.reshape(len(matrices),-1)
    gram=flat.conj()@flat.T
    # Only off-site values enter: the .5 I delta does not contribute.
    intensity=np.einsum('ax,ab,bx->x',kernels.conj(),gram,kernels,optimize=True).real.reshape(n,n,n)
    assert intensity.min()>-1e-12
    a=np.arange(n)<n//2;b=~a
    multiplicity=np.array([sum(bool(a[x] and b[(x-dx)%n]) for x in range(n))*n*n for dx in range(n)])
    cross_hs2=float(np.sum(multiplicity[:,None,None]*intensity))
    # Independent first-momentum matrix reconstruction and all discrete reality.
    sample=np.array([1,0,-1.],float)
    explicit=original.I-original.occupied(sample)
    reconstructed=.5*original.I
    for ev,proj in groups():
        reconstructed+=.5*(original.MASS+sum(original.GAMMA[j]*sample[j] for j in range(3)))@proj/np.sqrt(np.dot(sample,sample)+ev)
    error=float(np.max(abs(explicit-reconstructed)))
    # Cross-interface energy of pinching. Original H uses same Fourier modes.
    # h_r=sum_a Gamma_a D_a(r), r only along an axis. Only normal-axis
    # nonlocal hopping crosses the slab; onsite mass is unchanged.
    dx_kernel=np.fft.ifft(freq)
    energy=0.
    for dx in range(1,n):
        idx=dx*n*n
        pr=sum(kernels[j,idx]*matrices[j] for j in range(len(matrices)))
        hr=dx_kernel[dx]*original.GAMMA[0]
        energy += multiplicity[dx]*float(np.trace(hr@pr.conj().T).real)
    # -1/2 Tr H(P_pin-P); the two orientations cancel the Nambu1/2.
    assert energy>0 and cross_hs2>0 and error<2e-12
    return dict(points_per_axis=n,spatial_sites=n**3,Nambu_dimension=64*n**3,
        cross_block_HS_squared=cross_hs2,full_covariance_change_HS_squared=2*cross_hs2,
        HS_squared_per_N_squared=cross_hs2/n**2,
        pinching_energy_in_same_regulated_past_H=energy,energy_per_N_cubed=energy/n**3,
        mass_reconstruction_error=error)
def small_direct_check():
    n=3;dim=64;freq=np.fft.fftfreq(n,d=1/n)
    ps=[original.I-original.occupied(np.array([k,0,0.])) for k in freq]
    hs=[original.hamiltonian(np.array([k,0,0.])) for k in freq]
    pk=np.fft.ifft(np.array(ps),axis=0);hk=np.fft.ifft(np.array(hs),axis=0)
    p=np.block([[pk[(i-j)%n] for j in range(n)] for i in range(n)])
    h=np.block([[hk[(i-j)%n] for j in range(n)] for i in range(n)])
    pin=p.copy();pin[:dim,dim:]=0;pin[dim:,:dim]=0
    delta=pin-p
    charge=np.kron(np.eye(n),original.CHARGE)
    cross=fro2(p[:dim,dim:])
    direct=-.5*float(np.trace(h@delta).real)
    by_cross=float(np.trace(h[:dim,dim:]@p[dim:,:dim]).real)
    bounds=np.linalg.eigvalsh(pin)
    assert bounds.min()>-1e-12 and bounds.max()<1+1e-12
    assert abs(direct-by_cross)<1e-12 and abs(fro2(delta)-2*cross)<1e-12
    assert np.linalg.norm(charge@pin.conj()@charge+pin-np.eye(n*dim),2)<1e-12
    return dict(one_dimensional_regulator_sites=n,covariance_dimension=n*dim,
        product_covariance_min=float(bounds.min()),product_covariance_max=float(bounds.max()),
        Nambu_energy_formula_residual=abs(direct-by_cross),original_past_pinching_energy=direct,
        offdiagonal_HS_identity_residual=abs(fro2(delta)-2*cross))
def collar_integral_checks():
    # Leading flat-interface normal-symbol integral after tangential integration.
    # t=delta/z converts the remaining singular integral to a quadratic.
    points,weights=np.polynomial.legendre.leggauss(8)
    out=[];rank=64;upper=1.
    for delta in (.25,.125,.0625,.03125):
        lo=delta/upper
        t=(1+lo)/2+(1-lo)*points/2
        numerical=rank/(8*np.pi**3)*(1-lo)/2*np.sum(weights*(t-t*t))/delta**2
        exact=rank/(8*np.pi**3)*(1/(6*delta**2)-1/(2*upper**2)+delta/(3*upper**3))
        assert abs(numerical-exact)<1e-12
        out.append(dict(collar_width=delta,leading_symbol_integral_per_unit_area=float(exact),
            scaled_by_delta_squared=float(exact*delta**2),quadrature_residual=float(abs(numerical-exact))))
    return out
def run():
    rows=[calculate(n) for n in (5,9,13,17)]
    assert all(rows[i+1]['cross_block_HS_squared']>rows[i]['cross_block_HS_squared'] for i in range(len(rows)-1))
    assert all(rows[i+1]['pinching_energy_in_same_regulated_past_H']>rows[i]['pinching_energy_in_same_regulated_past_H'] for i in range(len(rows)-1))
    return dict(round=880,date='2026-10-06',fresh_numbered_groups=1,cumulative_numbered_groups=3665,
        all_checks_passed=True,original_mass_squared_blocks=len(groups()),
        original_Nambu_components=64,original_positive_modes_per_momentum=32,
        regulator_rows=rows,independent_position_matrix_check=small_direct_check(),
        continuum_leading_symbol_collar_checks=collar_integral_checks(),
        argument_scope='Across a smooth touching Cauchy interface, the original Dirac principal symbol makes the cross covariance non-Hilbert-Schmidt. The graded quasifree product of the same two reference marginals is not normal in the original global Fock folium. This blocks zero-collar reference factorization, not correlated regional composition, finite-mode preparation, or the entire Q/E program.',
        numerical_energy_is_original_auxiliary_past_not_actual_curved_future=True,
        continuum_obstruction_is_analytic_not_grid_extrapolation=True,
        pure_area_power_not_inferred_from_hard_Fourier_regulator=True,
        all_normal_nonGaussian_products_excluded=False,full_goal_completed=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    result=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
