"""847: original 64-component coefficients and the original ten-CAR code.

This checks a complete smooth time average in the original constant auxiliary
symbol. The inhomogeneous spacetime wave-packet existence argument is analytic.
It does not compute a single-point nonlinear gravitational response.
"""
from pathlib import Path
import argparse,itertools,json,sys
import numpy as np
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent;TARGET=HERE/'native_triple_source_probe_results.json'
sys.path.insert(0,str(HERE.parent/'805'))
import original_past_covariance_action as original
sys.path.insert(0,str(HERE.parent/'829'))
import majorana_code_source_bridge as code

def run():
    gamma,_,_,_,comp,sx,sz=code.code_data()
    sy=[code.mul((0,0,3),code.mul(sz[b],sx[b])) for b in range(5)]
    logical=code.mul((0,0,1),code.product(gamma[i] for i in (0,1,4,5,12,14)))
    detected=[]
    for axes in itertools.product(range(3),repeat=3):
        op=code.product((sx,sy,sz)[axis][block] for axis,block in zip(axes,(0,1,3)))
        kind,phase=comp(code.mul(logical,op))
        if kind=='scalar':detected.append((axes,phase))
    assert detected==[((2,2,1),2)],detected
    paulis=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1.])]
    assert max(np.max(abs(g[30:32,30:32]-p)) for g,p in zip(original.GAMMA,paulis))<1e-14
    # Exact dual-rail one-body dictionary, including the Y and Z signs.
    ann=np.array([[0,1],[0,0]],complex);z=paulis[2]
    cs=[np.kron(ann,np.eye(2)),np.kron(z,ann)];odd=[1,2]
    dictionary_error=0.
    for j,p in enumerate(paulis):
        dg=sum(p[a,b]*cs[a].conj().T@cs[b] for a in range(2) for b in range(2))
        expected=paulis[j]*(1 if j==0 else -1)
        dictionary_error=max(dictionary_error,float(np.max(abs(dg[np.ix_(odd,odd)]-expected))))
    assert dictionary_error<1e-14
    width=.2
    def quadrature(n):
        x,w=np.polynomial.legendre.leggauss(n);w=w*np.exp(-1/(1-x*x));w/=w.sum()
        return width*x,w
    def avg(radius,direction,axis,n=96):
        h=original.hamiltonian([0,0,direction*radius]);e,u=np.linalg.eigh(h)
        times,weights=quadrature(n)
        multiplier=np.einsum('t,tij->ij',weights,np.cos(times[:,None,None]*(e[:,None]-e[None,:])[None,:,:]/radius))
        # K/radius from -delta H: sign reverses for the -momentum partner.
        return u@((u.conj().T@(direction*original.GAMMA[axis])@u)*multiplier)@u.conj().T
    times,weights=quadrature(96);c=float(np.dot(weights,np.cos(2*times)))
    assert c>=np.cos(2*width)>0
    delta=.9*.6
    rows=[];reality=0.;quadrature_error=0.
    for radius in (2.,8.,32.,128.):
        plus=[avg(radius,1,j) for j in range(3)];minus=[avg(radius,-1,j) for j in range(3)]
        for ap,am in zip(plus,minus):
            reality=max(reality,float(np.max(abs(original.CHARGE@ap.conj()@original.CHARGE+am))))
        restricted=[a[30:32,30:32] for a in plus]
        coeffs=[np.array([np.trace(a@p).real/2 for p in paulis]) for a in restricted]
        selected=[coeffs[2],coeffs[2],coeffs[1]];blocks=[0,1,3]
        pieces=[]
        for b,cc in zip(blocks,selected):pieces.append([(sx[b],cc[0]),(sy[b],-cc[1]),(sz[b],-cc[2])])
        response=0.
        for terms in itertools.product(*pieces):
            op=code.product(p for p,_ in terms);kind,phase=comp(code.mul(logical,op))
            if kind=='scalar':
                assert phase in (0,2)
                response+=np.prod([value for _,value in terms])*(1 if phase==0 else -1)
        assert response>.8*c
        # The original logical axis on these three blocks detects only Z,Z,Y.
        assert abs(response-selected[0][2]*selected[1][2]*selected[2][1])<1e-13
        target=[c*paulis[0],c*paulis[1],paulis[2]]
        error=max(float(np.max(abs(a-b))) for a,b in zip(restricted,target))
        if radius==8.:
            quadrature_error=max(float(np.max(abs(avg(radius,1,j,64)-plus[j]))) for j in range(3))
        rows.append(dict(momentum=radius,selected_2_by_2_Pauli_coefficients=[a.tolist() for a in selected],
                         original_mass_time_average_error_from_principal=error,
                         normalized_triple_content_difference=float(delta*response)))
    assert reality<1e-11 and quadrature_error<1e-10
    assert rows[-1]['original_mass_time_average_error_from_principal']<rows[0]['original_mass_time_average_error_from_principal']/1000
    assert abs(rows[-1]['normalized_triple_content_difference']-delta*c)<1e-5
    return dict(round=847,all_checks_passed=True,fresh_test_groups=1,
        original_Nambu_dimension=64,momentum_partners_both_retained=True,
        original_all_mass_and_mixing_matrices_retained=True,
        original_code_CAR_modes=10,chosen_code_blocks=blocks,
        half_width_in_rescaled_time=width,smooth_bump_average_cosine=c,
        analytic_lower_bound_cosine=float(np.cos(2*width)),
        dual_rail_one_body_dictionary_residual=dictionary_error,
        full_Nambu_reality_residual=reality,independent_time_quadrature_error=quadrature_error,
        coefficient_rows=rows,
        frozen_principal_normalized_difference=delta*c,
        time_window_is_nonzero_and_smooth=True,
        numerical_model_is_original_inhomogeneous_spacetime=False,
        spatial_wave_packet_lifting_is_analytic=True,
        nonlinear_one_point_geometry_coefficient_computed=False,
        native_autonomous_source_measurement_proven=False,
        EFT_bandwidth_or_infinite_energy_limit_claimed=False)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();out=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert out==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(out,ensure_ascii=False,indent=2))
