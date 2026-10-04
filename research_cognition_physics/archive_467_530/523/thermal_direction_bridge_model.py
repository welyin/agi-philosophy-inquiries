"""Conditional thermal-reference to actual qubit-direction instrument bridge.

The continuum geometry, spherical UV cutoff, Gaussian source, classical record
processing and a complete qubit carrier are inputs. No unique FUCP dimension claim.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
TARGET=HERE/"thermal_direction_bridge_results.json"
I=np.eye(2,dtype=complex)
PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],dtype=complex)


def quad(order,lo,hi):
    x,w=np.polynomial.legendre.leggauss(order)
    return lo+(hi-lo)*(x+1)/2,w*(hi-lo)/2


def cq(ell,temp=1.):
    omega=np.sqrt(ell)
    return .5/(omega*np.tanh(omega/(2*temp))) if temp else .5/omega


def projected_mean(d,s,order=192):
    v,w=quad(order,0,1)
    return float(s*math.sqrt(2/math.pi)*np.sum(w*(1-v*v)**((d-1)/2)*np.exp(-s*s*v*v/2)))


def projected_three(s):
    return (1-1/s**2)*math.erf(s/math.sqrt(2))+math.sqrt(2/math.pi)*math.exp(-s*s/2)/s


def radial_polar(s,order=160):
    radius,wr=quad(order,0,s+12)
    angle,wa=quad(order,-1,1)
    dens=radius[:,None]**2*np.exp(-.5*(radius[:,None]-s)**2-radius[:,None]*s*(1-angle[None,:]))/math.sqrt(2*math.pi)
    weight=dens*wr[:,None]*wa[None,:]
    return float(np.sum(weight)),float(np.sum(weight*angle)),float(np.sum(weight*angle**2))


def one_minus_sinc(x):
    out=np.empty_like(x)
    small=abs(x)<1e-3
    z=x[small]**2
    out[small]=z/6-z*z/120+z*z*z/5040
    out[~small]=1-np.sin(x[~small])/x[~small]
    return out


def difference_variance(d,length,order=256,temp=1.,cutoff=2.,mu2=1.25,ratio=.5,nu=.25,b=1.):
    k,w=quad(order,0,cutoff)
    if d==3:
        angular=one_minus_sinc(k*length)
        weight=w*k*k/math.pi**2
    elif d==2:
        theta,wt=quad(order,0,math.pi)
        angular=2*np.sin(k[:,None]*length*np.cos(theta)[None,:]/2)**2 @ (wt/math.pi)
        weight=w*k/math.pi
    else:
        raise ValueError("Only the two numerical examples are integrated here")
    dr=float(weight@(cq(k*k,temp)*angular))
    dm=float(weight@(cq(k*k+mu2,temp)*angular))
    meter=2*(1+ratio*ratio)*nu*nu/(ratio*ratio)
    return dict(d=d,length=length,common=dr,massive=dm/ratio**2,meter=meter,
                coordinate_variance=(dr+dm/ratio**2+meter)/b**2)


def qubit_instrument(s,eta=.8):
    norm,a,mzz=radial_polar(s)
    assert abs(norm-1)<2e-11
    assert abs(a-projected_three(s))<2e-11
    moments=np.array([(1-mzz)/2,(1-mzz)/2,mzz])
    alpha=(math.sqrt((1+eta)/2)+math.sqrt((1-eta)/2))/2
    beta=(math.sqrt((1+eta)/2)-math.sqrt((1-eta)/2))/2
    vectors=[p.ravel(order="F") for p in PAULI]
    ident=I.ravel(order="F")
    common=alpha**2*np.outer(ident,ident.conj())
    common+=sum(beta**2*m*np.outer(v,v.conj()) for m,v in zip(moments,vectors))
    cross=alpha*beta*a*(np.outer(ident,vectors[2].conj())+np.outer(vectors[2],ident.conj()))
    jp,jm=common+cross,common-cross
    partial=lambda j:np.einsum('aiaj->ij',j.reshape((2,2,2,2),order="F"))
    expected=(I+eta*a*PAULI[2])/2
    assert np.max(abs(partial(jp).T-expected))<1e-11
    assert np.max(abs(partial(jp+jm)-I))<1e-11
    assert min(np.linalg.eigvalsh(jp)[0],np.linalg.eigvalsh(jm)[0])>-1e-12
    vals,vecs=np.linalg.eigh(expected)
    root=(vecs*np.sqrt(vals))@vecs.conj().T
    v=root.ravel(order="F")
    luders=np.outer(v,v.conj())
    gap=float(np.linalg.norm(jp-luders))
    assert gap>1e-5
    # Pointwise Kraus completeness, including an arbitrary complex orientation.
    direction=np.array([2.,-1.,3.]); direction/=np.linalg.norm(direction)
    spin=np.einsum('i,ijk->jk',direction,PAULI)
    kp,km=alpha*I+beta*spin,alpha*I-beta*spin
    assert np.max(abs(kp.conj().T@kp+km.conj().T@km-I))<1e-14
    return dict(snr=s,normalization_residual=abs(norm-1),mean_direction=a,
                second_axial_moment=mzz,choi_minimum=float(np.linalg.eigvalsh(jp)[0]),
                tp_residual=float(np.max(abs(partial(jp+jm)-I))),
                actual_vs_averaged_luders_choi_frobenius=gap)


def run():
    temp,cutoff,mu2,ratio,nu,b,eta=1.,2.,1.25,.5,.25,1.,.8
    # 1. Actual d=3 field covariance, two independent quadrature resolutions.
    rows=[]; quadrature_error=0.
    for length in (1.,8.,32.):
        low=difference_variance(3,length,128)
        high=difference_variance(3,length,256)
        quadrature_error=max(quadrature_error,abs(low["coordinate_variance"]-high["coordinate_variance"]))
        high["snr"]=length/math.sqrt(high["coordinate_variance"])
        high["ideal_antipodal_gap"]=eta*projected_three(high["snr"])
        rows.append(high)
    assert quadrature_error<2e-11
    common_bound=2*(temp*cutoff+cutoff**2/4)/math.pi**2
    massive_bound=2*cutoff**3*float(cq(mu2,temp))/(3*math.pi**2*ratio**2)
    meter_bound=2*(1+ratio**2)*nu**2/ratio**2
    sigma2_bound=(common_bound+massive_bound+meter_bound)/b**2
    assert all(row["coordinate_variance"]<=sigma2_bound for row in rows)
    # 2. Projected-normal formula independently from radial/angular integration.
    projection=[]
    for snr in (.1,1.,3.,8.):
        integral=projected_mean(3,snr)
        norm,polar,_=radial_polar(snr)
        assert abs(norm-1)<2e-11
        assert abs(integral-projected_three(snr))<2e-11 and abs(polar-integral)<2e-11
        projection.append(dict(snr=snr,laplace_integral=integral,closed_form=projected_three(snr),polar_integral=polar))
    # 3. A finite angular menu keeps a rigorous uniform antipodal gap.
    minimum_length=8.
    response_lower=1-6*sigma2_bound/minimum_length**2
    m,n=16,32
    cover=math.pi/(2*m)+math.pi/n
    lower=eta*(response_lower-cover)
    assert lower>0
    rng=np.random.default_rng(882)
    directions=rng.normal(size=(2048,3)); directions/=np.linalg.norm(directions,axis=1)[:,None]
    theta=np.arccos(directions[:,2]); phi=np.mod(np.arctan2(directions[:,1],directions[:,0]),2*math.pi)
    th=np.rint(theta*m/math.pi)*math.pi/m
    ph=np.rint(phi*n/(2*math.pi))*2*math.pi/n
    rounded=np.column_stack([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
    rounding_max=float(np.max(np.linalg.norm(directions-rounded,axis=1)))
    assert rounding_max<=cover
    # 4. Full post-measurement qubit channel is CP/TP, not averaged-effect Luders.
    instruments=[qubit_instrument(snr,eta) for snr in (1.,3.)]
    # 5. Dimension two defeats the implication direction contrast => absolute precision.
    lower_dim=[]
    for length in (8.,32.,128.):
        fine=difference_variance(2,length,512)
        coarse=difference_variance(2,length,384)
        assert abs(fine["coordinate_variance"]-coarse["coordinate_variance"])<1e-9
        fine["snr"]=length/math.sqrt(fine["coordinate_variance"])
        fine["direction_gap"]=eta*projected_mean(2,fine["snr"],384)
        lower_dim.append(fine)
    a=(temp+cutoff)/math.pi+float(cq(mu2,temp))*cutoff**2/(math.pi*ratio**2)+meter_bound
    log_coeff=2*temp/math.pi
    two_ratio_bound=(a+log_coeff*math.log(cutoff*minimum_length/2))/(b**2*minimum_length**2)
    two_gap_lower=eta*(1-4*two_ratio_bound)
    assert two_gap_lower>0
    assert lower_dim[-1]["common"]>lower_dim[0]["common"]
    assert all(x["direction_gap"]>=two_gap_lower for x in lower_dim)
    # 6. A larger carrier admits d=4 at fixed trace; it is not a complete qubit.
    gammas=np.array([np.kron(PAULI[0],I),np.kron(PAULI[1],I),
                     np.kron(PAULI[2],PAULI[0]),np.kron(PAULI[2],PAULI[1])])
    anticommutator_error=0.
    for i in range(4):
        for j in range(4):
            anticommutator_error=max(anticommutator_error,float(np.max(abs(gammas[i]@gammas[j]+gammas[j]@gammas[i]-(2*np.eye(4) if i==j else 0)))))
    assert anticommutator_error<1e-14
    four_gaps=[]; four_response=projected_mean(4,3.)
    for u in np.vstack([np.eye(4),np.ones((1,4))/2]):
        spin=np.einsum('i,ijk->jk',u,gammas)
        effect=(np.eye(4)+eta*four_response*spin)/2
        assert np.min(np.linalg.eigvalsh(effect))>=0 and np.max(np.linalg.eigvalsh(effect))<=1
        assert abs(np.trace(effect)-2)<1e-14
        four_gaps.append(float(np.linalg.norm(2*effect-np.eye(4),2)))
    assert max(abs(g-eta*four_response) for g in four_gaps)<1e-12
    cutoff_density=cutoff**3/(6*math.pi**2)
    return dict(date="2026-09-30",diagnostic_groups=6,failures=0,errors=0,
        numbered_round_created=False,numbered_test_increment=0,
        inputs=dict(temperature=temp,uv_cutoff=cutoff,massive_squared=mu2,ratio=ratio,
            readout_sigma=nu,reference_gradient=b,qubit_sharpness=eta,
            effective_continuum_geometry_is_input=True),
        diagnostics=dict(three_dimensional_readout=rows,field_quadrature_difference=quadrature_error,
            projected_normal=projection,
            uniform_finite_menu=dict(minimum_length=minimum_length,coordinate_variance_upper=sigma2_bound,
                response_lower=response_lower,cover_radius_upper=cover,antipodal_gap_lower=lower,
                angular_symbols=(m-1)*n+2,angular_bits=math.ceil(math.log2((m-1)*n+2)),
                random_rounding_crosscheck_max=rounding_max),
            qubit_instruments=instruments,
            lower_dimensional_counterexample=dict(examples=lower_dim,
                uniform_direction_gap_lower=two_gap_lower,absolute_difference_variance_is_unbounded=True),
            larger_carrier=dict(dimension=4,higher_dimensional_direction=4,
                clifford_residual=anticommutator_error,antipodal_gaps=four_gaps),
            resources=dict(cutoff_mode_density=cutoff_density,
                six_field_read_energy_increment=3*cutoff_density/(4*nu**2),
                reference_mean_gradient_energy_density=3*b*b/2,
                autonomous_qubit_control_energy_derived=False)),
        provenance_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
             for name in ("thermal_reference_dimension_model.py","thermal_reference_dimension_results.json")},
        scope=dict(conditional_three_dimensional_joint_model=True,fu_cp_unique_derivation=False,
            same_lattice_continuum_limit_proved=False,absolute_precision_separate_requirement=True,
            all_directions_continuum_is_input=True,finite_quantized_menu=True,
            exact_finite_record_budget_includes_all_hardware=False,
            field_source_independent_of_unknown_qubit=True,infinite_joint_confidence=False,
            finite_total_infinite_volume_energy=False,einstein_backreaction_solved=False))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-results",action="store_true")
    parser.add_argument("--check",action="store_true")
    args=parser.parse_args(); result=run()
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as f:
            f.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    if args.check:
        assert json.loads(TARGET.read_text("utf8"))==json.loads(json.dumps(result))
    print(json.dumps(result,ensure_ascii=False,indent=2))
