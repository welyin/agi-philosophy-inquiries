"""Change only the terminal heavy Majorana source at a common threshold; CP-sign inverse."""
import argparse
import json
import math
import runpy
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


def previous(name):
    return runpy.run_path(str(HERE/name))


def charged_at_threshold(fr, rep, point):
    """Transport the matched charged frame to a common scale just above the heavy cluster."""
    basis=np.array(rep["charged_left_MZ"]["real"])+1j*np.array(rep["charged_left_MZ"]["imag"])
    y=np.zeros((4,3,4),complex)
    y[0,:,:3]=np.diag(fr["central"]("up"))
    y[1,:,:3]=fr["low_ckm"](rep["gamma_MZ_deg"],fr["angles"]("direct"))@np.diag(fr["central"]("down"))
    y[2,:,:3]=basis*fr["central"]("charged_leptons")
    for a,susy,dt,grid in fr["segments"](.12):
        if a>=point:
            break
        assert not susy
        for j,m in enumerate(fr["fields"]["neutral_poles_GeV"]):
            if a==m:
                y[3,:,j]=fr["physics"]["Y0"][:,j]
        for j in range(0,len(grid)-1,2):
            f1=fr["beta"](y,grid[j],False)
            f2=fr["beta"](y+dt*f1/2,grid[j+1],False)
            f3=fr["beta"](y+dt*f2/2,grid[j+1],False)
            f4=fr["beta"](y+dt*f3,grid[j+2],False)
            y+=dt*(f1+2*f2+2*f3+f4)/6
    return np.linalg.svd(y[2,:,:3])[0][:,::-1]


def modes(m):
    vals,vectors=np.linalg.eigh(m.conj().T@m)
    return np.sqrt(np.maximum(vals,0)),vectors


def J(v):
    return (v[0,0]*v[1,1]*(v[0,1]*v[1,0]).conjugate()).imag


def inverse(z, h):
    """Sharp least terminal off-diagonal, varying all unobservable Majorana phases."""
    m=z@h@z.T
    sv,u=modes(m)
    ubar=(u[:,0]/u[:,0].conj())[:,None]*u.conj()
    columns=np.linalg.pinv(z)@ubar[:,1:].conj()
    np.testing.assert_allclose(z@columns,ubar[:,1:].conj(),atol=1e-13)
    coeff=sv[1:]*columns[0]*columns[1]
    phase=-np.exp(1j*(np.angle(coeff[0])-np.angle(coeff[1])))
    hp=(columns*(sv[1:]*[1.,phase]))@columns.T
    v2,up=modes(z@hp@z.T)
    assert abs(J(u))>1e-5
    np.testing.assert_allclose(v2[1:],sv[1:],rtol=1e-12)
    np.testing.assert_allclose(abs(up)**2,abs(u)**2,atol=2e-14)
    np.testing.assert_allclose(J(up),-J(u),atol=2e-14)
    low=abs(abs(coeff[0])-abs(coeff[1]))
    np.testing.assert_allclose(abs(hp[0,1]),low,atol=1e-18)
    e_coeff=sv[1:]*ubar[0,1:].conj()**2
    cross=e_coeff[1]*e_coeff[0].conjugate()
    cosine=(abs(m[0,0])**2-np.dot(abs(e_coeff),abs(e_coeff)))/(2*abs(cross))
    assert -1-1e-12<=cosine<=1+1e-12
    span=math.acos(np.clip(cosine,-1,1))
    phase_options=[np.exp(1j*(-np.angle(cross)+sign*span)) for sign in (-1,1)]
    alternatives=[]
    for rot in phase_options:
        candidate=(columns*(sv[1:]*[1.,rot]))@columns.T
        new_mass=z@candidate@z.T
        np.testing.assert_allclose(abs(new_mass[0,0]),abs(m[0,0]),rtol=2e-13,atol=1e-25)
        _,vectors=modes(new_mass)
        np.testing.assert_allclose(J(vectors),-J(u),atol=2e-14)
        alternatives.append(candidate)
    best=min(alternatives,key=lambda mat:abs(mat[0,1]))
    return hp,J(u),low,sv,abs(u)**2,best


def full_observables(d,c,h):
    heavy=np.block([[0*h,c],[c.T,h]])
    weak=np.column_stack([d,0*d])@np.linalg.inv(heavy)
    return weak@weak.conj().T,np.linalg.svd(heavy,compute_uv=False)


def propagated_inverse():
    """Invert the transported symmetric terminal tensor on the angle-passing lower-yd point."""
    model=previous("check.py")
    point=read("results.json")["cases"]["lower_down_k1_shared_third"]["roots"][1]
    assert all(point["windows_pass"])
    d,rep,track,*_=model["graph"](1,-1,"shared_third",.20)
    e=model["zer"]["unpack"](rep["charged_left_MZ"])
    pole=model["m"]
    v=model["V"]
    basis=np.array([[[1,0],[0,0]],[[0,1],[1,0]],[[0,0],[0,1]]])*1e-6
    hs=v.T@basis@v
    divided=pole[0]**2/(pole[0]**2-pole[1]**2)
    z=d/pole
    split=np.array([[0,divided],[divided,1]])
    tree=z@hs@z.T*1e9
    upper=z@(split*hs)@z.T*1e9
    # Independent resolvent on the mixed 1-keV insertion verifies the off-diagonal projector.
    mass=np.diag(pole)
    inv=np.diag(1/pole)
    radius=(pole[1]**2-pole[0]**2)/3
    mixed=hs[1]
    off=mass@mixed
    contour=np.zeros((3,3),complex)
    for angle in 2*np.pi*(np.arange(48)+.5)/48:
        step=radius*np.exp(1j*angle)
        top=(pole[1]**2+step)*np.eye(2)-mass@mass
        bottom=top-mixed.conj().T@mixed
        aux=np.linalg.solve(top-off@np.linalg.solve(bottom,off.conj().T),
                           -inv@mixed@inv+off@np.linalg.solve(bottom,inv))
        contour-=d@aux@d.T*1e9*step/48
    np.testing.assert_allclose(contour,upper[1],rtol=1e-12,atol=1e-17)
    response=e.T@model["lower_tensors"](np.array([tree,upper]),track)@e
    cols=response.reshape(3,9).T
    multipliers=np.array(point["H0_H6_multipliers"])
    h0=np.diag(model["H0"]).real*multipliers/1e-6
    orig=(cols@np.array([h0[0],0,h0[1]])).reshape(3,3)
    mass_eV,u=modes(orig)
    np.testing.assert_allclose(J(u),point["J"],atol=2e-12)
    np.testing.assert_allclose(mass_eV[1:],point["masses_eV"][1:],rtol=1e-10)
    ubar=(u[:,0]/u[:,0].conj())[:,None]*u.conj()
    targets=np.stack([s*np.outer(col.conj(),col.conj()).ravel()
                      for s,col in zip(mass_eV[1:],ubar[:,1:].T)],axis=1)
    terms,_,rank,_=np.linalg.lstsq(cols,targets,rcond=None)
    residual=np.linalg.norm(cols@terms-targets)/np.linalg.norm(targets)
    assert rank==3 and residual<1e-9, residual
    a=terms[1]
    phase=-np.exp(1j*(np.angle(a[0])-np.angle(a[1])))
    electron=targets[0]
    cross=electron[1]*electron[0].conjugate()
    cos_value=(abs(orig[0,0])**2-np.vdot(electron,electron).real)/(2*abs(cross))
    assert abs(cos_value)<=1+1e-10
    span=math.acos(np.clip(cos_value,-1,1))
    record=[]
    for rot in [phase]+[np.exp(1j*(-np.angle(cross)+sign*span)) for sign in (-1,1)]:
        source=terms@np.array([1,rot])
        low=(cols@source).reshape(3,3)
        sv,vec=modes(low)
        np.testing.assert_allclose(sv[1:],mass_eV[1:],rtol=2e-8)
        np.testing.assert_allclose(abs(vec)**2,abs(u)**2,atol=2e-8)
        np.testing.assert_allclose(J(vec),-J(u),atol=2e-8)
        record.append(dict(offdiag_keV=abs(source[1]),mbb_meV=1000*abs(low[0,0]),
                           J=J(vec),terminal_real_keV=source.real.tolist(),
                           terminal_imag_keV=source.imag.tolist()))
    for item in record[1:]:
        np.testing.assert_allclose(item["mbb_meV"],1000*abs(orig[0,0]),rtol=1e-9)
    return dict(original_J=J(u),original_mbb_meV=1000*abs(orig[0,0]),
                original_diagonal_keV=h0.tolist(),minimum_mixed=record[0],
                minimum_also_matching_mbb=min(record[1:],key=lambda r:r["offdiag_keV"]),
                relative_residual_reflected_basis=residual)


def compute():
    common=read("../1698/results.json")["shared_witness"]
    prev=previous("../1698/check.py")
    d,c,h=prev["neutral_graph"](read("../1698/inputs.json"))
    b,t=common["neutral_b"],common["neutral_t"]
    d,c,h=math.sqrt(b)*d,b*math.sqrt(t)*c,b*t*h
    cfg=read("../1703/inputs.json")
    benchmark=read("../1703/results.json")["cases"]["central_k2_shared_third"]
    phi=cfg["phase"]
    pre=np.exp(1j*phi*np.array([2,0,1]))
    oriented=np.exp(1j*np.array([0,2*phi,3*phi+np.angle(3+2*np.exp(4j*phi))]))
    leg=oriented/pre*np.array([1.,1.,benchmark["third_neutral_modulus"]])
    d=leg[:,None]*d
    fr=previous("../1702/check.py")
    fr["physics"]["Y0"]=(leg[:,None]*fr["physics"]["Y0"]).conj()
    rep=fr["solve"](2,"direct",seed=benchmark["gamma_quark_MZ_deg"],step=.12)
    scale=max(common["neutral_poles_GeV"])
    e=charged_at_threshold(fr,rep,scale)
    z=np.linalg.solve(c,d.T).T
    charged_source=e.T@z
    hp,j,lower,sv,p,both=inverse(charged_source,h)
    w,poles=full_observables(d,c,h)
    _,poles2=full_observables(d,c,hp)
    # Remove common BB† analytically; subtracting two O(10^-5) matrices loses this effect.
    old=z@h@np.linalg.inv(c)
    new=z@hp@np.linalg.inv(c)
    difference=new@new.conj().T-old@old.conj().T
    relative=np.linalg.norm(difference,2)/np.linalg.norm(w,2)
    dh=np.linalg.norm(hp-h,2)
    light_ee=lambda term:abs((charged_source@term@charged_source.T)[0,0])*1e12
    assert max(abs(poles2-poles))<=dh+1e-10
    assert lower>0 and relative<1e-12
    # The per-pole positive quadratic Higgs moment stays fixed when h changes.
    mH=np.linalg.norm(d@c.conj())**2
    port=np.column_stack([d,0*d])
    for end in (h,hp):
        heavy=np.block([[0*h,c],[c.T,end]])
        np.testing.assert_allclose(np.linalg.norm(port@heavy.conj().T)**2,mH,rtol=2e-15)
    ans=dict(common_threshold_GeV=scale,endpoint_gamma=rep["gamma_MZ_deg"],J=j,
             J_competitor=-j,matching_scale_mass_pair_eV=(sv[1:]*1e9).tolist(),
             mixing_moduli_squared=p.tolist(),sharp_minimum_abs_terminal_offdiag_GeV=lower,
             ratios_to_geometric_mean_original_diagonal=lower/math.sqrt(abs(h[0,0]*h[1,1])),
             terminal_matrix_original_GeV=h.real.tolist(),
             terminal_matrix_competitor_real_GeV=hp.real.tolist(),
             terminal_matrix_competitor_imag_GeV=hp.imag.tolist(),
             relative_change_of_weak_matrix_norm=relative,
             max_actual_heavy_pole_shift_GeV=max(abs(poles2-poles)),
             rigorous_heavy_pole_shift_bound_GeV=dh,
             relative_heavy_pole_shift_bound=dh/min(poles),
             positive_Higgs_threshold_moment_ratio=1.)
    ans["light_exchange_mbb_original_meV"]=light_ee(h)
    ans["light_exchange_mbb_at_minimum_offdiag_meV"]=light_ee(hp)
    ans["minimum_offdiag_also_preserving_mbb_GeV"]=abs(both[0,1])
    ans["light_exchange_mbb_after_additional_matching_meV"]=light_ee(both)
    ans["heavy_shift_bound_also_preserving_mbb_GeV"]=np.linalg.norm(both-h,2)
    assert abs(both[0,1])>=lower-1e-18
    ans["running_lower_down_k1_shared_third"]=propagated_inverse()
    return ans


def stable(x):
    if isinstance(x,dict):
        return {k:stable(v) for k,v in x.items()}
    if isinstance(x,list):
        return [stable(v) for v in x]
    return float(f"{x:.10g}")


if __name__=="__main__":
    cli=argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    result=stable(compute())
    dest=HERE/"independent_results.json"
    if cli.parse_args().check:
        assert result==read("independent_results.json")
        print("PASS: CP-reflected rank-two mass, sharp off-diagonal, light/heavy and current tests")
    else:
        with dest.open("x",encoding="utf8") as f:
            f.write(json.dumps(result,indent=2)+"\n")
