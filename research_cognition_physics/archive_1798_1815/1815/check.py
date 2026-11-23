"""One-loop CP checks and the exact finite-portal range of the first analytic invariant."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy

import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/name).read_bytes())


CFG=read("inputs.json")
BASE=runpy.run_path(str(HERE/CFG["gauge_running"]))
OLD=runpy.run_path(str(HERE/CFG["original_soft_model"]))
SPECTRA=read(CFG["gauge_spectra"])
DATA=read(CFG["low_central_masses_and_mixings"])
LOOP=16*math.pi**2


def dagger(x):
    return x.conj().T


def mixing(v):
    s13=abs(v[0,2])
    c13=math.sqrt(1-s13*s13)
    s12,s23=abs(v[0,1])/c13,abs(v[1,2])/c13
    c12,c23=np.sqrt(1-np.array([s12,s23])**2)
    j=(v[0,1]*v[1,2]*v[0,2].conj()*v[1,1].conj()).imag
    sin=j/(s12*c12*s23*c23*s13*c13**2)
    cos=(abs(v[1,0])**2-s12*s12*c23*c23-c12*c12*s23*s23*s13*s13)/(2*s12*c12*s23*c23*s13)
    norm=BASE["FLOW"]["mixing"](np.arcsin([s12,s13,s23]),math.atan2(sin,cos))
    np.testing.assert_allclose(abs(norm),abs(v),atol=3e-13)
    return norm


def run_to_mass(tag, step):
    data=BASE["inputs"](tag)
    point=SPECTRA[tag]["two_doublets_one_triplet"]
    pole=point["M_T_GeV"][0]
    y=data["y"]
    z=np.array([np.diag(y[:3]),data["v"]@np.diag(y[3:6]),np.diag(y[6:])])
    fields=BASE["thresholds"](point,data["MZ"])
    assert min(point["M_D_GeV"])>pole

    def derivative(t,m):
        g=BASE["couplings"](data,fields,t)
        matrices=m@m.conj().transpose(0,2,1)
        tr=np.trace(matrices,axis1=1,axis2=2).real@[3.,3.,1.]
        k=tr-BASE["FLOW"]["GAUGE_DAMPING"]@g
        delta=matrices.copy()
        delta[0]-=matrices[1]
        delta[1]-=matrices[0]
        return (1.5*delta@m+k[:,None,None]*m)/LOOP

    z=BASE["advance"](z,0,math.log(pole/data["MZ"]),fields,derivative,step)
    ul,up,_=np.linalg.svd(z[0])
    dl,down,_=np.linalg.svd(z[1])
    v=mixing(dagger(ul[:,::-1])@dl[:,::-1])
    return pole,up[::-1],down[::-1],v


def complete(h,z):
    real_gram=h+np.outer(z,z.conj())
    assert np.max(abs(real_gram.imag)) < 5e-17
    bare=np.linalg.cholesky(real_gram.real)
    xi=np.linalg.solve(bare,z)
    xi2=float(np.vdot(xi,xi).real)
    assert xi2<1
    np.testing.assert_allclose(bare@(np.eye(3)-np.outer(xi,xi.conj()))@bare.T,h,atol=2e-17)
    return bare,xi,math.sqrt(1-xi2)


def invariants(h,au,z):
    return [(np.vdot(z,mat@au@z)).imag for mat in [h,h@h]]


def null_top(h,au):
    phase=np.exp(-1j*np.angle(h[:,2]))
    h=phase[:,None]*h*phase.conj()
    rho=abs(h[0,1])
    z=np.array([math.sqrt(rho),-h[0,1].conj()/math.sqrt(rho),0.])
    bare,xi,m0=complete(h,z)
    tri=(h[0,1]*h[1,2]*h[2,0]).imag
    fixed=(au[1,1]-au[0,0])*tri
    for scale in [.53,1.,1.87]:
        # Preserve Im Z12 and H with an arbitrary real part of Z12.
        target12=complex(-scale*h[0,1].real,-h[0,1].imag)
        zz=np.array([math.sqrt(rho),target12.conjugate()/math.sqrt(rho),0.])
        complete(h,zz)
        i6,i8=invariants(h,au,zz)
        np.testing.assert_allclose(i8-(h[0,0]+h[1,1]).real*i6,fixed,rtol=1e-10,atol=1e-35)
    i6,i8=invariants(h,au,z)
    assert abs(i6)<1e-28 and abs(i8-fixed)<abs(fixed)*1e-10
    return dict(portal_trace=2*rho, charged_W_longitudinal_up_charm_top_conditional=[.5,.5,0],
                total_doublet_width_over_mass_asymptotic=2*rho/(16*math.pi),
                actual_i6_exact_zero=0.,i8=fixed, light_weak_triangle=tri,
                bare_real_portal_norm=float(np.linalg.norm(bare,2)),bare_over_heavy=m0)


def gram_geometry(h,au):
    b=h.imag
    eig,vec=np.linalg.eigh(1j*b)
    a=eig[-1]
    z0=math.sqrt(2*a)*vec[:,0]
    bare,xi,m0=complete(h,z0)
    plane=np.linalg.svd(b)[0][:,:2]
    projected=plane.T@(b@au-au@b)@plane/2
    slopes,axes=np.linalg.eigh(projected)
    c=slopes[1]
    np.testing.assert_allclose(slopes,[-c,c],atol=1e-18)
    i0,_=invariants(h,au,z0)
    p_zero=math.hypot(2*a,i0/c)
    axes=plane@axes
    orientation=np.sign(axes[:,0]@b@axes[:,1])

    def witness(p,angle):
        gap=math.sqrt(max(0,p*p-4*a*a))
        plus=math.sqrt((p+gap)/2)
        minus=a/plus
        ca,sa=math.cos(angle),math.sin(angle)
        x=plus*(ca*axes[:,0]+sa*axes[:,1])
        y=orientation*minus*(-sa*axes[:,0]+ca*axes[:,1])
        z=x+1j*y
        complete(h,z)
        observed=invariants(h,au,z)[0]
        expected=i0-c*gap*math.cos(2*angle)
        roundoff=40*np.finfo(float).eps*np.linalg.norm(h,2)*np.linalg.norm(au,2)*p
        np.testing.assert_allclose(observed,expected,rtol=1e-9,atol=roundoff)
        return observed
    cap=2*p_zero
    radius=c*math.sqrt(cap*cap-4*a*a)
    for ang in np.linspace(0,math.pi/2,13):
        witness(cap,ang)
    witness(p_zero,0 if i0>0 else math.pi/2)
    return dict(minimum_portal_trace=2*a,degree6_at_minimum=i0,coefficient_c=c,
                amount_at_3loop_unit_factor=i0/LOOP**3,
                minimum_trace_at_i6_zero_in_this_CP_frame=p_zero,
                total_width_over_mass_at_zero_threshold=p_zero/(16*math.pi),
                same_trace_example=cap, attainable_I6_at_same_trace=[i0-radius,i0+radius],
                real_bare_spectral_norm_min_trace=float(np.linalg.norm(bare,2)),
                original_soft_over_pole=float(np.linalg.norm(xi)),real_bare_mass_over_pole=m0)


def least_portal_all_phases(h,au):
    edges=np.array([h[0,1],h[1,2],h[2,0]])
    radii=abs(edges)
    flux=(np.angle(np.prod(edges))+math.pi/2)%math.pi-math.pi/2
    assert 0<abs(flux)<math.pi/4
    lo,hi=0.,min(radii)**2
    for _ in range(80):
        lagrange=(lo+hi)/2
        angles=.5*np.arcsin(lagrange/radii**2)
        if sum(angles)>abs(flux):
            hi=lagrange
        else:
            lo=lagrange
    angles=np.sign(flux)*.5*np.arcsin((lo+hi)/(2*radii**2))
    row_phases=np.r_[0.,np.cumsum(np.angle(edges[:2])-angles[:2])]
    u=np.exp(1j*row_phases)
    hr=u[:,None]*h*u.conj()
    eig,vec=np.linalg.eigh(1j*hr.imag)
    z=math.sqrt(-2*eig[0])*vec[:,0]
    complete(hr,z)
    trace=2*np.linalg.norm(radii*np.sin(angles))
    np.testing.assert_allclose(np.vdot(z,z).real,trace,rtol=2e-12)
    i6,i8=invariants(hr,au,z)
    tol=100*np.finfo(float).eps*np.linalg.norm(hr,2)*trace*np.linalg.norm(au,2)
    assert abs(i6)<tol
    for powers in (1,2,3):
        direction=np.diag(np.diag(au)**powers)
        b,p=hr.imag,hr.real
        derivative=8*sum(b[i,j]*p[i,j]*(direction[i,i]-direction[j,j])
                         for i,j in [(0,1),(1,2),(2,0)])
        assert abs(derivative)<8*tol
        np.testing.assert_allclose(derivative,-8*invariants(hr,direction,z)[0],atol=8*tol)
    eig=np.diag(au)
    differences=eig-np.roll(eig,-1)
    near={}
    for factor in [1.01,1.10,2.]:
        p=factor*trace
        curv_part=sum((radii*differences)**2)/4
        ellipse_part=(p*(max(eig)-min(eig))/4)**2
        near[str(factor)]=math.sqrt((curv_part+ellipse_part)*(p*p-trace*trace))
    return dict(minimum_portal_trace=trace,minimum_asymptotic_total_width_over_mass=trace/(16*math.pi),
                degree6_at_absolute_minimum_exact=0.,degree8_H_squared_at_absolute_minimum=i8,
                abs_I6_upper_if_width_at_most_multiple_of_minimum=near,
                weak_cycle_angle_mod_pi_rad=float(flux),optimal_angles_mod_pi_rad=angles.tolist())


def finite_loop(h,z,up,mass):
    """Test the exact two-vertex flavour products, without assigning any common loop kernel."""
    bare,xi,m0=complete(h,z)
    v=CFG["vev_GeV"]
    mat=np.block([[v/math.sqrt(2)*bare,np.zeros((3,1))],
                  [mass*xi.conj()[None,:],mass*np.array([[m0]])]])
    left,sv,_=np.linalg.svd(mat)
    ck=left[:3,::-1]
    sv=sv[::-1]
    k=dagger(ck)@ck
    radial=k*sv/v
    wprod=abs(ck)**2
    gold_l=(v/math.sqrt(2)*up)[:,None]*ck
    gold_r=-ck*sv
    vals=[radial*radial.T,(1j*radial)*(1j*radial).T,gold_l*gold_r.conj(),
          (k-.13*np.eye(4))*(.13*np.eye(4)).T,wprod]
    for val in vals:
        assert max(abs(val.imag.ravel()))<max(1,np.linalg.norm(val))*2e-15
    return float((ck[0,1]*ck[1,2]*ck[0,2].conj()*ck[1,1].conj()).imag)


def compute(tag):
    if tag=="MZ_1810_reference":
        v=OLD["empirical_matrix"](DATA)
        up,down=(np.array(DATA["yukawa"][key])[:,0] for key in ["up","down"])
        heavy=SPECTRA["common_MZ"]["two_doublets_one_triplet"]["M_T_GeV"][0]
    else:
        heavy,up,down,v=run_to_mass(tag,.12)
        _,up2,down2,v2=run_to_mass(tag,.24)
        for x,y in zip((up,down,abs(v)),(up2,down2,abs(v2))):
            np.testing.assert_allclose(x,y,rtol=2e-9,atol=2e-14)
    h=(v*down**2)@dagger(v)
    au=np.diag(up**2)
    ev,vec=np.linalg.eigh(1j*h.imag)
    z=math.sqrt(-2*ev[0])*vec[:,0]
    blocked=null_top(h,au)
    j=(v[0,1]*v[1,2]*v[0,2].conj()*v[1,1].conj()).imag
    tri=j*(down[1]**2-down[0]**2)*(down[2]**2-down[1]**2)*(down[2]**2-down[0]**2)
    np.testing.assert_allclose(blocked["light_weak_triangle"],tri,rtol=2e-10)
    # A CP basis with one complex link minimizes the portal cost to within the following bounds.
    radii=np.array([abs(h[i,jj]) for i,jj in [(0,1),(1,2),(2,0)]])
    sinflux=abs(tri)/np.prod(radii)
    min_trace_low=2*sinflux/np.linalg.norm(1/radii)
    min_trace_high=2*min(radii)*sinflux
    return dict(heavy_pole_GeV=heavy,up_yukawas=up.tolist(),down_yukawas=down.tolist(),light_J=j,
                leading_decoupling_factor_v2_over_M2=(CFG["vev_GeV"]/heavy)**2,
                tree_witness_J_at_300_GeV=finite_loop(h,z,up,300),
                tree_witness_J_at_3000_GeV=finite_loop(h,z,up,3000),
                all_up_row_phases_min_trace_bracket=[min_trace_low,min_trace_high],
                optimal_across_all_weak_row_phases=least_portal_all_phases(h,au),
                original_phase_convention=gram_geometry(h,au),top_portal_absent=blocked)


def verify(a,b):
    if isinstance(b,dict):
        assert a.keys()==b.keys()
        for k in b: verify(a[k],b[k])
    elif isinstance(b,list):
        assert len(a)==len(b)
        for x,y in zip(a,b): verify(x,y)
    else:
        assert math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-31),(a,b)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--check",action="store_true")
    observed={tag:compute(tag) for tag in ["MZ_1810_reference",*CFG["gauge_inputs"]]}
    path=HERE/"results.json"
    if parser.parse_args().check:
        verify(observed,json.loads(path.read_bytes()))
    else:
        path.write_bytes((json.dumps(observed,indent=2,allow_nan=False)+"\n").encode())
    print("1815 main",hashlib.sha256(path.read_bytes()).hexdigest()[:12])
