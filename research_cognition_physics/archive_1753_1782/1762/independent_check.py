"""Ordered charged-pole propagation on the measured two-dimensional flavor plane."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np


HERE = Path(__file__).resolve().parent
ARCH = HERE.parent
MASS = np.array([.00051099895,.1056583755,1.77693])
RC = math.sqrt(1+.0025/.0000754)
RLOW = math.sqrt(1+.00219/.0000773)
ANGLES = np.array([[.0803,.0899],[.2908,.3164],[.49,.63]])


def polar(a):
    left,s,right = np.linalg.svd(a,full_matrices=False)
    return left@right,s


def get_plane(ub,step):
    prior = runpy.run_path(str(ARCH/"1758"/"check.py"))
    _,d = prior["endpoint"](ub,step)
    eig,v = np.linalg.eigh(d["B"].conj().T@d["B"])
    A = d["L"]@((v/np.sqrt(eig))@v.conj().T)
    return A


def geometry(A):
    V,s = polar(A)
    electron = V[0].conj()
    pe = np.vdot(electron,electron).real
    el = electron/np.sqrt(pe)
    tag = V[2].conj()/np.linalg.norm(V[2])
    dark = np.array([-tag[1].conjugate(),tag[0].conjugate()])
    phase = np.vdot(el,dark)
    dark *= phase.conjugate()/abs(phase)
    basis = np.stack([el,dark],axis=1)
    assert 0<np.vdot(el,dark).real<1 and abs(np.vdot(el,dark).imag)<1e-14
    pmax = (1-math.sqrt(1-ANGLES[0,1]))/(2*pe)
    assert pe<(1+math.sqrt(1-ANGLES[0,1]))/2
    alpha = math.acos(abs(np.vdot(el,dark)))
    tilt = math.asin((s[0]-s[1])/sum(s))
    aligned = math.cos(max(0.,math.acos(math.sqrt(pmax))-alpha-tilt))**2
    isotropic = s[0]**2/sum(s*s)
    bound = max(aligned,isotropic)
    return V,basis,dict(electron_plane_probability=pe,
                      cone_opening_degrees=math.degrees(alpha),
                      max_transport_tilt_degrees=math.degrees(tilt),
                      allowed_max_electron_fraction_within_plane=pmax,
                      max_v_overlap_top_return=aligned,isotropic_max_weight=isotropic,
                      heaviest_weight_upper=bound,conditional_mass_ratio_upper=bound/(1-bound),
                      normalized_H_required_trace_distance=RLOW/(1+RLOW)-bound)


def charged_propagator(V,q):
    f = q*q/(q*q+MASS**2)
    return V.conj().T@(f[:,None]*V)


def compress_tests(A,V,basis):
    Rrequired = RLOW/(1+RLOW)
    h = A.conj().T@A
    kappa = np.linalg.cond(A)
    qmax = math.sqrt((MASS[2]**2-math.sqrt(RLOW)/kappa*MASS[0]**2)/(math.sqrt(RLOW)/kappa-1))
    rows = []
    G = np.eye(2)
    for q in [MASS[0],MASS[1],MASS[2],80.4,5129.332691]:
        F = charged_propagator(V,q)
        ri = np.linalg.inv(basis)@F@basis
        assert np.max(abs(ri.imag))<1e-13
        assert np.min(ri.real)>-1e-12
        vals = np.linalg.eigvalsh(F)
        cap = np.linalg.cond(A@F)**2
        bound = ((q*q+MASS[2]**2)/(q*q+MASS[0]**2))**2*kappa**2
        assert cap<=bound*(1+2e-12)
        G = F@G
        rows.append(dict(q_GeV=q,compressed_gain_ratio=vals[1]/vals[0],
                         arbitrary_seed_mass_ratio_max=cap,
                         spectrum_only_bound=bound))
    T = np.linalg.inv(basis)@(G@G.conj().T)@basis
    assert min(T.real.flat)>0 and max(abs(T.imag).flat)<1e-13
    return dict(required_weight=Rrequired,max_Q_floor_one_stage_by_mass_only_GeV=qmax,
                unfiltered_capacity=np.linalg.eigvalsh(h)[1]/np.linalg.eigvalsh(h)[0],
                one_stage=rows)


def signed_competitor(A,V):
    # One electron-line arm flips sign among three equal normalized paths.
    el = V[0].conj()/np.linalg.norm(V[0])
    P = np.outer(el,el.conj())
    U = np.kron(np.eye(3),np.eye(2))
    U = U.astype(complex)
    U[4:,4:] -= 2*P
    np.testing.assert_allclose(U.conj().T@U,np.eye(6),atol=1e-13)
    read = np.kron(np.ones((3,1))/np.sqrt(3),np.eye(2))
    F = read.conj().T@U@read
    np.testing.assert_allclose(F,np.eye(2)-(2/3)*P,atol=1e-14)
    filtered = A@F
    L,s,Rh = np.linalg.svd(filtered,full_matrices=False)
    h = (s[0]/s[1])**2
    x = math.sqrt((RC+1/RC-2)/(h+1/h-2))
    J0 = np.array([[x,1j*math.sqrt(1-x*x)],[1j*math.sqrt(1-x*x),x]])
    # Set the two eigenmasses with the seed mixing; leave its relative phase free.
    def light(phi):
        ph = np.diag([1.,np.exp(1j*phi)])
        J = Rh.conj().T@ph@J0@ph@Rh.conj()
        m = filtered@J@filtered.T
        l,ms,_ = np.linalg.svd(m)
        e,mu,_ = abs(l[:,0])**2
        obs = np.array([4*e*(1-e),abs(l[0,1])**2/(1-e),mu/(1-e)])
        return J,ms,obs
    # Numerically choose an interior admissible phase; no exhaustiveness claim needed.
    t = np.linspace(0,2*np.pi,4097)
    phase = np.zeros((len(t),2,2),complex)
    phase[:,0,0] = 1
    phase[:,1,1] = np.exp(1j*t)
    m = (L*s)@phase@J0@phase@(L*s).T
    xl,_,_ = np.linalg.svd(m)
    e = abs(xl[:,0,0])**2
    readings = np.column_stack([4*e*(1-e),abs(xl[:,0,1])**2/(1-e),abs(xl[:,1,0])**2/(1-e)])
    fraction = (readings-ANGLES[:,0])/(ANGLES[:,1]-ANGLES[:,0])
    margin = np.min(np.minimum(fraction,1-fraction),axis=1)
    best = int(np.argmax(margin))
    J,ms,obs = light(t[best])
    np.testing.assert_allclose(ms[0]/ms[1],RC,rtol=2e-12)
    np.testing.assert_allclose(J@J.conj().T,np.eye(2),atol=2e-14)
    np.testing.assert_allclose(obs,readings[best],atol=2e-12)
    ew,ev = np.linalg.eigh(F)
    archived = runpy.run_path(str(ARCH/"1761"/"check.py"))
    arcs = archived["phase_solutions"](A@ev[:,::-1],ew[1]/ew[0],RC)
    assert arcs["any_joint_angle_solution"]==bool(margin[best]>0)
    return dict(equal_arm_return_gains=np.linalg.eigvalsh(F),fitted_seed_phase_rad=t[best],
                mass_ratio=ms[0]/ms[1],angles=obs,
                central_frequency_all_angles_joint=arcs["any_joint_angle_solution"],
                minimum_fractional_margin_to_angle_windows=margin[best],
                seed_in_original_polar_frame=dict(real=J.real,imag=J.imag),
                electron_input_visible_return_probability=np.linalg.norm(F@el)**2)


def case(ub,step):
    A = get_plane(ub,step)
    V,basis,g = geometry(A)
    assert g["conditional_mass_ratio_upper"]<1.074
    return dict(ub=ub,geometry=g,positive_compressed_propagation=compress_tests(A,V,basis),
                explicit_three_arm_electron_sign=signed_competitor(A,V))


def tidy(z):
    if isinstance(z,dict): return {k:tidy(v) for k,v in z.items()}
    if isinstance(z,(list,tuple,np.ndarray)): return [tidy(v) for v in z]
    if isinstance(z,(float,np.floating)): return float(format(z,".10g"))
    return z


if __name__=="__main__":
    cli = argparse.ArgumentParser()
    cli.add_argument("--check",action="store_true")
    args = cli.parse_args()
    u = json.loads((ARCH/"1757"/"independent_results.json").read_bytes())["combined_central_solution"]["Vub"]
    out = tidy(dict(masses_GeV=MASS,required_mass_ratio_lower=RLOW,
                    inclusive=case(.0041,.25),profiled=case(u,.25)))
    outfile = HERE/"independent_results.json"
    if args.check:
        assert out==json.loads(outfile.read_bytes())
        print("1762 common charged-pole cone, arbitrary positive cascades and signed three-arm competitor verified")
    else:
        outfile.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n",encoding="utf-8")
