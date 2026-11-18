"""Independent round 1336 checks. One fixed n-dimensional left code in two
weak-doublet copies; the measured left weak current uses both copies.
Run --write once to save JSON; default execution recomputes and compares.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def adj(a):
    return a.conj().T


def norm(a):
    return np.linalg.norm(a, 2)


def haar(n, rng):
    q, r = np.linalg.qr(rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n)))
    phases = np.diag(r)/abs(np.diag(r))
    return q * phases


def joint_code(s, rng):
    n = len(s)
    s = np.diag(s)
    c = np.diag(np.sqrt(np.maximum(1 - np.diag(s)**2, 0)))
    v = adj(haar(n,rng))
    u0, u1 = haar(n,rng), haar(n,rng)
    return np.block([[u0@s@v, -u0@c], [u1@c@v, u1@s]])


def shape(matrix):
    left, sv, vh = np.linalg.svd(matrix)
    ratio = sv[1]/sv[0]
    return sv[0], ratio, left[:,0], vh[0].conj()


def mass_and_current(code, chi, t, heavy):
    n = len(chi)
    assert n == 2
    n = code.shape[0]//2
    # Original field columns: two weak doublets. Rows: singlet, heavy doublet.
    mat = np.vstack((t*np.kron(np.conj(chi)[None,:],np.eye(n)),
                     heavy * adj(code[:,n:])))
    right, masses, vh = np.linalg.svd(mat)
    keep = np.argsort(masses)[:n]
    return masses[keep], adj(vh)[:,keep], right[:,keep]


def conservative(ru, rd):
    denominator = np.sqrt(1-ru*ru)-rd
    if denominator <= 0:
        return None, None
    eta = (ru+rd)/denominator
    if 1-eta-ru <= 0:
        return eta, None
    return eta, min(1., eta/(1-eta-ru))


class Log:
    def __init__(self):
        self.count, self.excess = {}, {}

    def bound(self,k,lhs,rhs):
        excess = max(0.,float(lhs-rhs))
        if not np.isfinite(lhs) or excess>4e-10*max(1.,rhs):
            raise AssertionError((k,lhs,rhs))
        self.count[k]=self.count.get(k,0)+1
        self.excess[k]=max(excess,self.excess.get(k,0.))


def test_case(log, u, chi, t_up, t_down, heavy):
    n = len(u)//2
    w = u[:,:n]
    a = w[:n]
    ell = w[n:]
    b = np.conj(chi[0])*a + np.conj(chi[1])*ell
    a1, ru, yu, vu = shape(a)
    b1, rd, yd, vd = shape(b)
    c,s = abs(chi)
    eta_max, angular_bound = conservative(ru,rd)
    inp_angle = norm(vd-vu*np.vdot(vu,vd))
    out_angle = norm(yd-yu*np.vdot(yu,yd))
    log.bound("parseval",norm(adj(a)@a+adj(ell)@ell - np.eye(n)),0)
    if n>=3:
        constraints=np.stack((vu.conj(),vd.conj()))
        _,_,vh=np.linalg.svd(constraints,full_matrices=True)
        x=adj(vh)[:,-1]
        log.bound("weak_intersection",max(norm(a@x)-a1*ru,norm(b@x)-b1*rd),0)
        if c>1e-13 and eta_max is not None:
            log.bound("port_angle",s/(c*a1),eta_max)
        if angular_bound is not None:
            log.bound("shared_input_angle",inp_angle,angular_bound)
            log.bound("output_angle",out_angle,angular_bound)
    masses_u, ul, ur = mass_and_current(u,np.array([1.,0.]),t_up,heavy)
    masses_d, dl, dr = mass_and_current(u,chi,t_down,heavy)
    left = adj(ul)@dl
    # A positive Hankel-like contribution from heavy right doublets is possible.
    right_weak = adj(ur[n:])@dr[n:]
    weak = norm(left[:-1,-1])
    R = masses_u[-2]/masses_u[-1]
    alpha = t_up**2/(heavy**2 + t_up**2)
    for r,spectrum,t in ((ru,masses_u,t_up),(rd,masses_d,t_down)):
        ratio = spectrum[-2]/spectrum[-1]
        factor = np.sqrt(1+min(t/heavy,heavy/t)**2)
        log.bound("measured_ratio",r,ratio*factor)
    # Norm of the map from all weak light eigenvectors to W_perp.
    # The CS pairs ensure this equals the largest individual sine.
    up_weak_perp = norm(adj(u[:,n:])@ul[:,:-1])
    down_top_perp = norm(adj(u[:,n:])@dl[:,-1])
    log.bound("finite_left_two_subspaces",weak,inp_angle+up_weak_perp*down_top_perp)
    if t_up<heavy and t_down<heavy:
        shake = lambda t: min(1., t*t/(2*(heavy*heavy-t*t)))
        log.bound("decoupling_bound",weak,inp_angle + shake(t_up)*shake(t_down))
    if 2*alpha*ru**2 < 1:
        du = alpha*ru/(1-2*alpha*ru**2)
        log.bound("weak_mass_admixture_all_ratios",up_weak_perp,du)
        log.bound("full_left_all_ratios",weak,inp_angle+du)
    if 2*ru*ru<1:
        # alpha sqrt(1+min(t/M,M/t)^2)<=1.
        log.bound("admixture_in_observed_ratio",up_weak_perp,R/(1-2*ru**2))
    return dict(ru=float(ru),rd=float(rd),input_angle=float(inp_angle),
                left_mixed_column=float(weak),right_current_norm=float(norm(right_weak)),
                shape_angle_bound=angular_bound)


def calculate():
    rng = np.random.default_rng(1336)
    log=Log()
    cases=0
    witness=None
    for n in (3,4,7,11):
        for a0 in (.18,.61,.96):
            spectrum=a0*np.r_[1.,np.geomspace(.016,.0008,n-1)]
            for slope in (.001,.007,.019,.43):
                phase=np.exp(1j*.71)
                chi=np.array([np.sqrt(1-slope**2),phase*slope])
                u=joint_code(spectrum,rng)
                for t_up,t_down,heavy in ((.2,.3,7.),(.12,.41,1.),(.998,1.003,1.),
                                          (2.,.55,1.),(11.,13.,1.),(.02,.025,1.)):
                    row=test_case(log,u,chi,t_up,t_down,heavy)
                    cases+=1
                    if n==3 and a0==.61 and slope==.007 and t_up==.998:
                        witness=row
        for _ in range(5):
            u=haar(2*n,rng)
            chi=haar(2,rng)[:,0]
            test_case(log,u,chi,.7,.9,1.)
            cases+=1
    # n=2: two rank-one high directions can be orthogonal even for a
    # nonorthogonal pair of ports. The proof's shared weak vector is absent.
    w2=np.array([[1.,0.],[0.,0.],[0.,0.],[0.,1.]])
    u2=np.column_stack((w2,np.eye(4)[:,[1,2]]))
    chi=np.array([.025,np.sqrt(1-.025**2)])
    two=test_case(log,u2,chi,.08,.12,1.)
    if two["input_angle"]<.999:
        raise AssertionError("n=2 exception lost")
    # Historical common-scale central target, not an experimental exclusion.
    R_u,R_d=.619/171.7,.055/2.89
    ru,rd=np.sqrt(2)*np.array([R_u,R_d])
    e,ang=conservative(ru,rd)
    portable=ang+ru/(1-2*ru**2)
    observed=ang+R_u/(1-2*ru**2)
    if observed>=.0416:
        raise AssertionError("common-scale diagnostic")
    return dict(status="PASS",round=1336,shape_family_cases=cases,
                counts=log.count,max_excess=log.excess,witness=witness,
                two_dimensional_exception=two,central_benchmark=dict(
                R_u=R_u,R_d=R_d,rho_u_bound=ru,rho_d_bound=rd,
                eta_bound=e,sine_bound=ang,all_scale_ratio_bound=portable,
                all_scale_observed_ratio_bound=observed,
                historical_left_target=.0416))


def compare(x,y):
    if isinstance(y,dict):
        if x.keys()!=y.keys():
            raise AssertionError("keys differ")
        for k in y:
            compare(x[k],y[k])
    elif isinstance(y,float):
        np.testing.assert_allclose(x,y,atol=3e-10,rtol=2e-8)
    elif x!=y:
        raise AssertionError((x,y))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    p=Path(__file__).with_name("independent_results.json")
    ans=calculate()
    if args.write:
        p.write_text(json.dumps(ans,indent=2)+"\n",encoding="utf-8")
    else:
        compare(ans,json.loads(p.read_text(encoding="utf-8")))
    print(json.dumps(dict(status=ans["status"],tests=sum(ans["counts"].values()),
                          cases=ans["shape_family_cases"],bound=ans["central_benchmark"])))
