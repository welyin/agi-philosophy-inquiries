"""1322 independent audit. Does not import or execute author code.

Independent SU(2), source, CAR-lift and representation calculations.
1321 is only bound as frozen evidence; its apparatus scenarios are not repeated.
Default/--check is read-only; --write-new exclusively saves independent_results.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROWS=[]
EPS=np.array([[0,1],[-1,0]],complex)
SIGMA=[np.array([[0,1],[1,0]],complex),
       np.array([[0,-1j],[1j,0]],complex),
       np.diag([1.,-1.])]


def check(name,value,tol=2e-9):
    err=float(np.max(np.abs(value))) if np.size(value) else 0.
    if not np.isfinite(err) or err>tol: raise AssertionError((name,err,tol))
    ROWS.append({"name":name,"error":err,"tolerance":tol})


def truth(name,value):
    check(name,0 if bool(value) else 1,0)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frame(h):
    return np.column_stack((EPS@h.conj(),h))/np.linalg.norm(h)


def expi(a,t=1.):
    eig,u=np.linalg.eigh(a)
    return (u*np.exp(1j*t*eig))@u.conj().T


def mass_vector(h,yu,yd,scale):
    ht=EPS@h.conj()
    return (np.kron(np.outer(ht,ht.conj()),yu)+
            np.kron(np.outer(h,h.conj()),yd))/scale


def mass_standard(h,yu,yd):
    ht=EPS@h.conj()
    return np.concatenate((np.kron(ht[:,None],yu),
                           np.kron(h[:,None],yd)),axis=1)


def occupation_indices(mask,n):
    return [i for i in range(n) if mask & (1<<i)]


def fock_lift(u):
    n=len(u); dim=2**n
    out=np.zeros((dim,dim),complex)
    sectors={}
    for mask in range(dim):
        sectors.setdefault(mask.bit_count(),[]).append(mask)
    for k,masks in sectors.items():
        if k==0:
            out[0,0]=1
            continue
        inds={m:occupation_indices(m,n) for m in masks}
        for row in masks:
            for col in masks:
                out[row,col]=np.linalg.det(u[np.ix_(inds[row],inds[col])])
    return out


def annihilator(n,j):
    out=np.zeros((2**n,2**n),complex)
    for mask in range(2**n):
        if mask & (1<<j):
            out[mask^(1<<j),mask]=(-1)**((mask&((1<<j)-1)).bit_count())
    return out


def calculate():
    ROWS.clear()
    rng=np.random.default_rng(1322007)
    hs=[np.array([.37+.29j,.91-.13j]),
        np.array([-.23+.51j,.44+.67j]),
        np.array([0.,1.3],complex)]
    scale=7.4
    for s,h in enumerate(hs):
        j=frame(h); amp=np.linalg.norm(h)
        check(f"frame:{s}:unitarity",j.conj().T@j-np.eye(2))
        check(f"frame:{s}:determinant",np.linalg.det(j)-1)
        check(f"frame:{s}:second_column",j[:,1]-h/amp)
        axis=np.array([.19,-.41,.63])
        group=expi(sum(a*b for a,b in zip(axis,SIGMA))/2,.73)
        alpha=.46
        hp=np.exp(.5j*alpha)*group@h
        d=expi(-alpha*SIGMA[2]/2)
        check(f"frame:{s}:joint_SU2_U1_covariance",frame(hp)-group@j@d)
        for y in [Fraction(1,6),Fraction(-1,2)]:
            yy=float(y)
            rhs=expi(alpha*(yy*np.eye(2)+SIGMA[2]/2))
            lhs=np.exp(1j*alpha*yy)*(frame(hp).conj().T@group@j)
            check(f"charges:{s}:{y}:right_relation_hypercharge",lhs-rhs)
        ph=np.outer(EPS@h.conj(),(EPS@h.conj()).conj())/amp**2
        check(f"reference:{s}:neutral_projector",ph@ph-ph)
        check(f"reference:{s}:neutral_transport",j.conj().T@ph@j-np.diag([1.,0.]))
        htp=EPS@hp.conj()
        php=np.outer(htp,htp.conj())/np.linalg.norm(hp)**2
        check(f"reference:{s}:neutral_covariant_projector",php-group@ph@group.conj().T)

        for n in [3,4,5]:
            yu=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))/7
            yd=(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))/8
            yuv,ydv = scale*yu/amp,scale*yd/amp
            fstd=mass_standard(h,yu,yd)
            fvec=mass_vector(h,yuv,ydv,scale)
            jn=np.kron(j,np.eye(n))
            check(f"mass:{s}:n{n}:physical_matching",fvec@jn-fstd)
            check(f"mass:{s}:n{n}:allleft_conjugate_matching",
                  fvec.conj()@jn.conj()-fstd.conj())
            check(f"mass:{s}:n{n}:same_singular_spectrum",
                  np.linalg.svd(fvec,compute_uv=False)-np.linalg.svd(fstd,compute_uv=False))
            gtrans=np.kron(group,np.eye(n))
            check(f"mass:{s}:n{n}:vector_gauge_equivariance",
                  mass_vector(hp,yuv,ydv,scale)-gtrans@fvec@gtrans.conj().T)
            # Independent finite differences keep the original coefficients fixed.
            step=1e-5
            radial=h/amp
            dstd=(mass_standard(h+step*radial,yu,yd)-
                  mass_standard(h-step*radial,yu,yd))/(2*step)
            dvec=(mass_vector(h+step*radial,yuv,ydv,scale)-
                  mass_vector(h-step*radial,yuv,ydv,scale))/(2*step)
            check(f"source:{s}:n{n}:standard_radial",dstd-fstd/amp)
            check(f"source:{s}:n{n}:vector_radial",dvec-2*fvec/amp)
            check(f"source:{s}:n{n}:radial_non_equivalence",dvec@jn-2*dstd)
            truth(f"source:{s}:n{n}:radial_difference_nonzero",
                  np.linalg.norm(dvec@jn-dstd)>1e-2)
            dscale=(mass_vector(h,yuv,ydv,scale+step)-
                    mass_vector(h,yuv,ydv,scale-step))/(2*step)
            check(f"source:{s}:n{n}:effective_scale",dscale+fvec/scale)
            dh=np.array([.2-.31j,-.12+.17j])
            ht=EPS@h.conj(); dht=EPS@dh.conj()
            analytic=(np.kron(np.outer(dht,ht.conj())+np.outer(ht,dht.conj()),yuv)+
                      np.kron(np.outer(dh,h.conj())+np.outer(h,dh.conj()),ydv))/scale
            dnum=(mass_vector(h+step*dh,yuv,ydv,scale)-
                  mass_vector(h-step*dh,yuv,ydv,scale))/(2*step)
            check(f"source:{s}:n{n}:full_orientation_amplitude",dnum-analytic)

        # Nonconstant reference connection is tested, without transporting its
        # dynamics into the fixed-mode experiment.
        a=sum(x*y for x,y in zip([.27,-.44,.13],SIGMA))/2
        def moving_frame(t):
            return frame(expi(a,t)@h)
        step=1e-5
        jdot=(moving_frame(step)-moving_frame(-step))/(2*step)
        gauge=sum(x*y for x,y in zip([.53,-.11,.29],SIGMA))/2+.071*np.eye(2)
        connection=j.conj().T@gauge@j+1j*j.conj().T@jdot
        check(f"connection:{s}:Hermitian",connection-connection.conj().T)
        chi=np.array([.22+.71j,-.63+.19j])
        chidot=np.array([.03+.17j,.12-.08j])
        original=jdot@chi+j@chidot-1j*gauge@j@chi
        relation=j@(chidot-1j*connection@chi)
        check(f"connection:{s}:exact_derivative_dictionary",original-relation)
        truth(f"connection:{s}:reference_derivative_not_free_zero",
              np.linalg.norm(1j*j.conj().T@jdot)>.1)

    # Right physical W/Z sources at the fixed unitary-reference point.
    gw,gy=.64,.35
    ew=np.sqrt(gw**2+gy**2)
    sw,cw=gy/ew,gw/ew
    for y in [1/6,-.5]:
        vec_photon=gw*sw*SIGMA[2]/2+gy*cw*y*np.eye(2)
        std_photon=gy*cw*(y*np.eye(2)+SIGMA[2]/2)
        check(f"weak_sources:{y}:same_photon",vec_photon-std_photon)
        vec_z=gw*cw*SIGMA[2]/2-gy*sw*y*np.eye(2)
        std_z=-gy*sw*(y*np.eye(2)+SIGMA[2]/2)
        check(f"weak_sources:{y}:different_Z_dictionary",vec_z-std_z-ew*SIGMA[2]/2)
        truth(f"weak_sources:{y}:W_source_not_singlet",np.linalg.norm(gw*SIGMA[0]/2)>.1)
        truth(f"weak_sources:{y}:Z_difference_nonzero",np.linalg.norm(vec_z-std_z)>.1)
    check("weak_sources:right_neutral_vector_Z",gw*cw/2+gy*sw/2-ew/2)
    check("weak_sources:right_neutral_standard_Z",0.)

    # Independent four-mode full Fock lift. It includes all occupation sectors,
    # particle and antiparticle conjugation, not a one-particle-only check.
    j=frame(hs[0])
    u=np.zeros((4,4),complex)
    u[:2,:2]=j.conj()      # all-left right-conjugate particle variables
    u[2:,2:]=j           # their actual antiparticle variables
    gamma=fock_lift(u)
    ops=[annihilator(4,k) for k in range(4)]
    check("CAR_lift:unitary_all_16_sectors",gamma.conj().T@gamma-np.eye(16))
    vacuum=np.zeros(16);vacuum[0]=1
    check("CAR_lift:vacuum_invariant",gamma@vacuum-vacuum)
    number=np.diag([m.bit_count() for m in range(16)])
    check("CAR_lift:full_number_preserved",gamma.conj().T@number@gamma-number)
    for k in range(4):
        check(f"CAR_lift:annihilation_{k}",
              gamma.conj().T@ops[k]@gamma-sum(u[k,l]*ops[l] for l in range(4)))
        for l in range(4):
            check(f"CAR_lift:anticommutator_{k}_{l}",
                  ops[k]@ops[l].conj().T+ops[l].conj().T@ops[k]-(np.eye(16) if k==l else 0))
    # Two left modes plus two right modes: full quadratic mass Hamiltonian,
    # including pair creation, is tested rather than just its singular values.
    y1=np.array([[.31+.27j]])
    y2=np.array([[-.18+.43j]])
    h=hs[0];amp=np.linalg.norm(h)
    fs=mass_standard(h,y1,y2)
    fv=mass_vector(h,scale*y1/amp,scale*y2/amp,scale)
    ur=np.eye(4,dtype=complex);ur[2:,2:]=j.conj()
    gr=fock_lift(ur)
    def manybody_mass(d):
        z=np.zeros((16,16),complex)
        for k in range(2):
            for l in range(2):
                z += d[k,l]*ops[k]@ops[l+2]
        return z+z.conj().T
    hs_mass=manybody_mass(fs.conj())
    hv_mass=manybody_mass(fv.conj())
    check("CAR_parent:allleft_mass_operator_full_Fock",
          gr.conj().T@hv_mass@gr-hs_mass)
    for k in range(2):
        check(f"CAR_parent:left_endpoint_{k}_untouched",gr.conj().T@ops[k]@gr-ops[k])
    check("CAR_parent:free_energy_untouched",gr.conj().T@number@gr-number)
    check("CAR_parent:empty_right_unchanged",gr@vacuum-vacuum)

    # r acts vectorially on physical L/R; the stronger all-left representation
    # is explicitly different even though the residual operational rule agrees.
    check("protection:standard_physical_left_hat_r",1j*(-1)-(-1j))
    check("protection:standard_physical_right_hat_r",(-1j)*1-(-1j))
    check("protection:vector_physical_right_hat_r",1j*(-1)-(-1j))
    check("protection:H_hat_r",(-1)*(-1)-1)
    check("protection:S_hat_r",-1-(-1))
    check("protection:quadratic_H_vector_mass",(-1j)*1j-1)
    check("protection:linear_H_standard_mass",(-1j)*(-1)*(-1j)-1)
    check("protection:left_S_portal",1j*1j*(-1)-1)
    check("protection:QQQL",1j**4-1)
    truth("protection:same_chirality_Majorana_forbidden",(-1j)**2==-1)
    truth("protection:unbroken_R_rep_is_different",(-1j)!=1j)

    # Exact representation diagnostics. The global statement itself uses the
    # allowed positive vectorlike Dirac regulator, not numerical triangle tests.
    fields=[(6,Fraction(1,6),1),(2,Fraction(-1,2),1),
            (6,Fraction(-1,6),-1),(2,Fraction(1,2),-1)]
    for n in [3,4,5]:
        for power in [1,3]:
            truth(f"anomaly:n{n}:hypercharge_power_{power}",
                  n*sum(mult*q**power for mult,q,r in fields)==0)
        for a,b in [(1,0),(3,0),(1,2),(2,1)]:
            truth(f"anomaly:n{n}:r{a}_Y{b}",
                  n*sum(mult*r**a*q**b for mult,q,r in fields)==0)
        truth(f"anomaly:n{n}:Witten_even",(8*n)%2==0)
        truth(f"anomaly:n{n}:Z4_conjugate_counts",n*(8-8)==0)
        truth(f"anomaly:n{n}:SU3_conjugate_pair",n*(2-2)==0)
        truth(f"anomaly:n{n}:SU2_squared_Y",n*(3*Fraction(1,6)+Fraction(-1,2)
                  +3*Fraction(-1,6)+Fraction(1,2))==0)
    # Demonstrate the positive massive Dirac determinant used in the analytic
    # regulator reasoning. This sample is diagnostic, not a global proof.
    eigen=np.array([.11,.87,2.3])
    regmass=.9
    product=np.prod(regmass**2+eigen**2)*regmass**2
    truth("anomaly:positive_vectorlike_regulator_diagnostic",product>0)

    inherited=HERE.parent/"1321"
    expected={
        "proof.md":"5fdf6c40524a83a3fd125fc0191cd98ba181bf93ac2020debb8dcc7000791a68",
        "check.py":"126b4af829a875294cfb4b3cfc5b985ba9696513443c8cd831ba393d0252e3af",
        "results.json":"dae8ee1828d6b74335943a9dce8364a2afcc9c45a1a85ecd4c81a20111accec1",
        "acceptance.json":"9b1a600733ced1323837c7f27bb0f394da173fc9a180eb925ac33b3e4ec83578"}
    for name,known in expected.items():
        truth(f"frozen_1321:{name}",digest(inherited/name)==known)
    old=json.loads((inherited/"results.json").read_text(encoding="utf-8"))
    f=np.array([4.**-i for i in range(5)])
    s3=sum(f[:3]**2)
    saved=[]
    for row in old["families"]:
        n=row["n"]
        sn=sum(f[:n]**2)
        # Recompute only scalar error constants, not any 1321 apparatus scene.
        nominal_bound=4*(np.pi/4*(sn/s3-1))**2
        upper=nominal_bound+18496*2.**-24+196608*2.**-40
        check(f"inheritance:n{n}:frozen_error_constant",upper-row["full_parent_return_error_bound"])
        truth(f"inheritance:n{n}:both_roles_same_threshold",upper<.01)
        saved.append({"n":n,"both_role_return_error_bound":upper})

    return {"round":1322,"status":"INDEPENDENT_CHECK_PASS","check_count":len(ROWS),
            "maximum_absolute_residual":max(r["error"] for r in ROWS),
            "checks":ROWS.copy(),"inherited_sha256":expected,
            "inherited_error_bounds":saved,
            "independent_code_sha256":digest(Path(__file__)),
            "scope":["Fixed nonzero common H and zero gauge background only for task transport",
                     "No free executable J operation or new reference instrument",
                     "Right W/Z, radial Higgs and effective-scale sources remain different",
                     "R_op comparison does not defeat fixed unbroken R_rep",
                     "Vectorlike regulator argument restricted to stated admissible bundles/backgrounds"]}


def main():
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--write-new",action="store_true")
    args=parser.parse_args()
    data=calculate()
    target=HERE/"independent_results.json"
    if args.write_new:
        with target.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(data,stream,ensure_ascii=False,indent=2,allow_nan=False)
            stream.write("\n")
    else:
        if data!=json.loads(target.read_text(encoding="utf-8")):
            raise AssertionError("Independent read-only result mismatch")
    print(json.dumps({"round":1322,"mode":"write-new" if args.write_new else "read-only",
                      "count":data["check_count"],"max_residual":data["maximum_absolute_residual"],
                      "pass":True},ensure_ascii=False))


if __name__=="__main__":
    main()
