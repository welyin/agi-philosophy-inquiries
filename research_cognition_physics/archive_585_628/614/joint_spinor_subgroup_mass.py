"""614: explicit original subgroup, CAR mass and Spin(10) carrier dictionary.

The group decomposition is standard. The checked interface retains the exact
old masses/charges and distinguishes a carrier from an enlarged gauge theory.
"""
import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
import joint_fermion_gauss_completion as old
import joint_original_mass_spinor_bridge as bridge

HERE=Path(__file__).resolve().parent
TARGET=HERE/"joint_spinor_subgroup_mass_results.json"
STATES=[tuple(c) for n in (0,2,4) for c in itertools.combinations(range(5),n)]
MASKS=[sum(1<<j for j in c) for c in STATES]
EPS=np.array([[0,1],[-1,0]],complex)

def op(a):return float(np.linalg.norm(a,2))

def exterior(U):
    out=np.zeros((16,16),complex)
    for i,A in enumerate(STATES):
        for j,B in enumerate(STATES):
            if len(A)==len(B):
                out[i,j]=np.linalg.det(U[np.ix_(A,B)]) if A else 1.
    return out

def carrier(C,W,z):
    U=np.zeros((5,5),complex);U[:3,:3]=z**-2*C;U[3:,3:]=z**3*W
    return U

def left_rep(C,W,z):
    out=old.representation(C,W,z)[::2,::2].copy()
    for name,sl in old.SLICES.items():
        if name not in bridge.LEFT:
            a,b=sl.start//2,sl.stop//2
            out[a:b,a:b]=out[a:b,a:b].conj()
    return out

def dictionary():
    J=np.zeros((16,16),complex)
    def put(col,state,value=1):
        J[STATES.index(tuple(state)),col]=value
    for a in range(3):
        for alpha in range(2):put(2*a+alpha,(a,3+alpha))
    # colour anti-fundamental in wedge^2 C^3.
    for a in range(3):put(6+a,[j for j in range(3) if j!=a],(-1)**a)
    # full Hodge dual in wedge^4 C^5.
    for a in range(3):put(9+a,[j for j in range(5) if j!=a],(-1)**a)
    for alpha in range(2):
        for beta in range(2):
            if EPS[beta,alpha]:
                a=3+beta
                put(12+alpha,[j for j in range(5) if j!=a],(-1)**a*EPS[beta,alpha])
    put(14,(3,4));put(15,())
    assert op(J.conj().T@J-np.eye(16))==0
    return J

def ph_matrix():
    U=np.zeros((32,32),complex);V=np.zeros_like(U)
    for name,sl in old.SLICES.items():
        if name in bridge.LEFT:U[sl,sl]=np.eye(sl.stop-sl.start)
        else:V[sl,sl]=np.kron(np.eye((sl.stop-sl.start)//2),EPS)
    return np.block([[U,V],[V.conj(),U.conj()]])

def bdg(phi):
    h,d=old.mass_matrices(phi)
    return np.block([[h,d],[d.conj().T,-h.T]])

def subgroup_check():
    J=dictionary();rng=np.random.default_rng(614)
    errors=[];higgs=[]
    for _ in range(18):
        C=old.gauge.group_exp(rng.normal(size=8),3)
        W=old.gauge.group_exp(rng.normal(size=3),2);z=np.exp(1j*rng.normal())
        U=carrier(C,W,z);R=left_rep(C,W,z)
        errors.append(op(exterior(U)@J-J@R))
        H=rng.normal(size=2)+1j*rng.normal(size=2)
        higgs.append(float(np.linalg.norm(U@np.r_[np.zeros(3),H]
                                   -np.r_[np.zeros(3),z**3*W@H])))
    root=np.exp(1j*np.pi/3)
    center=op(carrier(np.exp(2j*np.pi/3)*np.eye(3),-np.eye(2),root)-np.eye(5))
    charges=np.array([-2,-2,-2,3,3])
    exterior_charges=np.array([sum(charges[list(s)]) for s in STATES])
    old_charges=np.array([1]*6+[-4]*3+[2]*3+[-3]*2+[6,0])
    charge_error=op(J.conj().T@np.diag(exterior_charges)@J-np.diag(old_charges))
    assert max(errors+[center,charge_error]+higgs)<3e-13
    return dict(intertwiner_error=max(errors),center_Z6_error=center,charge_dictionary_error=charge_error,
        Higgs_subspace_error=max(higgs),wedge_degrees=[0,2,4],internal_dimension=16,
        all_left_charges=old_charges.tolist(),no_new_gauge_bosons_required_by_restriction=True)

def mass_check():
    J=dictionary();A=np.kron(J,np.eye(2));W=ph_matrix()
    K=np.block([[A,np.zeros_like(A)],[np.zeros_like(A),A.conj()]])@W
    assert op(K@K.conj().T-np.eye(64))<1e-14
    rng=np.random.default_rng(6141);spectrum=[];cov=[];blocks=[];source=[]
    for _ in range(12):
        phi=rng.normal(size=5)*.3
        B=bdg(phi);N=W@B@W.conj().T;D=N[:32,32:]
        M=np.array([[np.trace(EPS.conj().T@D[2*i:2*i+2,2*j:2*j+2])/2
                     for j in range(16)] for i in range(16)])
        blocks.extend([op(N[:32,:32]),op(D-np.kron(M,EPS)),op(M-M.T)])
        target=K@B@K.conj().T
        spectrum.append(float(np.max(abs(np.linalg.eigvalsh(B)-np.linalg.eigvalsh(target)))))
        C=old.gauge.group_exp(rng.normal(size=8),3)
        V=old.gauge.group_exp(rng.normal(size=3),2);z=np.exp(1j*rng.normal())
        H=z**3*V@(phi[:2]+1j*phi[2:4]);newphi=np.r_[H.real,H.imag,phi[4]]
        Rt=np.kron(exterior(carrier(C,V,z)),np.eye(2))
        RtN=np.block([[Rt,np.zeros_like(Rt)],[np.zeros_like(Rt),Rt.conj()]])
        newtarget=K@bdg(newphi)@K.conj().T
        cov.append(op(newtarget-RtN@target@RtN.conj().T))
        for a in range(5):
            step=np.eye(5)[a]*1e-6
            derivative=(bdg(phi+step)-bdg(phi-step))/(2e-6)
            actual=(K@bdg(phi+step)@K.conj().T-K@bdg(phi-step)@K.conj().T)/(2e-6)
            source.append(op(actual-K@derivative@K.conj().T))
    assert max(blocks+spectrum+cov)<5e-13 and max(source)<1e-9
    B=bdg(bridge.PHI);N=W@B@W.conj().T
    nu_pair=np.trace(EPS.conj().T@N[30:32,62:64])/2
    return dict(canonical_error=op(W@W.conj().T-np.eye(64)),
        maximal_mass_block_error=max(blocks),BdG_spectrum_error=max(spectrum),
        original_mass_covariance_error=max(cov),common_scalar_source_error=max(source),
        transformed_neutral_pair=[float(nu_pair.real),float(nu_pair.imag)],
        same_original_complex_Y=True,number_of_CAR_modes_preserved=32,
        finite_mass_interface_not_full_lattice_measure=True)

def clifford():
    annih=[]
    for j in range(5):
        a=np.zeros((32,32),complex)
        for m in range(32):
            if m&(1<<j):
                sign=(-1)**((m&((1<<j)-1)).bit_count())
                a[m^(1<<j),m]=sign
        annih.append(a)
    gamma=[]
    for a in annih:gamma.extend([a+a.conj().T,-1j*(a-a.conj().T)])
    generators=[(-.5j*gamma[a]@gamma[b])[np.ix_(MASKS,MASKS)]
                for a,b in itertools.combinations(range(10),2)]
    return generators,annih

def weyl_dimension(highest):
    rho=[4,3,2,1,0];dimension=Fraction(1)
    lam=[Fraction(x) for x in highest]
    for i,j in itertools.combinations(range(5),2):
        dimension*=((lam[i]+rho[i])**2-(lam[j]+rho[j])**2)/Fraction(rho[i]**2-rho[j]**2)
    assert dimension.denominator==1
    return dimension.numerator

def enlargement_check():
    generators,annih=clifford()
    B=np.zeros((16,16),complex);B[0,0]=1
    # No new fermionic spin dimension: these are internal Spin(10) matrices.
    residuals=[T@B+B@T.T for T in generators]
    constraints=np.column_stack([np.r_[x.real.ravel(),x.imag.ravel()] for x in residuals])
    rank=int(np.linalg.matrix_rank(constraints,tol=1e-12))
    assert rank==21
    X=np.diag([2*len(s)-5 for s in STATES])
    defect=X@B+B@X.T
    assert op(defect+10*B)==0
    su5=[]
    for a,b in itertools.combinations(range(5),2):
        A=annih[a].conj().T@annih[b]
        su5.extend([(A+A.conj().T)[np.ix_(MASKS,MASKS)],
                    (-1j*(A-A.conj().T))[np.ix_(MASKS,MASKS)]])
    for a in range(4):
        t=annih[a].conj().T@annih[a]-annih[a+1].conj().T@annih[a+1]
        su5.append(t[np.ix_(MASKS,MASKS)])
    assert len(su5)==24
    su5_res=max(op(T@B+B@T.T) for T in su5)
    assert su5_res==0
    dim16=weyl_dimension(["1/2"]*4+["-1/2"])
    dim126=weyl_dimension([1,1,1,1,-1])
    assert (dim16,dim126)==(16,126)
    return dict(full_Spin10_generator_count=45,Majorana_stabilizer_constraint_rank=rank,
        stabilizer_dimension=45-rank,SU5_Majorana_error=su5_res,
        extra_U1_Majorana_charge=-10,full_group_singlet_Ward_defect=op(defect),
        exact_Weyl_dimension_spinor=dim16,exact_Weyl_dimension_symmetric_extremal_square=dim126,
        gauging_full_carrier_with_unchanged_real_singlet_fails=True,
        subgroup_restriction_has_no_such_Majorana_obstruction=True,
        no_physical_126_scalar_required_for_subgroup_only=True)

def run():
    deps=("research_note_531.md","research_note_598.md","research_note_611.md","research_note_613.md",
          "joint_fermion_gauss_completion.py","joint_original_mass_spinor_bridge.py")
    return dict(round=614,tests_run=3,failures=0,errors=0,
        subgroup=subgroup_check(),mass=mass_check(),enlargement=enlargement_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_group_global_quotient_preserved=True,
            original_finite_CAR_mass_and_sources_preserved=True,
            Spin10_carrier_distinguished_from_physical_gauge_enlargement=True,
            full_Spin10_unchanged_scalar_branch_rejected=True,
            chiral_measure_locality_reconstruction_GR_and_observations_still_open=True))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write-results",action="store_true")
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open("x",encoding="utf8",newline="\n") as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    else:assert json.loads(TARGET.read_text("utf8"))==result
    print(json.dumps(result,ensure_ascii=False))
