"""1275: full Gibbs two-read calibration; one scoped author group.
Reuses only the frozen 1237 carrier construction. --write is exclusive.
"""
from __future__ import annotations
import os
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import argparse
import hashlib
import importlib.util
import json
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "1237" / "check.py"
TOL = 3e-10

def need(ok, label):
    if not bool(ok):
        raise AssertionError(label)

def norm(a):
    return float(np.linalg.norm(a))

def qstr(q):
    return f"{q.numerator}/{q.denominator}"

def exact_certificate():
    # In units 4 Z_C, a_C has numerator 2+2C and Re c(t_m) has 1.
    difference = [F(2)-F(1), F(2)]
    quarter_denominator = [F(4)/4, F(8)/4]
    need(difference == quarter_denominator, "all C: a_C-Re c(tm)=1/4")
    b = F(2, 5)  # exp(-2 beta_star e)
    p = F(2)/(2+5*b)
    w = 2*b / ((2+5*b)*(1+b))
    p11 = p-w/2
    joint = [[1-2*p+p11, p-p11], [p-p11, p11]]
    bare = [[F(1, 2), F(0)], [F(0), F(1, 2)]]
    tv = sum(abs(joint[a][b]-bare[a][b]) for a in range(2) for b in range(2))/2
    robust = tv-F(1, 1000)-F(1, 1000)
    need(w == F(1,7) and p == F(1,2), "thermal block weights")
    need(joint == [[F(3,7),F(1,14)],[F(1,14),F(3,7)]], "exact joint table")
    need(tv == F(1,7) and robust == F(493,3500), "TV arithmetic")
    return dict(beta_star_boltzmann=qstr(b),p_B=qstr(p),ST_weight_each=qstr(w),
                joint_00_01_10_11=[[qstr(x) for x in row] for row in joint],
                bare_00_01_10_11=[[qstr(x) for x in row] for row in bare],
                total_variation=qstr(tv),
                epsilon_actual=qstr(F(1,1000)),epsilon_bare=qstr(F(1,1000)),
                robust_lower_bound=qstr(robust),
                quarter_identity_polynomial_coefficients=[qstr(x) for x in difference])

def analytic(beta, g, angle):
    z = beta*np.sqrt(2)*g
    ch, sh = np.cosh(z), np.sinh(z)
    zc = 1+2*ch
    p = 2/(2+5*np.exp(-2*beta))
    w = 2*np.exp(-2*beta)/((2+5*np.exp(-2*beta))*(1+np.exp(-2*beta)))
    ac = (1+ch)/(2*zc)
    weights = np.array([(ch/8+1/4)/zc,
                        (1+np.exp(-z))/(8*zc),(1+np.exp(z))/(8*zc),
                        np.exp(-z)/(16*zc),np.exp(z)/(16*zc)])
    frequencies = np.array([0,1,-1,2,-2])
    ct = np.sum(weights*np.exp(1j*frequencies*angle))
    p11 = p-2*w*ac+2*w*ct.real
    joint = np.array([[1-2*p+p11,p-p11],[p-p11,p11]])
    energy = 2*w*np.sqrt(2)*g*sh/zc
    return joint, energy, p, p-2*w*ac+2*w*ct, weights

def calculate():
    spec = importlib.util.spec_from_file_location("frozen_1237_carrier_1275",OLD)
    old = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(old)
    _, wmap, permutations, _ = old.carrier()
    hb = np.array([4.,4.,6.,6.,6.,6.,6.])
    h0diag = (hb[:,None]+np.array([0.,2.])).ravel()
    h0 = np.diag(h0diag)
    source, target = np.zeros((14,2),complex), np.zeros((14,2),complex)
    source[1,0],source[3,1]=1,1
    target[4::2]=wmap
    sproj, tproj = source@source.conj().T,target@target.conj().T
    u1=np.eye(14)-sproj-tproj+source@target.conj().T+target@source.conj().T
    p14=np.diag(np.repeat([1.,1.,0.,0.,0.,0.,0.],2))
    p=np.kron(np.eye(3),p14)
    projs=[np.eye(42)-p,p]
    mask=np.diag(p).real
    dmask=mask[:,None]*mask[None,:]+(1-mask[:,None])*(1-mask[None,:])
    ac=np.array([[0.,1.,0.],[1.,0.,1.],[0.,1.,0.]])
    fullh0=np.kron(np.eye(3),h0)
    beta_star=float(np.log(2.5)/2)
    betas=[.2,beta_star,.9]
    angles=[0.,np.pi/4,np.pi/2,np.pi,2*np.pi]
    g=1.
    residuals={k:0. for k in ["U1_unitarity","gate_energy","path_decoration",
        "transformed_projector","Gibbs_marginal","joint_table","row_marginal",
        "column_marginal","normalization","nonselective_bare_energy",
        "measurement_energy","complex_correlation","detailed_balance","block_spectrum"]}
    sums=np.zeros((len(betas),len(angles),2,2))
    energy_sums=np.zeros(len(betas))
    def report(key, value):
        residuals[key]=max(residuals[key],float(value))
    report("U1_unitarity",norm(u1.conj().T@u1-np.eye(14)))
    report("gate_energy",norm(u1@h0-h0@u1))
    expected_spectrum=np.sort(np.concatenate([h0diag-np.sqrt(2)*g,h0diag,h0diag+np.sqrt(2)*g]))
    for perm in permutations:
        u2b=np.zeros((7,7),complex)
        u2b[:2,:2]=np.eye(2)
        u2b[2:,2:]=perm
        u2=np.kron(u2b,np.eye(2))
        h=fullh0.astype(complex).copy()
        h[14:28,:14]=g*u1
        h[:14,14:28]=g*u1.conj().T
        h[28:42,14:28]=g*u2
        h[14:28,28:42]=g*u2.conj().T
        d=np.zeros((42,42),complex)
        d[:14,:14]=np.eye(14)
        d[14:28,14:28]=u1
        d[28:,28:]=u2@u1
        report("gate_energy",norm(u2@h0-h0@u2))
        report("path_decoration",norm(d.conj().T@h@d-(fullh0+g*np.kron(ac,np.eye(14)))))
        q=p14-sproj+tproj
        pt=np.kron(np.eye(3),q)
        pt[:14,:14]+=sproj-tproj
        report("transformed_projector",norm(d.conj().T@p@d-pt))
        vals,vec=np.linalg.eigh(h)
        report("block_spectrum",norm(vals-expected_spectrum))
        us=[(vec*np.exp(-1j*vals*angle/(np.sqrt(2)*g)))@vec.conj().T for angle in angles]
        for bi,beta in enumerate(betas):
            weights=np.exp(-beta*(vals-vals.min()))
            weights/=weights.sum()
            rho=(vec*weights)@vec.conj().T
            marg=np.einsum("cbwcdw->bd",rho.reshape(3,7,2,3,7,2))
            rb=np.exp(-beta*(hb-hb.min()));rb/=rb.sum()
            report("Gibbs_marginal",norm(marg-np.diag(rb)))
            dephased=dmask*rho
            energy=float(np.trace(h@(dephased-rho)).real)
            report("nonselective_bare_energy",abs(np.trace(fullh0@(dephased-rho))))
            expected_energy=analytic(beta,g,0)[1]
            report("measurement_energy",abs(energy-expected_energy))
            energy_sums[bi]+=energy
            for ti,(angle,u) in enumerate(zip(angles,us)):
                joint=np.zeros((2,2))
                for first in range(2):
                    post=projs[first]@rho@projs[first]
                    evolved=u@post@u.conj().T
                    for last in range(2):
                        joint[first,last]=np.trace(projs[last]@evolved).real
                expected,_,pb,cp,sw=analytic(beta,g,angle)
                report("joint_table",norm(joint-expected))
                report("normalization",abs(joint.sum()-1))
                report("row_marginal",norm(joint.sum(axis=1)-[1-pb,pb]))
                report("column_marginal",norm(joint.sum(axis=0)-[1-pb,pb]))
                measured_c=np.trace(rho@u.conj().T@p@u@p)
                report("complex_correlation",abs(measured_c-cp))
                report("detailed_balance",abs(sw[1]/sw[2]-np.exp(-beta*np.sqrt(2)*g)))
                report("detailed_balance",abs(sw[3]/sw[4]-np.exp(-2*beta*np.sqrt(2)*g)))
                sums[bi,ti]+=joint
    count=len(permutations)
    sums/=count
    energy_sums/=count
    need(count==720,"all original permutation flag blocks")
    need(max(residuals.values())<TOL,"all matrix checks")
    rows=[]
    for bi,beta in enumerate(betas):
        rows.append(dict(beta=beta,g=g,e=1.,initial_measurement_energy=float(energy_sums[bi]),
            angle_samples=[dict(omega_t=float(angle),joint=sums[bi,ti].tolist()) for ti,angle in enumerate(angles)]))
    matched=sums[1,2]
    bare=np.diag([.5,.5])
    tv=float(np.sum(np.abs(matched-bare))/2)
    need(abs(tv-1/7)<TOL,"matched finite gap")
    return dict(schema="research_round_1275_v1",status="PASS",author_groups=1,
        source_construction=dict(old_carrier="1237/check.py:carrier",
            kinematic_dimension=144,physical_boundary_dimension=7,work_dimension=2,
            flag_dimension=720,program_dimension=3,full_matrix_block_dimension=42,
            direct_Gibbs_blocks=count,inherited_carrier_recounted=False),
        exact=exact_certificate(),diagnostics=rows,residuals=residuals,
        matched_reference=dict(beta_star=beta_star,omega=np.sqrt(2),time=np.pi/(2*np.sqrt(2)),
            joint=matched.tolist(),isolated_boundary_joint=bare.tolist(),
            total_variation=tv,initial_measurement_energy=float(energy_sums[1])),
        scope=dict(full_interacting_Gibbs=True,full_Gibbs_KMS_refuted=False,
            continuous_parameters_proved_analytically=True,numerical_samples_are_proof=False,
            ideal_Luders_contract=True,actual_intermediate_instrument_certified=False,
            actual_reference_preparation_certified=False,measurement_free=False,
            first_measurement_energy_is_minimum_work=False,
            R_marginal_unchanged=True,
            unknown_input_reversible_preservation_certified=False,
            input_correlations_may_change=True,
            fixed_table_for_arbitrary_input=False,
            source_certificate_for_instrument=False,observer_time_certified=False,
            Lorentz_area_certified=False,geometry_generation=False,
            user_conjecture_refuted=False))

def compare(a,b,path="results"):
    if isinstance(a,dict):
        need(isinstance(b,dict) and a.keys()==b.keys(),path+" keys")
        for k in a:compare(a[k],b[k],path+"."+k)
    elif isinstance(a,list):
        need(isinstance(b,list) and len(a)==len(b),path+" length")
        for j,(x,y) in enumerate(zip(a,b)):compare(x,y,path+"["+str(j)+"]")
    elif isinstance(a,float):
        need(isinstance(b,(int,float)) and np.isclose(a,b,rtol=1e-10,atol=1e-10),path+" float")
    else:
        need(a==b,path+" exact")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--write",action="store_true")
    args=parser.parse_args()
    result=calculate()
    paths=[HERE/"admission.md",HERE/"proof.md",HERE/"check.py",OLD,
           OLD.with_name("proof.md"),OLD.with_name("results.json"),HERE.parent/"1274/proof.md"]
    result["source_sha256"]={p.relative_to(HERE.parent).as_posix():
        hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    target=HERE/"results.json"
    if args.write:
        with target.open("x",encoding="utf-8",newline="\n") as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write("\n")
        mode="CREATED_EXCLUSIVE"
    else:
        compare(result,json.loads(target.read_text(encoding="utf-8")))
        mode="READ_ONLY_MATCH"
    print(json.dumps(dict(status="PASS",mode=mode,author_groups=1,
        total_variation=result["matched_reference"]["total_variation"],
        first_measurement_energy=result["matched_reference"]["initial_measurement_energy"],
        maximum_residual=max(result["residuals"].values())),ensure_ascii=False))
if __name__=="__main__":
    main()
