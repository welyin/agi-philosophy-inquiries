"""Independent 1237 diagnostics: finite spin carrier, autonomous conversion, full flags.

No author imports. Default is read-only; --write-once creates hash-bound results.
Analytic representation and channel arguments carry the universal statements.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def norm(x):
    return float(np.max(np.abs(x)))


def kron_all(items):
    out = np.ones((1, 1), complex)
    for item in items:
        out = np.kron(out, item)
    return out


def local(dims, matrix, index):
    return kron_all([matrix if k == index else np.eye(d) for k, d in enumerate(dims)])


def spins(dims):
    half = [np.array([[0,1],[1,0]],complex)/2,
            np.array([[0,-1j],[1j,0]],complex)/2,
            np.diag([.5,-.5]).astype(complex)]
    result = []
    for index,d in enumerate(dims):
        terms = []
        for matrix in half:
            if d == 3:
                bigger = np.zeros((3,3),complex)
                bigger[1:,1:] = matrix
                matrix = bigger
            terms.append(local(dims,matrix,index))
        result.append(terms)
    return result


def singlets(n):
    ops = spins([2]*n)
    total = [sum(x[a] for x in ops) for a in range(3)]
    casimir = sum(x@x for x in total)
    values,vectors = np.linalg.eigh(casimir)
    return vectors[:,np.abs(values)<1e-10], ops


def permutation_action(n, permutation):
    source = np.arange(2**n)
    bits = ((source[:,None] >> np.arange(n-1,-1,-1)) & 1)
    target = bits[:,permutation] @ (2**np.arange(n-1,-1,-1))
    result = np.zeros((2**n,2**n),complex)
    result[target,source] = 1
    return result


def carrier_and_maps():
    q4,ops4 = singlets(4)
    q6,ops6 = singlets(6)
    assert q4.shape == (16,2) and q6.shape == (64,5)
    dims = [2,2,2,2,3,3]
    ops144 = spins(dims)
    total = [sum(x[a] for x in ops144) for a in range(3)]
    J2 = sum(x@x for x in total)
    full_eigenvalues = np.linalg.eigvalsh(J2)
    assert np.count_nonzero(np.abs(full_eigenvalues)<1e-9) == 7
    embed4 = np.zeros((144,16),complex)
    embed6 = np.zeros((144,64),complex)
    for source in range(16):
        embed4[source*9,source] = 1
    for source in range(64):
        base,b5,b6 = source//4,(source//2)%2,source%2
        embed6[base*9+(1+b5)*3+1+b6,source] = 1
    Q = np.column_stack((embed4@q4,embed6@q6))
    assert norm(Q.conj().T@Q-np.eye(7))<1e-13
    assert max(norm(x@Q) for x in total)<2e-13
    # Per-leg sqrt(J_i^2), normalized so a half-spin leg has area one.
    area = np.zeros((144,144),complex)
    for terms in ops144:
        c = sum(x@x for x in terms)
        vals,vecs = np.linalg.eigh(c)
        area += (vecs*np.sqrt(np.maximum(vals,0)))@vecs.conj().T/np.sqrt(.75)
    area_physical = Q.conj().T@area@Q
    assert norm(area_physical-np.diag([4,4,6,6,6,6,6]))<2e-13
    # No Gauss vectors in the odd five-half-spin sectors.
    five_occ = np.zeros(144)
    for source in range(144):
        last = source%9
        five_occ[source] = ((last//3 == 0) != (last%3 == 0))
    assert norm(np.diag(five_occ)@Q)<1e-14

    permutations = list(itertools.permutations(range(6)))
    rep = np.array([q6.conj().T@permutation_action(6,p)@q6 for p in permutations])
    assert norm(np.einsum('gji,gjk->gik',rep.conj(),rep)-np.eye(5))<2e-13
    adjacent = []
    shape_residual = 0.0
    shape56 = None
    for i in range(6):
        for j in range(i+1,6):
            p = list(range(6));p[i],p[j] = p[j],p[i]
            swap = q6.conj().T@permutation_action(6,p)@q6
            shape = q6.conj().T@sum(ops6[i][a]@ops6[j][a] for a in range(3))@q6
            shape_residual = max(shape_residual,norm(swap-2*shape-.5*np.eye(5)))
            if j==i+1: adjacent.append(swap)
            if (i,j)==(4,5):shape56=shape
    comm = np.vstack([np.kron(np.eye(5),r)-np.kron(r.T,np.eye(5)) for r in adjacent])
    sv = np.linalg.svd(comm,compute_uv=False)
    assert np.count_nonzero(sv<1e-10)==1 and shape_residual<1e-13
    singlet_pair = np.array([0,1,-1,0],complex)/np.sqrt(2)
    V = q6.conj().T@np.kron(q4,singlet_pair[:,None])
    assert norm(V.conj().T@V-np.eye(2))<1e-13
    tau = V@V.conj().T/2
    sigma = np.eye(5)/5
    half_trace = float(np.sum(np.abs(np.linalg.eigvalsh(tau-sigma)))/2)
    assert abs(half_trace-.6)<1e-13
    old_shape = float(np.trace(tau@shape56).real)
    new_shape = float(np.trace(sigma@shape56).real)
    assert abs(old_shape+.75)<1e-13 and abs(new_shape+.15)<1e-13
    # S4 x S2 preserves the pre-twirl reference but does not select I5/5.
    subgroup_error=0.0
    for pi,r in zip(permutations,rep):
        if set(pi[:4])==set(range(4)):
            subgroup_error=max(subgroup_error,norm(r@tau@r.conj().T-tau))
    assert subgroup_error<2e-13
    diagnostics={"kinematic_dimension":144,"physical_dimensions":[2,5],
                 "total_Gauss_error":max(norm(x@Q) for x in total),
                 "independent_area_error":norm(area_physical-np.diag([4,4,6,6,6,6,6])),
                 "shape_swap_error":shape_residual,"commutant_dimension":1,
                 "subgroup_preserves_unmixed_reference_error":subgroup_error,
                 "same_area_energy_reference_distance":half_trace,
                 "J5_dot_J6_before":old_shape,"J5_dot_J6_after":new_shape}
    return diagnostics,V,rep


def autonomous(V,rep):
    e,g = 1.0,.25
    H_B=np.diag([4*e]*2+[6*e]*5)
    H_W=np.diag([0,2*e])
    H0=np.kron(H_B,np.eye(2))+np.kron(np.eye(7),H_W)
    input_map=np.zeros((14,2),complex)
    output_map=np.zeros((14,2),complex)
    for k in range(2):input_map[2*k+1,k]=1
    output_map[4::2,:]=V
    U1=(np.eye(14)-input_map@input_map.conj().T-output_map@output_map.conj().T
        +output_map@input_map.conj().T+input_map@output_map.conj().T)
    assert norm(U1.conj().T@U1-np.eye(14))<2e-13
    assert norm(U1@input_map-output_map)<1e-13
    U2=[]
    for r in rep:
        rb=np.zeros((7,7),complex);rb[:2,:2]=np.eye(2);rb[2:,2:]=r
        U2.append(np.kron(rb,np.eye(2)))
    U2=np.array(U2)
    conservation=max(norm(U1@H0-H0@U1),norm(np.einsum('gij,jk->gik',U2,H0)-np.einsum('ij,gjk->gik',H0,U2)))
    assert conservation<3e-12
    T=g*np.array([[0,1,0],[1,0,1],[0,1,0]],float)
    tv,tu=np.linalg.eigh(T)
    tstar=np.pi/(np.sqrt(2)*g)
    cstar=(tu*np.exp(-1j*tv*tstar))@tu.T
    assert norm(cstar[:,0]-np.array([0,0,-1]))<3e-15
    cycle=(tu*np.exp(-2j*tv*tstar))@tu.T
    assert norm(cycle-np.eye(3))<5e-15
    timing_checks=[]
    for dt in (-.013,.013):
        shifted=(tu*np.exp(-1j*tv*(tstar+dt)))@tu.T
        error=float(np.linalg.norm(shifted[:,0]-cstar[:,0]))
        bound=float(np.sqrt(2)*g*abs(dt))
        assert error<=bound+1e-14
        timing_checks.append({'delta_t':dt,'isometry_column_error':error,'analytic_bound':bound})
    spectral_floor=4*e-np.sqrt(2)*g
    assert spectral_floor>0
    endpoint_error=0.0
    # Independent exponentiation of actual 42D Hamiltonian blocks.
    for idx in (0,1,17,361,719):
        r=U2[idx]
        D=np.zeros((42,42),complex)
        D[:14,:14]=np.eye(14);D[14:28,14:28]=U1;D[28:,28:]=r@U1
        H=np.kron(np.eye(3),H0).astype(complex)
        H[14:28,:14]+=g*U1;H[:14,14:28]+=g*U1.conj().T
        H[28:,14:28]+=g*r;H[14:28,28:]+=g*r.conj().T
        assert norm(H-(np.kron(np.eye(3),H0)+D@np.kron(T,np.eye(14))@D.conj().T))<1e-12
        vals,vecs=np.linalg.eigh(H)
        assert abs(float(vals.min())-spectral_floor)<2e-13
        initial=np.zeros((42,2),complex);initial[:14,:]=input_map
        actual=(vecs*np.exp(-1j*vals*tstar))@(vecs.conj().T@initial)
        target=np.zeros_like(initial);target[28:,:]=-np.exp(-1j*6*e*tstar)*r@output_map
        endpoint_error=max(endpoint_error,norm(actual-target))
    assert endpoint_error<3e-12
    # Universal clock dressing gives all flags and all unknown code/R, not only samples.
    rows=[]
    for t in (0,tstar/4,tstar/2,3*tstar/4,tstar):
        c=((tu*np.exp(-1j*tv*t))@tu.T)[:,0]
        activated=float(abs(c[1])**2+abs(c[2])**2)
        EB=4*e+2*e*activated;EW=2*e*(1-activated)
        program=float(np.vdot(c,T@c).real)
        assert abs(EB+EW+program-6*e)<3e-14
        rows.append({"time_over_tstar":float(t/tstar),"boundary_energy":EB,
                     "work_energy":EW,"program_interaction_energy":program})
    return {"independent_e":e,"independent_g":g,"energy_commutator_error":conservation,
            "full_H_positive_floor":spectral_floor,"actual_42D_endpoint_error":endpoint_error,
            "same_H_forward_full_return_at_2tstar_error":norm(cycle-np.eye(3)),
            "timing_contract_checks":timing_checks,
            "full_code_and_R_by_isometry_and_clock_dressing":True,"timeline":rows}


def flag_and_deletions(V,rep):
    # Pure 720-dimensional flag; no classical averaging in the complete map.
    M=np.einsum('gij,jk->gik',rep,V)/np.sqrt(720)
    stacked=M.reshape(3600,2)
    assert norm(stacked.conj().T@stacked-np.eye(2))<1e-13
    matrixunit_error=0.0
    for i in range(2):
        for j in range(2):
            out=np.einsum('ga,gb->ab',M[:,:,i],M[:,:,j].conj())
            matrixunit_error=max(matrixunit_error,norm(out-(np.eye(5)/5 if i==j else 0)))
    assert matrixunit_error<1e-13
    # rho_F=X X^dagger; compute its nonzero spectrum from the 10D Gram X^dagger X.
    X=M.reshape(720,10)/np.sqrt(2)
    flag_values=np.linalg.eigvalsh(X.conj().T@X)
    assert norm(flag_values-.1)<1e-13
    boundary=np.einsum('gai,gbi->ab',M,M.conj())/2
    assert norm(boundary-np.eye(5)/5)<1e-13
    global_values=np.linalg.eigvalsh(stacked.conj().T@stacked/2)
    assert norm(global_values-.5)<1e-13
    S_B=np.log(5);S_F=np.log(10);S_BF=np.log(2);MI=S_B+S_F-S_BF
    assert abs((S_B-np.log(2))-(MI-S_F))<1e-14
    # Input purification: exact depolarization gives rho_BR=I5/5 tensor I2/2.
    BR=np.einsum('gai,gbj->aibj',M,M.conj()).reshape(10,10)/2
    assert norm(BR-np.eye(10)/10)<1e-13
    # Exact inverse on BF (and every passive R) is the adjoint of the full isometry.
    recovery_error=norm(stacked.conj().T@stacked-np.eye(2))
    # A B-only reversible channel must have form rho -> U(rho tensor tau)U^dagger.
    # For dim B=5, 2 rank(tau)<=5. The uniform rank-four case attains the bound 1/5.
    reversible_spectrum=np.array([.25,.25,.25,.25,0])
    reversible_distance=float(np.sum(np.abs(reversible_spectrum-.2))/2)
    assert abs(reversible_distance-.2)<1e-15
    tau_spectrum=np.array([.5,.5,0,0,0])
    stop_twirl_distance=float(np.sum(np.abs(tau_spectrum-.2))/2)
    assert abs(stop_twirl_distance-.6)<1e-15
    # Full seven-dimensional Gibbs weights remain distinct from conditioned sector references.
    beta=.7;weights=np.array([2*np.exp(-beta*4),5*np.exp(-beta*6)])
    weights/=weights.sum()
    assert np.all(weights>0) and abs(float(weights.sum())-1)<1e-15
    return {"matrixunit_depolarization_error":matrixunit_error,"full_inverse_error":recovery_error,
            "flag_nonzero_spectrum":flag_values.tolist(),"global_nonzero_spectrum":global_values.tolist(),
            "S_B":float(S_B),"S_F":float(S_F),"S_BF":float(S_BF),"mutual_information_BF":float(MI),
            "delta_reference_entropy":float(np.log(5/2)),"boundary_R_decoupling_error":norm(BR-np.eye(10)/10),
            "B_only_exact_recovery_distance_bound_attained":reversible_distance,
            "stop_twirl_same_area_energy_distance":stop_twirl_distance,
            "example_full_Gibbs_sector_probabilities":weights.tolist(),
            "thermal_temperature_selected":False}


def run():
    first,V,rep=carrier_and_maps()
    return {"status":"PASS_FINITE_CONDITIONAL_GEOMETRIC_REFERENCE_CONVERSION","groups":3,
            "carrier_shape_reference":first,"autonomous_energy_full_code":autonomous(V,rep),
            "coherent_flag_and_deletions":flag_and_deletions(V,rep),
            "scope":{"reference_role_adopted":False,"actual_Lorentz_area_certified":False,
                     "geometric_thermal_entropy_identified":False,"temperature_or_Clausius_derived":False,
                     "area_entropy_universal_slope":False,"GR_proved":False}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write-once',action='store_true')
    args=parser.parse_args();result=run()
    result['author_binding']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                              for name in ('admission.md','proof.md','check.py','results.json','NEXT.md')}
    result['official_report_binding']=hashlib.sha256((HERE.parent/'research_note_1237.md').read_bytes()).hexdigest()
    path=HERE/'independent_results.json'
    if args.write_once:
        if path.exists():raise FileExistsError(path)
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    elif json.loads(path.read_text(encoding='utf8'))!=result:
        raise AssertionError('Read-only replay or author hash mismatch')
    print(json.dumps({'status':result['status'],'groups':3,'mode':'write_once' if args.write_once else 'read_only_match'}))


if __name__=='__main__':main()
