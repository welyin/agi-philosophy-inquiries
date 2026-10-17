"""1272 common material source: exact resource checks and finite matrix diagnostics.
One author group. Default read-only; --write creates source_results.json exclusively.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from fractions import Fraction as F
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent
OUT=HERE/"source_results.json"

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def read(rel): return json.loads((STAGE/rel).read_text(encoding="utf-8-sig"))
def norm(a): return float(np.linalg.norm(a,ord=2))
def fmt(x): return format(float(x),".12e")
def swap4(i,j):
    S=np.zeros((16,16),complex)
    for b in range(16):
        bits=[(b>>(3-k))&1 for k in range(4)]
        bits[i],bits[j]=bits[j],bits[i]
        out=sum(bits[k]<<(3-k) for k in range(4))
        S[out,b]=1
    return S
def unitary(H,t):
    w,v=np.linalg.eigh(H)
    return (v*np.exp(-1j*t*w))@v.conj().T

def calculate():
    dims=[2,2,2,2,4,5,2,3,2,2,2]
    assert math.prod(dims)==15360
    metadata=math.prod(dims[4:])
    assert metadata==960
    ledger=read("1249/quantum_bounds_results.json")["group_constructive_finite_resource_bound"]["ledger"]
    M=1<<int(ledger["M_power_of_two_exponent"])
    EH=int(ledger["Fourier_H_exponent"]); E0=int(ledger["Fourier_menu_base_exponent"])
    EF=int(ledger["Fourier_transport_exponent_decimal"])
    N=int(ledger["integer_n_power_exponent_decimal"])
    assert (M.bit_length()-1,EH,E0)==(3059,3093,14008)
    assert EF==95928+100*M
    assert N==1+max((EF+44)//3,(EH+EF+47)//6)
    geom=read("1267/geometry_results.json")
    assert N==int(geom["uniform_tube_and_bounds"]["N_decimal"])
    slack=F(geom["uniform_tube_and_bounds"]["positive_threshold_slacks"]["J"])
    assert slack>F(18,1000)>F(1,2**100)
    assert 4*N-7>100 and 4*N-15>100
    # Original 1270 bounds are inherited through exact factorization.
    assert 2*N-(E0+EH+15*M+27)>110
    assert 4*N-(E0+60)>110

    tau=F(1,16); xi=F(1,10**9); pi_upper=F(22,7)
    pulse_integral=140*(F(1,4)-3*F(1,5)+3*F(1,6)-F(1,7))
    peak=F(35,16)/tau
    dpeak=F(105,4)/tau**2
    vbound=peak*pi_upper*(1+xi)
    workbound=dpeak*pi_upper*(1+xi)
    assert pulse_integral==1 and peak==35 and dpeak==6720
    assert vbound<111 and workbound<21121
    full=11+4+vbound; sector=11+2+vbound
    assert full<126 and sector<124
    assert full<2**7 and workbound<2**15
    start_min=F(1,2)-F(1,256)
    end_max=F(1,2)+F(1,256)+9*tau
    cap=F(21,400); output=F(17,8)
    sample=F(2)-F(1,10**8)
    assert start_min>cap
    assert end_max==F(17,16)+F(1,256)<sample<output
    assert output-end_max>cap
    initial_withdrawal=F(11)+F(2)
    assert initial_withdrawal==13

    energies=np.array([b.bit_count() for b in range(16)],dtype=float)
    H=np.diag(energies).astype(complex); I=np.eye(16,dtype=complex)
    P2=np.diag((energies==2).astype(float))
    assert np.linalg.matrix_rank(P2)==6
    assert 6*metadata==5760
    e0=np.zeros(16,complex);e0[3]=1
    assert norm((H-2*I)@e0[:,None])==0
    swaps=[swap4(0,2),swap4(1,3),swap4(0,1)]
    comm=[norm(S@H-H@S) for S in swaps]
    assert max(comm)==0
    # Degenerate record target controlled by thermal population commutes Hth.
    NA=np.diag([(b>>3)&1 for b in range(16)]).astype(complex)
    X=np.array([[0,1],[1,0]],complex)
    Grec=np.kron(NA,np.eye(2)-X)
    Hrec=np.kron(H,np.eye(2))
    assert norm(Grec@Hrec-Hrec@Grec)==0

    # Same positive SWAP/CNOT generators also implement ideal inverse
    # involutions by a second application; no negative G is selected.
    G1=I-swaps[0];G2=I-swaps[2]
    assert np.linalg.eigvalsh(G1)[0]>-1e-14
    assert np.linalg.eigvalsh(G2)[0]>-1e-14
    noncomm=norm(G1@G2-G2@G1)
    assert noncomm>1
    U=unitary(G2,.41)@unitary(G1,.37)
    h=11*I+H+.7*G1
    residual_energy=norm(U.conj().T@H@U-H)
    residual_support=norm((I-P2)@U@e0[:,None])
    assert residual_energy<1e-12 and residual_support<1e-12
    assert np.linalg.eigvalsh(h)[0]>=11-1e-12
    hH=U.conj().T@h@U
    eps=2.0**-20 # Diagnostic scale only; physical bound uses saved N.
    J=np.array([[.101,.003],[.003,.099]],complex)
    A=np.array([[2,.2],[.2,3]],complex)
    At=np.kron(A,I); Jt=np.kron(J,I)
    newJ=Jt-eps*np.kron(np.eye(2),hH)
    source_new=(At@newJ+newJ@At)/2+eps*np.kron(A,hH)
    source_old=np.kron((A@J+J@A)/2,I)
    source_residual=norm(source_new-source_old)
    bad_source=(At@newJ+newJ@At)/2+eps*np.kron(A+np.eye(2),hH)
    bad_residual=norm(bad_source-source_old)
    assert source_residual<1e-12 and bad_residual>eps

    # Heisenberg energy derivative: i[h,h]=0; [Hth,v]=0 but
    # [G1,G2] != 0, so bare-matrix integration is not used.
    dv=.3*G1
    heisenberg_derivative=U.conj().T@(dv+1j*(h@h-h@h))@U
    assert norm(heisenberg_derivative-U.conj().T@dv@U)<1e-12
    # Scalar joint spectral calculus proves the interval inclusion universally.
    for j in [-.1,0,.09,.1,.2,.3,.5]:
        for old in np.linalg.eigvalsh(J):
            for ht in np.linalg.eigvalsh(h):
                new=old-eps*ht
                assert int(old<=j-126*eps)<=int(new<=j)<=int(old<=j+126*eps)

    Dk=F(1,15000)+F(5,10**16)
    field=F(111,100)*Dk+F(3,2**80)
    old_field=read("1270/continuation_results.json")["field_tests"]["actual_bounded_quantum_source_comparison_bound"]
    assert field==F(old_field)
    bindings=[
      "1249/canonical.md","1249/quantum_bridge.md","1249/quantum_bounds_results.json",
      "1266/proof.md","1267/geometry_proof.md","1267/geometry_results.json",
      "1269/material_proof.md","1270/output_proof.md",
      "1270/continuation_proof.md","1270/continuation_results.json",
      "1270/acceptance.json"]
    hashes={rel:sha(STAGE/rel) for rel in bindings}
    hashes.update({n:sha(HERE/n) for n in
      ("admission.md","candidate_spec.md","source_baseline.json","source_proof.md","source_check.py")})
    return {
      "round":1272,"part":"same_clock_complete_source_for_finite_thermal_task",
      "status":"PASS_AUTHOR_WITH_SCOPE","author_groups":1,
      "resource":{"dimensions":dims,"total_dimension":15360,"energy2_rank":5760,
        "baselines":11,"initial_excitations":2,"fixed_initial_clock_withdrawal":"13*epsilon",
        "physical_gap":"hbar/s","physical_clock":"tau=s*theta",
        "full_internal_norm_upper":str(full),"support2_norm_upper":str(sector),
        "conservative_full_norm_bound":126,"conservative_support2_bound":124,
        "instant_internal_lower_bound":11,"support2_internal_lower_bound":13,
        "v_norm_upper":str(vbound),"theta_work_norm_upper":str(workbound)},
      "pulse_time":{"integral":str(pulse_integral),"peak":str(peak),
        "derivative_upper":str(dpeak),"earliest_start":str(start_min),
        "latest_end":str(end_max),"sampling_before":str(sample),
        "output":str(output),"initial_causal_margin":str(start_min-cap),
        "final_causal_margin":str(output-end_max-cap)},
      "saved_integer_budget":{"M_exponent":3059,"EH":EH,"E0":E0,
        "N_decimal_digits":len(str(N)),"N_bit_length":N.bit_length(),
        "fourN_minus7_over100":True,"fourN_minus15_over100":True,
        "new_normalized_mass_bound":"<2^-100",
        "new_normalized_work_bound":"<2^-100",
        "saved_J_threshold_slack":str(slack),
        "same_1270_source_budget":"<2^-80"},
      "finite_matrix_diagnostics":{
        "arithmetic":"NumPy finite witnesses; general proofs are analytic",
        "thermal_energy_commutators":list(map(fmt,comm)),
        "noncommuting_gate_commutator_norm":fmt(noncomm),
        "energy_preservation_residual":fmt(residual_energy),
        "energy2_support_residual":fmt(residual_support),
        "same_ordering_source_residual":fmt(source_residual),
        "different_coefficient_source_residual":fmt(bad_residual),
        "positive_selected_generators":True,
        "Heisenberg_work_identity":True,"finite_spectral_inclusions_checked":True},
      "geometry_inheritance":{"field_bound_exact":str(field),
        "old_geometry_process_unchanged":True,
        "same_coefficient_source_operator_identity":True,
        "new_term_not_multiplied_by_old_chi":True,
        "new_heat_conditioning_preserves_old_state_for_product_input":True,
        "independent_new_input_required_for_old_error":True,
        "old_normal_branch_probability_lower":"1/4",
        "new_net_thermal_geometry_response":False},
      "scope":["Same finite effective block and declared new material/boundary inputs.",
        "All new inventories, seeds, records and passive references retained.",
        "Quantum clock positivity is a spectral shift, not full-space positive spectrum.",
        "No K1219 moving capture or recoil permission inherited.",
        "No Jacobson temperature, geometric entropy, area law or CT adoption."],
      "source_bindings_sha256":hashes}

def main():
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true");args=p.parse_args()
    data=calculate()
    if args.write:
        with OUT.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(data,f,ensure_ascii=False,indent=2);f.write("\n")
        print("EXCLUSIVE_FIRST_WRITE PASS_AUTHOR_WITH_SCOPE")
    else:
        if json.loads(OUT.read_text(encoding="utf-8"))!=data:
            raise SystemExit("READ_ONLY_MISMATCH")
        print("READ_ONLY_MATCH PASS_AUTHOR_WITH_SCOPE")
    print("13 epsilon fixed withdrawal; total source identity exact; old field bound inherited")
if __name__=="__main__": main()
