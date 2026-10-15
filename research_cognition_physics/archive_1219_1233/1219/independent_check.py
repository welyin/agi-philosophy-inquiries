#!/usr/bin/env python3
"""Independent 1219 finite joint-block audit; never imports author code.
Default replays and checks the saved JSON. --write creates a new result only.
Numerical checks support, but do not replace, the declared operator identities.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse, hashlib, json, math
import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
TOL = 3e-12

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def assert_close(a, b, label):
    err = float(np.max(np.abs(np.asarray(a)-np.asarray(b))))
    if err > TOL:
        raise AssertionError((label, err))
    return err

def bits(k, n=7):
    return [(k >> (n-1-j)) & 1 for j in range(n)]

def ind(b):
    k = 0
    for v in b:
        k = 2*k + v
    return k

def permutation(fn, n=7):
    dim = 2**n
    p = np.zeros((dim,dim), complex)
    for k in range(dim):
        p[ind(fn(bits(k,n))),k] = 1
    return p

def traced_prob(vec, keep, n=7):
    p = np.abs(vec)**2
    out = np.zeros([2]*len(keep))
    for k, val in enumerate(p):
        b=bits(k,n)
        out[tuple(b[j] for j in keep)] += val
    return out

def bare_source():
    vals=[]
    for k in range(128):
        a,da,b,db,rm,re,j=bits(k)
        ma,mb=3+a,3+b
        vals.append((a,b,(5/3)*ma,(5/3)*mb,
                     (4/3)*ma*(1-2*da),(4/3)*mb*(1-2*db)))
    return np.asarray(vals)

def prepare_seed_stock(d, excited_labels):
    # State order seed, target h, stock h; target 0, stock 1.
    n=4*d
    v=np.zeros(n,complex)
    p=np.zeros((n,n),complex)
    charges=np.zeros(n)
    for s in range(d):
        v[4*s+1]=1/math.sqrt(d)
        for a in range(2):
            for b in range(2):
                col=4*s+2*a+b
                aa,bb=(b,a) if s in excited_labels else (a,b)
                p[4*s+2*aa+bb,col]=1
                charges[col]=(3+a)+(3+b)
    out=p@v
    r=np.zeros((2,2),complex)
    for s in range(d):
        for b in range(2):
            u=np.array([out[4*s+2*a+b] for a in range(2)])
            r += np.outer(u,u.conj())
    assert_close(p.conj().T@p,np.eye(n),"preparation unitary")
    assert_close(charges[:,None]*p,p*charges[None,:],"preparation full rest energy")
    # Both target and stock have the same prepared velocity, so all four P commute.
    assert_close(np.vdot(out,charges*out),7,"target plus stock energy")
    return r

def check_full_preparation_extension(d, excited_labels, seed_direction):
    # Target and stock are each four-dimensional h/direction types.
    dim=16*d
    gate=np.zeros((dim,dim),complex)
    energy=np.zeros(dim)
    momentum=np.zeros(dim)
    for z in range(d):
        for t in range(4):
            for c in range(4):
                a,da=divmod(t,2)
                b,db=divmod(c,2)
                col=16*z+4*t+c
                tt,cc=(2*b+da,2*a+db) if z in excited_labels and da==db else (t,c)
                gate[16*z+4*tt+cc,col]=1
                energy[col]=(5/3)*((3+a)+(3+b)+1)
                momentum[col]=(4/3)*((3+a)*(1-2*da)+(3+b)*(1-2*db)+(1-2*seed_direction))
    assert_close(gate.conj().T@gate,np.eye(dim),"full preparation extension")
    assert_close(energy[:,None]*gate,gate*energy[None,:],"full preparation P0")
    assert_close(momentum[:,None]*gate,gate*momentum[None,:],"full preparation Px")
    return dim


def calculate():
    I=np.eye(128,dtype=complex)
    S=permutation(lambda x: x[2:4]+x[0:2]+x[4:])
    Sh=permutation(lambda x: [x[2],x[1],x[0],x[3]]+x[4:])
    R=permutation(lambda x: x[:4]+[x[4]^x[1],x[5]^x[1],x[6]])
    J=permutation(lambda x: x[:6]+[x[6]^x[0]])
    U=(I-1j*S)/math.sqrt(2)
    src=bare_source()
    P0=src[:,2]+src[:,3]
    Px=src[:,4]+src[:,5]
    assert_close(S@S,I,"full type swap involution")
    assert_close(U.conj().T@U,I,"partial exchange unitary")
    for gate in [S,U,R,J]:
        assert_close(P0[:,None]*gate,gate*P0[None,:],"complete P0")
        assert_close(Px[:,None]*gate,gate*Px[None,:],"complete Px")
    assert_close(src[:,2]**2-src[:,4]**2,(3+src[:,0])**2,"A mass shell")
    assert_close(src[:,3]**2-src[:,5]**2,(3+src[:,1])**2,"B mass shell")
    # Source slot basis is reordered as (logical 0, logical 1, vac).
    # The physical old logical input is TWO-dimensional. The vac block is
    # an algebraic extension, not an extra granted source preparation.
    # Initial record embedding retains all 16 material columns.
    B=np.zeros((128,16),complex)
    for j in range(16):
        B[8*j,j]=1
    Vs=[]
    max_gram=0.
    for prog in range(3):
        mid=I if prog==0 else R
        undo=R if prog==2 else I
        blocks=[J@U@undo@mid@(U if q==1 else I)@B for q in range(3)]
        V=np.zeros((384,48),complex)
        for q,block in enumerate(blocks):
            V[128*q:128*(q+1),16*q:16*(q+1)]=block
        max_gram=max(max_gram,assert_close(V.conj().T@V,np.eye(48),"unknown Q/material isometry"))
        Vs.append(blocks[1])
    # Common finite free phase commutes with all declared response gates.
    free=np.exp(-0.37j*P0)
    assert_close(free[:,None]*U,U*free[None,:],"common free phase")

    prep_dims=[check_full_preparation_extension(4,{3},0),
               check_full_preparation_extension(5,{3,4},1)]
    ra=prepare_seed_stock(4,{3})
    rb=prepare_seed_stock(5,{3,4})
    assert_close(ra,np.diag([.75,.25]),"A exact finite preparation")
    assert_close(rb,np.diag([.6,.4]),"B exact finite preparation")

    qa,qb=F(1,4),F(2,5)
    expected=[
        [[0.,.6],[0.,.4]],
        [[.375,.3],[.125,.2]],
        [[0.,.6],[0.,.4]]
    ]
    outputs=[]
    initial=np.array([float(qa),float(qb),
                      (5/3)*(3+float(qa)),(5/3)*(3+float(qb)),
                      (4/3)*(3+float(qa)),-(4/3)*(3+float(qb))])
    for prog,V in enumerate(Vs):
        table=np.zeros((2,2))
        means=np.zeros(6)
        pointer_mismatch=0.
        for a in range(2):
            for b in range(2):
                w=float((qa if a else 1-qa)*(qb if b else 1-qb))
                vec=V[:,ind([a,0,b,1])]
                table += w*traced_prob(vec,[6,1])
                prob=np.abs(vec)**2
                means += w*(prob@src)
                pointer_mismatch += w*sum(prob[k] for k in range(128) if bits(k)[0]!=bits(k)[6])
        assert_close(table,expected[prog],"joint thermal pointer and port")
        assert_close(pointer_mismatch,0,"entity pointer copy")
        delta=means-initial
        assert_close(delta[[2,4]]+delta[[3,5]],[0,0],"A/B complete exchange")
        want=[.15,.075,.15][prog]
        assert_close(delta[0],want,"rest thermal exchange")
        assert_close(delta[2],(5/3)*want,"laboratory energy exchange")
        assert_close(delta[4],[-133/15,-133/30,-133/15][prog],"actual A momentum exchange")
        outputs.append({
            "program":prog,
            "joint_pointer_rows_direction_columns":np.round(table,12).tolist(),
            "delta_h_A":round(float(delta[0]),12),
            "delta_lab_E_A":round(float(delta[2]),12),
            "delta_Px_A":round(float(delta[4]),12),
            "reflection_probability":round(float(table[:,1].sum()),12),
            "pointer1_and_reflection":round(float(table[1,1]),12)
        })
    # Destroy the intermediate records' coherence by an additional environment.
    # Undoing the declared records alone then cannot restore the original interference.
    rho=np.zeros((128,128),complex)
    for a in range(2):
        for b in range(2):
            w=float((qa if a else 1-qa)*(qb if b else 1-qb))
            v=R@U@B[:,ind([a,0,b,1])]
            rho += w*np.outer(v,v.conj())
    for k in range(128):
        for l in range(128):
            if bits(k)[4]!=bits(l)[4]:
                rho[k,l]=0
    finish=J@U@R
    rho=finish@rho@finish.conj().T
    leaked_reflection=float(sum(rho[k,k].real for k in range(128) if bits(k)[1]))
    assert_close(leaked_reflection,.5,"additional environment defeats record-only undo")

    residue=Sh.conj().T@(Px[:,None]*Sh)-np.diag(Px)
    norm=float(np.linalg.norm(residue,2))
    assert_close(norm,8/3,"internal-only swap momentum defect")
    defect_mean=0.
    for a in range(2):
        for b in range(2):
            w=float((qa if a else 1-qa)*(qb if b else 1-qb))
            k=ind([a,0,b,1,0,0,0])
            assert_close(residue[k,k],(8/3)*(b-a),"input-sector defect formula")
            defect_mean+=w*float(residue[k,k].real)
    assert_close(defect_mean,.4,"mean momentum deletion defect")
    assert_close(Sh.conj().T@(P0[:,None]*Sh)-np.diag(P0),0,"internal-only energy still conserved")

    eta=F(1,1000)
    joint_gap=F(1,5)
    lower=joint_gap-2*eta
    assert lower==F(99,500)
    assert 2*F(7)*F(5,3)+2*F(5,3)+4==F(92,3)
    assert F(6,5)*F(5,3)==2
    assert F(6,5)*F(4,3)==F(8,5)
    # No conditional selection is needed; conditional on reflection, the excitation
    # probability happens to be 2/5 for both programs and is not this witness.
    assert_close(outputs[0]["pointer1_and_reflection"]/outputs[0]["reflection_probability"],.4,"conditional value 0")
    assert_close(outputs[1]["pointer1_and_reflection"]/outputs[1]["reflection_probability"],.4,"conditional value 1")

    sources=[
        "archive_1159_1160/1160/proof.md",
        "archive_1106_1119/1118/proof.md",
        "archive_1106_1119/1119/proof.md",
        "archive_1201_1218/1205/proof.md"
    ]
    return {
        "schema":"1219-independent-joint-thermal-port-v1",
        "status":"PASS",
        "all_checks_passed":True,
        "check_groups":6,
        "independent_code_sha256":sha(Path(__file__).resolve()),
        "source_sha256":{p:sha(STAGE.parent/p) for p in sources},
        "groups":{
            "full_type_four_momentum_and_unknown_input":{
                "pass":True,"material_dimension":16,
                "with_record_dimension":128,
                "actual_unknown_input_dimension":2,
                "source_slot_dimension":3,
                "source_slot_basis":["logical0","logical1","vac"],
                "mother_extension_columns":48,
                "actual_input_columns":32,
                "max_gram_error_rounded_12":round(max_gram,12),
                "scope":"All 32 legal logical/material columns plus 16 mathematical vac columns; old 288-dimensional clock dressing is an inherited analytic identity, not reconstructed here."
            },
            "finite_pure_seed_and_excited_stock":{
                "pass":True,"seed_dimensions":[4,5],
                "full_preparation_mother_dimensions":prep_dims,
                "initial_target_h":0,"initial_stock_h":1,
                "marginals":[[.75,.25],[.6,.4]],
                "rest_energy_each_target_plus_stock":7,
                "seed_rest_masses":[1,1],
                "four_stationary_record_program_rest_masses":[1,1,1,1],
                "full_added_prepared_energy":"92/3",
                "scope":"Explicit independent preparation witness with retained seed/stock; no externally supplied Gibbs mixed state and no proof of microscopic manufacture."
            },
            "joint_pointer_port_and_exchange":{"pass":True,"programs":outputs},
            "retained_records_and_future":{
                "pass":True,"extra_environment_then_record_only_undo_reflection":round(leaked_reflection,12),
                "scope":"Deleting inaccessible correlations is a different channel, not a harmless state simplification."
            },
            "internal_only_swap_deletion":{
                "pass":True,"operator_Px_defect_norm":"8/3",
                "input_mean_Px_defect":"2/5",
                "total_lab_energy_defect":0
            },
            "finite_ports_and_conditional_error_budget":{
                "pass":True,
                "proper_flight_time":"6/5","coordinate_flight_time":2,
                "outgoing_x":["8/5","-8/5"],"port_separation":"16/5",
                "joint_event_gap":"1/5","hypothetical_each_program_total_error":"1/1000",
                "gap_lower_bound":"99/500",
                "reflection_conditioned_excitation_both_programs":"2/5",
                "scope":"Error values are a conditional budget, not experimental calibration; joint unconditioned event is the witness."
            }
        },
        "scope_limits":[
            "The internal h references are Gibbs; directed motion, pure seeds/stocks, clock and records are resources, not a Gibbs state of the moving apparatus.",
            "K1219 is an extra effective finite contact permission; this code does not show old H_D alone triggers the new gates.",
            "No independent field reference, acceleration response, Unruh relation, area law or GR is certified."
        ]
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--write",action="store_true")
    args=ap.parse_args()
    result=calculate()
    path=HERE/"independent_results.json"
    if args.write:
        if path.exists():
            raise SystemExit("Refusing to overwrite independent_results.json")
        path.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    saved=json.loads(path.read_text(encoding="utf-8"))
    if saved!=result:
        raise SystemExit("Saved independent result differs from replay")
    print(json.dumps({"status":"PASS","all_checks_passed":True,"results_match":True,
                      "check_groups":result["check_groups"]},ensure_ascii=False))

if __name__=="__main__":
    main()
