"""Independent round 1253 audit; no imports from author or legacy calculations."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
RESULT = HERE / 'independent_results.json'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def check(a, label, tolerance=3e-10):
    value = float(np.linalg.norm(a))
    assert value < tolerance, (label, value)
    return value

def tensor(*factors):
    result = np.array([1], dtype=complex)
    for factor in factors:
        result = np.kron(result, factor)
    return result

def calculate():
    # Independent tensor order: head, a, pointer, memory, L, old_q, physical_R.
    i2, i3 = np.eye(2), np.eye(3)
    y = np.array([[0,-1j],[1j,0]])
    x = np.array([[0,1],[1,0]])
    z, one = i2[:,0], i2[:,1]
    vac = i3[:,0]
    b = -y
    ap, am = (i2+b)/2, (i2-b)/2
    ur = tensor(ap,i2)-1j*tensor(am,y)
    um = tensor(np.outer(z,z),i2)-1j*tensor(np.outer(one,one),x)
    swap = np.zeros((9,9),complex)
    for l in range(3):
        for r in range(3):
            swap[3*r+l,3*l+r]=1
    gates = [tensor(np.eye(8),swap), tensor(ur,i2,i3,i3), tensor(i2,um,i3,i3)]
    prefixes=[np.eye(72,dtype=complex)]
    for gate in gates:
        prefixes.append(gate@prefixes[-1])
    prefix_unitarity=max(check(p.conj().T@p-np.eye(72),'prefix unitary') for p in prefixes)
    # The physical S is independently assembled in four 216-dimensional blocks.
    swap_lr=np.zeros((216,216),complex)
    for a in range(8):
        for l in range(3):
            for q in range(3):
                for r in range(3):
                    col=((a*3+l)*3+q)*3+r
                    row=((a*3+r)*3+q)*3+l
                    swap_lr[row,col]=1
    pm=tensor(np.eye(4),np.outer(one,one),i3,i3,i3)
    physical_swap=[]
    commutator_record=[]
    for prefix in prefixes:
        v=tensor(prefix,i3)
        s=v@swap_lr@v.conj().T
        physical_swap.append(s)
        commutator_record.append(check(pm@s-s@pm,'physical record commutator'))
        check(s@s-np.eye(216),'physical involution')
    commutator_h=[]
    for j,gate in enumerate(gates):
        g=tensor(gate,i3)
        commutator_h.append(check(physical_swap[j+1]@g-g@physical_swap[j], 'off-diagonal H commutator'))
    # No author Material is imported. Auxiliary factors are restored below.
    hc=1.5*np.eye(4)
    for j in range(3):
        hc[j+1,j]=hc[j,j+1]=math.sqrt((j+1)*(3-j))/2
    he,hv=np.linalg.eigh(hc)
    mstar=17+2*math.sqrt(19)
    v9=np.eye(9)-swap
    mass=tensor(hc,np.eye(9))+(mstar+1)*np.eye(36)+tensor(np.eye(4),v9)
    me,mv=np.linalg.eigh(mass)
    check(np.array([me[0]-(mstar+1), me[-1]-(mstar+6)]),'mass endpoints')
    def phi(age):
        return np.exp(-1j*mstar*age)*(hv@(np.exp(-1j*age*he)*(hv.conj().T@np.eye(4)[:,0])))
    aux=tensor(np.array([1,1j])/math.sqrt(2),z,z)
    def encoding(decoded):
        values=decoded.reshape(4,3,3,-1)
        blocks=[]
        for j in range(4):
            raw=np.einsum('a,lrk,b->albrk',aux,values[j],vac).reshape(216,-1)
            blocks.append(tensor(prefixes[j],i3)@raw)
        return np.vstack(blocks)
    delta, star=math.pi/4,math.pi/2
    small_j=np.column_stack([tensor(phi(delta),i3[:,k],vac) for k in range(3)])
    physical_j=encoding(small_j)
    check(physical_j.conj().T@physical_j-np.eye(3),'initial qutrit isometry')
    def evolution(s):
        return encoding(mv@(np.exp(-1j*s*me)[:,None]*(mv.conj().T@small_j)))
    target=np.exp(-1j*star)*encoding(np.column_stack([tensor(phi(delta+star),vac,i3[:,k]) for k in range(3)]))
    all_qutrit_endpoint=check(evolution(star)-target,'all-qutrit true-R endpoint')
    mask=np.tile(np.repeat(np.arange(8)%2,27),4)
    sample_rows=[]
    for s in [math.pi/4,star-F(4,75),star,star+F(4,75),7*math.pi/12]:
        s=float(s)
        output=evolution(s)
        probability=math.sin((delta+s)/2)**6
        clock=check(output.conj().T@(mask[:,None]*output)-probability*np.eye(3),'clock all qutrit')
        nonvac=output[:,1:3].reshape(288,3,2)
        erasure=math.cos(s)**2
        matrix_unit_errors=[]
        for a in range(2):
            for bindex in range(2):
                reduced=np.einsum('dr,ds->rs',nonvac[:,:,a],nonvac[:,:,bindex].conj())
                expected=(1-erasure)*np.outer(i3[:,a+1],i3[:,bindex+1])
                if a==bindex:
                    expected=expected+erasure*np.outer(vac,vac)
                matrix_unit_errors.append(check(reduced-expected,'all-input R channel'))
        # The normalized Choi state retains a passive reference.
        joint=nonvac/math.sqrt(2)
        choi=np.einsum('dri,dsj->risj',joint,joint.conj()).reshape(6,6)
        ideal=np.zeros(6,complex)
        ideal[2]=ideal[5]=1/math.sqrt(2)
        residual=choi-np.outer(ideal,ideal.conj())
        reference_distance=float(np.sum(abs(np.linalg.eigvalsh(residual)))/2)
        check(np.array([reference_distance-erasure]),'reference erasure distance')
        centre=mstar+3.5
        full_distance=float(np.linalg.norm(output-np.exp(-1j*centre*(s-star))*target,2))
        assert full_distance<=2.5*abs(s-star)+3e-10
        sample_rows.append(dict(s=s,clock_probability=probability,clock_residual=clock,
            all_input_R_residual=max(matrix_unit_errors),reference_half_trace_distance=reference_distance,
            complete_isometry_distance=full_distance))
    code=small_j[:,1:3]
    spectral_residuals=[]
    for j,count in enumerate([1,3,4,4,3,1]):
        selected=mv[:,abs(me-(mstar+1+j))<1e-9]
        coefficients=selected.conj().T@code
        spectral_residuals.append(check(coefficients.conj().T@coefficients-count/16*np.eye(2),'source spectrum'))
    interaction=check(code.conj().T@tensor(np.eye(4),v9)@code-np.eye(2),'full interaction energy')
    mean=mstar+3.5
    check(code.conj().T@mass@code-mean*np.eye(2),'full mean')
    check(code.conj().T@((mass-mean*np.eye(36))@(mass-mean*np.eye(36)))@code-1.75*np.eye(2),'full variance')
    algebra=dict(prefix_unitarity_residual=prefix_unitarity,
        record_commutator_residuals=commutator_record, mass_commutator_residuals=commutator_h,
        all_qutrit_endpoint_residual=all_qutrit_endpoint, spectral_weights=[1,3,4,4,3,1],
        spectral_denominator=16,spectral_residuals=spectral_residuals,
        interaction_energy_residual=interaction,mass_min=float(me[0]),mass_max=float(me[-1]),
        sample_rows=sample_rows)
    # Exact analytic certificates, independent of floating matrix diagnostics.
    lower=F(3,8)
    assert F(5,3)**2<3
    endpoint_lower=F(3,64)*(7+4*F(5,3))
    assert endpoint_lower==F(41,64) and endpoint_lower>lower
    # d[p'(a)]/da=(3/2) sin^4(a/2) [5-6 sin^2(a/2)].
    assert F(1,2)<F(5,6)<F(11,12)
    eta,n=F(1,100),30000
    eps=2*eta
    time=eps/lower
    whole=F(5,2)*time
    receiver=time*time
    assert (time,whole,receiver)==(F(4,75),F(2,15),F(16,5625))
    exponent=2*n*eta*eta
    exp_lower=sum((exponent**j)/math.factorial(j) for j in range(19))
    assert exponent==6 and exp_lower>400
    upper_failure=F(2)/exp_lower
    assert upper_failure<F(1,200)
    assert 19<F(9,2)**2 # mstar+6=23+2 sqrt(19)<32.
    assert 32*(n+1)==960032
    budget=dict(derivative_lower=str(lower),right_endpoint_derivative_strict_lower=str(endpoint_lower),
        p_prime_critical_sine_squared='5/6',eta=str(eta),samples=n,
        probability_error=str(eps),time_error=str(time),complete_output_bound=str(whole),
        receiver_bound=str(receiver),hoeffding_exponent=str(exponent),
        failure_strict_upper='1/200',failure_event='acceptance AND time_error_exceeds_4/75',
        accepted_conditional_probability_claim=False,composite_inventory_strict_upper=960032,
        external_facilities_included=False)
    author=json.loads((HERE/'results.json').read_text('utf-8-sig'))
    author_sha=sha(HERE/'check.py')
    assert author['code_sha256']==author_sha
    return dict(status='PASS',independent_groups=2,author_code_sha256=author_sha,
        independent_code_sha256=sha(Path(__file__)),group1_independent_original_gates_and_source=algebra,
        group2_exact_error_and_joint_coverage=budget,
        scope='Given curved fixed-worldline composite; initial boundary and terminal effects are explicit inputs; no dynamic geometry.')

def compare(saved,current):
    if isinstance(current,dict):
        assert saved.keys()==current.keys()
        for key in current: compare(saved[key],current[key])
    elif isinstance(current,list):
        assert len(saved)==len(current)
        for a,b in zip(saved,current): compare(a,b)
    elif isinstance(current,float):
        assert math.isclose(saved,current,rel_tol=2e-8,abs_tol=3e-9),(saved,current)
    else:
        assert saved==current,(saved,current)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-once',action='store_true')
    args=parser.parse_args()
    result=calculate()
    if args.write_once:
        with RESULT.open('x',encoding='utf-8') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
        print('FIRST_WRITE_PASS_1253_INDEPENDENT: 2 groups')
    else:
        compare(json.loads(RESULT.read_text('utf-8-sig')),result)
        print('READ_ONLY_MATCH_1253_INDEPENDENT: 2 groups')
