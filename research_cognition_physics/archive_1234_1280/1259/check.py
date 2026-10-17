"""1259 finite program-record gate and same-action resource certificate.

Default is read-only reconstruction/comparison. --write exclusively creates
results.json once. The normal-state, source and reference quantifiers are proved
in proof.md; no spatial grid, ODE, field truncation or device simulation is run.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent.parent
SOURCES = [
    HERE / 'proof.md',
    HERE.parent / '1247' / 'proof.md',
    HERE.parent / '1249' / 'canonical.md',
    HERE.parent / '1255' / 'proof.md',
    HERE.parent / '_admission' / 'local_terminal_task_contract_after1256.md',
    HERE.parent / '_admission' / 'finite_terminal_geometry_gate_after1256.md',
    HERE.parent / '_admission' / 'local_collective_record_connection_after1251.md',
]
TOL = 5e-12


def norm(x):
    return float(np.linalg.norm(x, ord=2))


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [encode(v) for v in x]
    return x


def unitary(g):
    values, vectors = np.linalg.eigh(g)
    return (vectors * np.exp(-1j * values)) @ vectors.conj().T


def operators():
    l0 = np.zeros((3, 3), complex)
    l1 = np.zeros((3, 3), complex)
    l0[0, 2], l0[2, 0] = 1j, -1j
    l1[1, 2], l1[2, 1] = 1j, -1j
    angles = [(math.pi/2, 0.0), (0.0, math.pi/2),
              (math.pi/5, 0.0), (0.0, 0.0)]
    return l0, l1, angles


def gate_and_data_reference():
    intervals = [(F(15,16), F(17,16)), (F(31,16), F(33,16))]
    shifts = [F(4), F(8), F(0)]  # flag 0, 1, blank
    images = sorted((l+d, r+d) for d in shifts for l,r in intervals)
    assert all(r-l == F(1,8) for l,r in images)
    assert all(images[i][1] < images[i+1][0] for i in range(5))
    assert images[0][0] == F(15,16) and images[-1][1] == F(161,16)

    l0, l1, angles = operators()
    nk, nd, nf = 4, 4, 3
    ni, no = nk*nd, nk*nd*nf
    generator = np.zeros((no,no),complex)
    injection = np.zeros((no,ni),complex)
    expected = np.zeros_like(injection)
    coefficients = []
    for k,(a,b) in enumerate(angles):
        generator[k*nd*nf:(k+1)*nd*nf,k*nd*nf:(k+1)*nd*nf] = np.kron(np.eye(nd),a*l0+b*l1)
        m = [math.sin(a),math.sin(b),math.cos(a)*math.cos(b)]
        coefficients.append(m)
        for d in range(nd):
            injection[(k*nd+d)*nf+2,k*nd+d] = 1
            for f in range(nf):
                expected[(k*nd+d)*nf+f,k*nd+d] = m[f]
    u = unitary(generator)
    actual = u@injection
    unitary_error = norm(u.conj().T@u-np.eye(no))
    column_error = norm(actual-expected)
    isometry_error = norm(actual.conj().T@actual-np.eye(ni))
    recovery_error = norm(u.conj().T@actual-injection)
    kraus_error = max(norm(actual[a::nf,:]-np.diag(np.repeat(np.array(coefficients)[:,a],nd))) for a in range(nf))
    assert max(unitary_error,column_error,isometry_error,recovery_error,kraus_error) < TOL
    assert norm(generator) <= math.pi/2+TOL

    # Choi purification of every K/data matrix unit, including a 16-dimensional R.
    out = (actual/math.sqrt(ni)).reshape(nk,nd,nf,ni)
    initial = (np.eye(ni)/math.sqrt(ni)).reshape(nk,nd,ni)
    data_r_out = np.einsum('kdfr,kefs->dres',out,out.conj()).reshape(nd*ni,nd*ni)
    data_r_in = np.einsum('kdr,kes->dres',initial,initial.conj()).reshape(nd*ni,nd*ni)
    data_reference_error = norm(data_r_out-data_r_in)
    assert data_reference_error < TOL

    # Decoded B data is invariant. Bare B shifts with the flag and is not fixed.
    bdata = np.array([0.97,1.03,1.97,2.03])
    decoded_b = np.kron(np.eye(nk),np.kron(np.diag(bdata),np.eye(nf)))
    bare_b = decoded_b+np.kron(np.eye(ni),np.diag([4.,8.,0.]))
    decoded_commutator = norm(generator@decoded_b-decoded_b@generator)
    bare_commutator = norm(generator@bare_b-bare_b@generator)
    assert decoded_commutator < TOL and bare_commutator > 6.0
    return dict(status='PASS', intervals=images, represented_input_dimension=ni,
                output_dimension=no, reference_dimension=ni,
                flag_columns=coefficients, generator_norm=norm(generator),
                unitary_error=unitary_error, isometry_error=isometry_error,
                all_kraus_error=kraus_error, column_error=column_error,
                inverse_recovery_error=recovery_error,
                data_reference_marginal_error=data_reference_error,
                decoded_B_commutator=decoded_commutator,
                bare_B_commutator=bare_commutator)


def polynomial_derivative(p):
    return [F(i)*p[i] for i in range(1,len(p))]


def polynomial_value(p,x):
    return sum(c*x**i for i,c in enumerate(p))


def source_and_pulse():
    t = F(1,128)
    # g(theta)=p(u)/T, u=(theta-theta_a)/T.
    p = [F(0),F(0),F(0),F(140),F(-420),F(420),F(-140)]
    dp = polynomial_derivative(p)
    ddp = polynomial_derivative(dp)
    integral = sum(c/F(i+1) for i,c in enumerate(p))
    assert integral == 1
    for polynomial in (p,dp,ddp):
        assert polynomial_value(polynomial,F(0)) == 0
        assert polynomial_value(polynomial,F(1)) == 0
    peak = polynomial_value(p,F(1,2))/t
    assert peak == F(35,16)/t == 280
    eps = F(135005697,4*10**13)
    interaction_coefficients = [eps*c/t for c in p]
    j_coefficients = [-c for c in interaction_coefficients]
    mass_change = [j+v for j,v in zip(j_coefficients,interaction_coefficients)]
    jdot = [c/t for c in polynomial_derivative(j_coefficients)]
    vdot = [c/t for c in polynomial_derivative(interaction_coefficients)]
    assert mass_change == [0]*len(p)
    assert [a+b for a,b in zip(jdot,vdot)] == [0]*len(dp)
    assert polynomial_value(j_coefficients,F(0)) == 0
    assert polynomial_value(j_coefficients,F(1)) == 0

    # Independent finite same-mass check: no Einstein tensor used to define source.
    l0,l1,angles = operators()
    gen = np.zeros((12,12),complex)
    for k,(a,b) in enumerate(angles):
        gen[3*k:3*k+3,3*k:3*k+3] = a*l0+b*l1
    kvals = np.array([-0.3,-0.1,0.2,0.4])
    h0 = np.diag(np.repeat(0.71+0.125+0.5*kvals*kvals,3))
    jin = np.eye(12)*0.01
    m0 = jin+h0
    comm_h0 = norm(gen@h0-h0@gen)
    comm_jin = norm(gen@jin-jin@gen)
    worst_mass = 0.0
    for x in [F(0),F(1,4),F(1,2),F(3,4),F(1)]:
        v = float(eps*polynomial_value(p,x)/t)*gen
        jt = jin-v
        worst_mass = max(worst_mass,norm(jt+h0+v-m0))
    assert max(comm_h0,comm_jin,worst_mass) < TOL
    return dict(status='PASS', pulse_polynomial_u=p, pulse_area=integral,
                endpoint_C2_jets_zero=True, duration=t, maximum_g=peak,
                interaction_mass_coefficients=interaction_coefficients,
                J_compensation_coefficients=j_coefficients,
                total_mass_change_polynomial=mass_change,
                source_identity='T_c=nu_c*M_total*U_c*U_c, from same worldline action',
                H0_G_commutator=comm_h0,Jin_G_commutator=comm_jin,
                finite_mass_identity_error=worst_mass,
                clock_scope='finite internal-time contract; no full quantum-clock preparation or constraint theorem')


def resource_and_failure():
    v,n,hbar,s = F(1000,513**3),F(1),F(1),F(200000)
    eps = hbar/(s*s*n*v)
    assert eps == F(135005697,4*10**13)
    gate_time,flag_time = F(1,128),F(1,128)
    peak_g = F(35,16)/gate_time
    generator_upper = F(11,7)  # strict: pi/2<11/7
    v_upper = eps*peak_g*generator_upper
    assert v_upper == 440*eps < F(3,2000)
    j_floor = F(1,100)-v_upper
    assert j_floor > F(17,2000)
    processing = gate_time+flag_time
    travel_upper = F(35,342)
    access_upper = travel_upper+processing
    access_margin = F(1,8)-access_upper
    assert processing == F(1,64)
    assert access_upper == F(1291,10944)
    assert access_margin == F(77,10944) > 0
    latest_start = F(67,32)
    latest_end_upper = latest_start+access_upper
    assert latest_end_upper == F(24205,10944) < F(71,32) < F(9,4)

    l0,_,_ = operators()
    u = unitary(math.pi/2*l0)
    blank = np.array([0.,0.,1.],complex)
    wrong = np.array([0.,1.,0.],complex)
    p_correct = abs(u@blank)**2
    p_wrong = abs(u@wrong)**2
    tv = float(np.sum(abs(p_correct-p_wrong))/2)
    assert abs(tv-1) < TOL

    # Coherent nonblank input: normalized projection is only a comparator.
    delta = F(1,10000)
    psi = math.sqrt(1-float(delta))*blank+math.sqrt(float(delta))*wrong
    rho = np.outer(psi,psi.conj())
    rho_blank = np.outer(blank,blank.conj())
    trace_distance = float(sum(abs(np.linalg.eigvalsh(rho-rho_blank)))/2)
    assert abs(trace_distance-math.sqrt(float(delta))) < TOL
    p_actual = abs(u@psi)**2
    tv_delta = float(np.sum(abs(p_actual-p_correct))/2)
    assert tv_delta <= math.sqrt(float(delta))+TOL
    return dict(status='PASS', n=n,hbar=hbar,physical_scale=s,cell_volume=v,
                epsilon=eps,gate_time=gate_time,flag_time=flag_time,
                maximum_g=peak_g,generator_strict_upper=generator_upper,
                interaction_strict_upper=v_upper,
                interaction_upper_target=F(3,2000),
                input_J_floor=F(1,100),remaining_J_strict_lower=j_floor,
                J_lower_target=F(17,2000),processing_time=processing,
                travel_strict_upper=travel_upper,access_strict_upper=access_upper,
                access_margin=access_margin,latest_finish_strict_upper=latest_end_upper,
                deletion_correct_probability=p_correct.tolist(),
                deletion_nonblank_probability=p_wrong.tolist(),deletion_TV=tv,
                comparator_delta=delta,projection_trace_distance=trace_distance,
                measured_delta_TV=tv_delta,conservative_blank_error=F(1,50),
                scope='conditional finite block, no 513^3 grid or inherited local-K classification certificate')


def calculate():
    hashes = {str(p.relative_to(RESEARCH)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCES}
    return encode(dict(round=1259,status='PASS',author_groups=3,
                       source_sha256=hashes,
                       code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       numeric_comparison=dict(absolute_tolerance=TOL,relative_tolerance=1e-11),
                       gate_data_reference=gate_and_data_reference(),
                       pulse_same_source=source_and_pulse(),
                       finite_resources_and_deletion=resource_and_failure()))


def compare(a,b,path='result'):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for key in a: compare(a[key],b[key],path+'.'+key)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)): compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert isinstance(b,(int,float)) and math.isclose(a,b,rel_tol=1e-11,abs_tol=TOL),(path,a,b)
    else:
        assert a==b,(path,a,b)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write',action='store_true',help='exclusive first results save')
    args=parser.parse_args()
    result=calculate()
    output=HERE/'results.json'
    if args.write:
        with output.open('x',encoding='utf-8',newline='\n') as stream:
            json.dump(result,stream,ensure_ascii=False,indent=2)
            stream.write('\n')
        print('WRITE_ONCE_PASS: 1259; 3 author groups; results.json created exclusively')
    else:
        compare(result,json.loads(output.read_text(encoding='utf-8')))
        print('READ_ONLY_MATCH: 1259; 3 author groups; no files changed')


if __name__=='__main__':
    main()