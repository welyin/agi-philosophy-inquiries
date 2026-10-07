"""Working874: finite-CAR positive-root constraint kernel, not a field measure."""
from pathlib import Path
import argparse,json
import numpy as np
HERE=Path(__file__).resolve().parent
TARGET=HERE/'constraint_clock_kernel_probe_results.json'
def annihilation(i):
    c=np.zeros((4,4),complex)
    for bits in range(4):
        if bits>>i&1:c[bits^(1<<i),bits]=(-1)**((bits&((1<<i)-1)).bit_count())
    return c
def lift(h):
    c=[annihilation(i) for i in range(2)]
    return sum(h[i,j]*c[i].conj().T@c[j] for i in range(2) for j in range(2))
def run():
    bg=json.loads((HERE.parent/'861/magnetic_reduced_hamiltonian_results.json').read_text('utf-8'))['original_859_background_clock_reduction'][1]
    a=bg['a'];v=bg['actual_clock_speed'];qstar=v*v/(4*a);hstar=v/(2*a)
    time=1/hstar
    sigma_y=np.array([[0,-1j],[1j,0]],complex)
    angle=.61;rotate=np.cos(angle)*np.eye(2)-1j*np.sin(angle)*sigma_y
    hs=[qstar*np.diag([.08,.12]),qstar*rotate@np.diag([.04,.13])@rotate.conj().T]
    nodes,weights=np.polynomial.legendre.leggauss(64)
    # Unit compact kernel: 3/4 (1-z^2) on [-1,1].
    weights=.75*weights*(1-nodes*nodes);assert abs(sum(weights)-1)<1e-14
    def kernel(B,width=0.,include_FP=True):
        ev,Q=np.linalg.eigh(B)
        if width==0:
            speed=np.sqrt(v*v-4*a*ev)
            f=2*ev/(v+speed)
            val=np.exp(-1j*time*f)
            if not include_FP:val*=v/speed
        else:
            source=ev[:,None]-width*qstar*nodes[None,:]
            speed=np.sqrt(v*v-4*a*source)
            f=2*source/(v+speed)
            integrand=np.exp(-1j*time*f)
            # In z=C(pi), d pi=d z/C_pi. Correct FP cancels it.
            # The wrong version is normalized by common v for comparison.
            if not include_FP:integrand*=v/speed
            val=integrand@weights
        return (Q*val)@Q.conj().T
    B=[lift(h) for h in hs]
    exact=kernel(B[1])@kernel(B[0])
    assert np.linalg.norm(exact.conj().T@exact-np.eye(4))<5e-15
    wrong=kernel(B[1],include_FP=False)@kernel(B[0],include_FP=False)
    wrong_defect=float(np.linalg.norm(wrong.conj().T@wrong-np.eye(4),2))
    assert wrong_defect>.1
    checks=[]
    for width in (.04,.02,.01):
        got=kernel(B[1],width)@kernel(B[0],width)
        err=float(np.linalg.norm(got-exact,2))
        checks.append(dict(delta_kernel_width_fraction=width,ordered_kernel_error=err))
    assert checks[-1]['ordered_kernel_error']<checks[0]['ordered_kernel_error']/12
    order_defect=float(np.linalg.norm(exact-kernel(B[0])@kernel(B[1]),2))
    assert order_defect>1e-4
    m=.5;r=.8
    # Fock basis bits 00,10,01,11; finite preparation from870.
    rho=np.zeros((4,4),complex);rho[0,0]=rho[3,3]=(1-r)/2
    rho[1:3,1:3]=np.array([[r,m],[m,r]])/2
    assert np.linalg.eigvalsh(rho).min()>0 and abs(np.trace(rho)-1)<1e-15
    exact_trace=np.trace(rho@exact)
    sample=kernel(B[1],.01)@kernel(B[0],.01)
    trace_error=float(abs(np.trace(rho@sample)-exact_trace))
    # A source derivative at fixed histories must include the root/Jacobian.
    eps=1e-4
    source_dir=qstar*lift(np.diag([1.,0.]))
    def history(s,width=0):
        return kernel(B[1],width)@kernel(B[0]+s*source_dir,width)
    source_exact=(history(eps)-history(-eps))/(2*eps)
    source_errors=[]
    for width in (.04,.02,.01):
        approx=(history(eps,width)-history(-eps,width))/(2*eps)
        source_errors.append(float(np.linalg.norm(approx-source_exact,2)))
    assert source_errors[-1]<source_errors[0]/12
    return dict(working_round=874,formal_reports=873,cumulative_numbered_groups=3658,
        fresh_numbered_groups=0,all_probe_checks_passed=True,
        original861_alpha=a,original861_clock_speed=v,
        exact_root_unitarity_error=float(np.linalg.norm(exact.conj().T@exact-np.eye(4))),
        compact_delta_regularization=checks,source_insertion_errors=source_errors,
        noncommuting_history_order_defect=order_defect,
        omitted_clock_FP_unitarity_defect=wrong_defect,
        receiver_r=r,receiver_m=m,receiver_history_trace_error=trace_error,
        scope='Finite ordered CAR spectral constraint-kernel calibration with inherited bosonic body coefficients; not a full field measure or original reduced quantization.',
        full_field_measure_or_quantization_proved=False)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    data=run()
    if args.write:
        assert not TARGET.exists()
        TARGET.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert data==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(data,ensure_ascii=False,indent=2))
