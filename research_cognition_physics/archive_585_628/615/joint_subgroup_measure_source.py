"""615: exact restricted overlap measure and source on a declared flat slice.

One site in four Euclidean directions; antiperiodic fermions in direction 4.
The original subgroup, 16 charges, m0=1 and mature S9 auxiliary prescription
are retained. This is not a Hamiltonian reconstruction or an interacting SM.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
import joint_spinor_subgroup_mass as old
import joint_chiral_fibre_source as spin

HERE=Path(__file__).resolve().parent
TARGET=HERE/'joint_subgroup_measure_source_results.json'
Q5=np.array([-2,-2,-2,3,3])
Q=np.array([sum(Q5[list(s)]) if s else 0 for s in old.STATES])
COEFFICIENTS=[Fraction((k+1)*(k+2)*(9-k),990) for k in range(9)]

def op(a):return float(np.linalg.norm(a,2))

def pfaffian(matrix):
    """Skew elimination with simultaneous row/column pivoting."""
    a=np.array(matrix,dtype=complex,copy=True);value=1+0j
    assert len(a)%2==0 and np.max(abs(a+a.T))<2e-12
    for k in range(0,len(a),2):
        j=k+1+int(np.argmax(abs(a[k,k+1:])))
        if abs(a[k,j])<1e-28:return 0j
        if j!=k+1:
            a[[k+1,j],:]=a[[j,k+1],:]
            a[:,[k+1,j]]=a[:,[j,k+1]]
            value=-value
        pivot=a[k,k+1];value*=pivot
        x=a[k,k+2:].copy();y=a[k+1,k+2:].copy()
        a[k+2:,k+2:]+=(np.outer(y,x)-np.outer(x,y))/pivot
    return value

def internal_pairing():
    _,ann=old.clifford();gamma=[]
    for a in ann:gamma.extend([a+a.conj().T,-1j*(a-a.conj().T)])
    C=np.eye(32,dtype=complex)
    for j in (1,3,5,7,9):C=C@gamma[j]
    return [(C@g)[np.ix_(old.MASKS,old.MASKS)] for g in gamma]

T=internal_pairing()
B=1j*spin.G5@spin.GAMMA[1]@spin.GAMMA[3]
eig,vec=np.linalg.eigh(spin.G5)
VP=vec[:,eig>.5];VM=vec[:,eig<-.5]
EPS=VP.T@B@VP

def frames(theta):
    u=np.zeros((64,32),complex);v=np.zeros_like(u)
    D=np.zeros((64,64),complex);H=np.zeros_like(D)
    for i,q in enumerate(Q):
        x=np.pi+q*theta
        rot=np.cos(x/2)*np.eye(4)+1j*np.sin(x/2)*spin.GAMMA[3]
        X=-np.cos(x)*np.eye(4)+1j*np.sin(x)*spin.GAMMA[3]
        u[4*i:4*i+4,2*i:2*i+2]=rot@VP
        v[4*i:4*i+4,2*i:2*i+2]=rot@VM
        D[4*i:4*i+4,4*i:4*i+4]=(np.eye(4)+X)/2
        H[4*i:4*i+4,4*i:4*i+4]=spin.G5@X
    return u,v,D,H

def pairing(theta,E):
    u=frames(theta)[0]
    return u.T@np.kron(sum(e*t for e,t in zip(E,T)),B)@u

def determinant(theta):
    _,v,D,_=frames(theta)
    bar=np.kron(np.eye(16),VP.conj().T)
    return np.linalg.det(bar@D@v)

def ab(theta):
    return (float(np.cos(theta)**2),float(np.cos(1.5*theta)**2),
            float(-np.sin(2*theta)),float(-1.5*np.sin(3*theta)))

def measure(theta):
    a,b,da,db=ab(theta);z=0.;dz=0.
    for k,c in enumerate(COEFFICIENTS):
        n=8-k;c=float(c);z+=c*a**k*b**n
        if k:dz+=c*k*a**(k-1)*b**n*da
        if n:dz+=c*n*a**k*b**(n-1)*db
    return z,dz

def physical(theta):
    return float(np.prod(np.cos(Q*theta/2)**2))

def full(theta):return physical(theta)*measure(theta)[0]

def matrix_check():
    rng=np.random.default_rng(615)
    errors=[];det_errors=[];gaps=[];orientation=[]
    base=np.eye(10)[0]
    pf0=pfaffian(pairing(0.,base));det0=determinant(0.)
    u0,v0,_,_=frames(0.);jac0=np.linalg.det(np.column_stack((u0,v0)))
    assert abs(pf0-1)<3e-13 and abs(det0-1)<3e-13
    for theta in (0.,.17,.31,np.pi/3,np.pi/2,2.14):
        u,v,D,H=frames(theta);g5=np.kron(np.eye(16),spin.G5)
        P=(np.eye(64)-H)/2
        errors.extend([op(u@u.conj().T-P),op(H@H-np.eye(64)),
            op(D@g5+g5@D-2*D@g5@D),op(u.conj().T@v)])
        gaps.append(float(min(abs(np.linalg.eigvalsh(H)))))
        orientation.append(abs(np.linalg.det(np.column_stack((u,v)))/jac0-1))
        det_errors.append(abs(determinant(theta)/det0-physical(theta)))
        for _ in range(5):
            E=rng.normal(size=10);E/=np.linalg.norm(E);r=float(E[:6]@E[:6])
            A=pairing(theta,E);a,b,_,_=ab(theta)
            Te=sum(e*t*(np.cos(theta) if j<6 else np.cos(1.5*theta))
                   for j,(e,t) in enumerate(zip(E,T)))
            expected=-np.kron(Te,EPS)
            value=pfaffian(A)/pf0;analytic=(a*r+b*(1-r))**8
            errors.extend([op(A-expected),abs(value-analytic),abs(pfaffian(A)**2-np.linalg.det(A))])
    assert max(errors+det_errors+orientation)<2e-12 and min(gaps)>1-1e-13
    # An actual singular auxiliary configuration, despite an exact Wilson gap.
    singular=pairing(np.pi/3,np.eye(10)[6])
    singular_norm=op(singular);assert singular_norm<8e-15
    return dict(max_pair_projection_GW_residual=float(max(errors)),
        max_physical_determinant_residual=float(max(det_errors)),
        frame_Jacobian_variation=float(max(orientation)),Wilson_gap_min=min(gaps),
        individual_auxiliary_singular_norm=singular_norm,
        periodic_spatial_one_site_AP_temporal=True,original_internal_dimension=16,
        original_four_Euclidean_directions_and_m0_are_inputs=True)

def source_check():
    assert sum(COEFFICIENTS)==1
    nodes,weights=np.polynomial.legendre.leggauss(6)
    r=(nodes+1)/2;weights=weights/2*12*r**2*(1-r)
    rows=[];errors=[]
    for theta in (.07,.17,.31,np.pi/3,np.pi/2,2.14):
        a,b,da,db=ab(theta);z,dz=measure(theta)
        quad=float(np.dot(weights,(a*r+b*(1-r))**8))
        deriv=float(np.dot(weights,8*(a*r+b*(1-r))**7*(da*r+db*(1-r))))
        errors.extend([abs(z-quad),abs(dz-deriv)])
        row=dict(theta=theta,auxiliary_factor=z,auxiliary_log_derivative=dz/z)
        if theta<.4:
            bare_log=float(-np.sum(Q*np.tan(Q*theta/2)))
            total_log=bare_log+dz/z
            step=1e-6
            fd=(np.log(full(theta+step))-np.log(full(theta-step)))/(2*step)
            derivative_error=abs(fd-total_log)
            assert derivative_error<2e-7
            row.update(physical_determinant=physical(theta),full_factor=full(theta),
                       bare_log_derivative=bare_log,full_log_derivative=total_log,
                       finite_difference_error=derivative_error,
                       omitted_source_if_auxiliary_discarded=abs(dz/z))
        rows.append(row)
    z1,d1=measure(np.pi/3);z2,d2=measure(np.pi/2)
    assert abs(z1-1/720896)<1e-18 and abs(d1/z1+16*np.sqrt(3))<1e-12
    assert abs(z2-1/14080)<1e-18 and abs(d2/z2-24)<1e-12
    assert max(errors)<1e-13
    return dict(exact_coefficients=[str(x) for x in COEFFICIENTS],
        integration_crosscheck_error=max(errors),rows=rows,
        auxiliary_strictly_positive_for_all_real_theta_analytic=True,
        full_factor_nonzero_guaranteed_window='abs(theta)<pi/6',
        individual_zero_not_integrated_zero=True,
        source_convention='log derivative; induced action -log Z has opposite sign')

def condition_check():
    # The subgroup preserves the split R6 + R4, but has no invariant vector.
    theta=.17
    V5=old.carrier(np.eye(3),np.eye(2),np.exp(1j*theta))
    real10=np.block([[V5.real,-V5.imag],[V5.imag,V5.real]])
    smallest=float(np.linalg.svd(real10-np.eye(10),compute_uv=False)[-1])
    assert smallest>.3
    # No hypercharge-singlet in 5 + conjugate(5); same for conjugate Spin9.
    fixed_rank=int(np.linalg.matrix_rank(real10-np.eye(10),tol=1e-12));assert fixed_rank==10
    r0=Fraction(3,5)
    variance=Fraction(3*2,5*5*6);difference=Fraction(28)*variance*Fraction(25,16)
    assert variance==Fraction(1,25) and difference==Fraction(7,4)
    rows=[]
    for theta in (.03,.06,.17,.31):
        a,b,da,db=ab(theta);c=.6*a+.4*b;dc=.6*da+.4*db
        z,dz=measure(theta);alternative=c**8
        assert z>alternative
        rows.append(dict(theta=theta,S9_factor=z,fixed_split_factor=alternative,
            log_source_difference=dz/z-8*dc/c,
            quartic_difference_ratio=(z-alternative)/theta**4))
    assert abs(rows[0]['quartic_difference_ratio']-1.75)<.03
    return dict(original_subgroup_no_fixed_vector_rank=fixed_rank,
        flat_holonomy_smallest_distance_from_fixed_vector=smallest,
        cannot_inherit_Spin9_sufficient_proof_by_group_inclusion=True,
        invariant_measures_same_second_moment=str(r0/6),
        sphere_split_variance=str(variance),fixed_split_variance='0',
        leading_partition_difference_theta4=str(difference),rows=rows,
        alternative_is_measure_choice_not_proved_full_chiral_theory=True,
        uniform_S9_choice_is_fixed_in_main_candidate=True)

def run():
    deps=('research_note_606.md','research_note_610.md','research_note_612.md','research_note_613.md',
          'research_note_614.md','joint_spinor_subgroup_mass.py','joint_chiral_fibre_source.py')
    return dict(round=615,tests_run=3,failures=0,errors=0,matrix=matrix_check(),
        integrated_source=source_check(),conditions=condition_check(),
        dependency_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in deps},
        scope=dict(original_subgroup_and_mature_auxiliary_measure_used=True,
            exact_finite_slice_integral_and_common_source=True,
            one_site_zero_mass_background_no_scalar_integration=True,
            no_full_locality_reconstruction_or_GR_claim=True,
            quantum_measure_selection_remains_declared_input=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--write-results',action='store_true')
    args=p.parse_args();result=run()
    if args.write_results:
        with TARGET.open('x',encoding='utf8',newline='\n') as f:f.write(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    else:assert json.loads(TARGET.read_text('utf8'))==result
    print(json.dumps(result,ensure_ascii=False))
