"""855: locality defect of total-word compression and tensor closure repair.

The source coefficients are dimensionless diagnostic pulse parameters,
not the physical normalized sources of the continuum theory.
"""
from pathlib import Path
import argparse,json,math
from fractions import Fraction
import numpy as np
HERE=Path(__file__).resolve().parent;TARGET=HERE/'locality_preserving_task_compression_results.json'

def ladder(n):
    q=np.zeros((n,n))
    for j in range(n-1):q[j,j+1]=q[j+1,j]=j+1
    return q

def expi(h,t):
    e,v=np.linalg.eigh(h);return (v*np.exp(1j*t*e))@v.conj().T

def flip_coeff(q,n):
    # P(B record flips)= (1-<0|cos(2 b q)|0>)/2.
    return -Fraction((-1)**n*2**(2*n),2*math.factorial(2*n))*int(np.linalg.matrix_power(q.astype(np.int64).astype(object),2*n)[0,0])

def run():
    rows=[]
    for d in (1,2,3):
        n=d+2;q=ladder(n);eye=np.eye(n)
        full_a=np.kron(q,eye);full_b=np.kron(eye,q)
        tri=[i*n+j for i in range(n) for j in range(n) if i+j<=d]
        rect=[i*n+j for i in range(d+1) for j in range(d+1)]
        def comp(ids):
            w=np.eye(n*n)[:,ids]
            return w,w.T@full_a@w,w.T@full_b@w
        wt,a,b=comp(tri);wr,ar,br=comp(rect)
        comm=a@b-b@a;commr=ar@br-br@ar
        assert np.max(abs(comm))>0 and np.max(abs(commr))==0
        # All record matrix units act on a separate two-qubit factor.
        # The common projection therefore commutes with the full finite
        # record algebra; this defect is not a failure to close that algebra.
        vac=np.zeros(n*n);vac[0]=1
        # Check every ordered local-generator word up to the declared d.
        import itertools
        worst_tri=worst_rect=0.0
        for k in range(d+1):
            for word in itertools.product((0,1),repeat=k):
                vf=vac.copy();vt=wt.T@vac;vr=wr.T@vac
                for z in word:
                    vf=(full_a,full_b)[z]@vf
                    vt=(a,b)[z]@vt;vr=(ar,br)[z]@vr
                worst_tri=max(worst_tri,float(np.max(abs(vf-wt@vt))))
                worst_rect=max(worst_rect,float(np.max(abs(vf-wr@vr))))
        assert worst_tri==0 and worst_rect==0
        # The exact commutator is a cutoff-boundary effect: it annihilates
        # every environment basis state of total occupation below d.
        interior=[k for k,j in enumerate(tri) if j//n+j%n<d]
        assert np.max(abs(comm[:,interior]))==0
        # Analytic non-signalling violation coefficient a^2 b^(2d):
        # a creates one A excitation; B then has only d-1 levels available.
        coefficient=flip_coeff(ladder(d),d)-flip_coeff(ladder(d+1),d)
        formula=Fraction((-1)**d*2**(2*d-1)*math.factorial(d)**2,math.factorial(2*d))
        assert coefficient==formula and coefficient!=0
        # Actual record probabilities from the finite joint unitary.
        # Local interaction is X_record tensor q_environment. Its record
        # value equals environment occupation parity when starting at zero.
        # The resulting joint sector is unitarily equivalent to a,b above.
        aa=.43;bb=.37
        v0=wt.T@vac;vr0=wr.T@vac
        vt=expi(b,bb)@expi(a,aa)@v0
        vt0=expi(b,bb)@v0
        vr=expi(br,bb)@expi(ar,aa)@vr0
        vrbase=expi(br,bb)@vr0
        mask=np.array([j%n%2==1 for j in tri]);maskr=np.array([j%n%2==1 for j in rect])
        pt=float(sum(abs(vt[mask])**2));pt0=float(sum(abs(vt0[mask])**2))
        pr=float(sum(abs(vr[maskr])**2));pr0=float(sum(abs(vrbase[maskr])**2))
        assert abs(pt-pt0)>1e-7 and abs(pr-pr0)<2e-14
        # Independent joint record-sector calculation with X_A and X_B.
        X=np.array([[0,1],[1,0]],float);I=np.eye(2)
        A=np.kron(np.kron(X,I),a);B=np.kron(np.kron(I,X),b)
        v=np.zeros(4*len(tri));v[0]=1
        vv=expi(B,bb)@expi(A,aa)@v
        PB=np.kron(np.kron(I,np.diag([0.,1.])),np.eye(len(tri)))
        direct=float(np.vdot(vv,PB@vv).real)
        assert abs(direct-pt)<2e-14
        rows.append(dict(declared_word_order=d,triangular_environment_dimension=len(tri),
            rectangular_environment_dimension=len(rect),full_record_factor_dimension=4,
            word_amplitude_error_triangular=worst_tri,word_amplitude_error_rectangular=worst_rect,
            commutator_max_entry_triangular=float(np.max(abs(comm))),
            commutator_max_entry_rectangular=float(np.max(abs(commr))),
            commutator_vanishes_on_all_lower_occupation_states=True,
            mixed_probability_monomial=dict(a_power=2,b_power=2*d,total_degree=2*d+2),
            mixed_probability_coefficient_exact=str(coefficient),
            mixed_probability_coefficient_formula_exact=str(formula),
            pulse_parameters=[aa,bb],triangular_B_flip_with_A=pt,triangular_B_flip_without_A=pt0,
            spurious_influence_triangular=pt-pt0,spurious_influence_rectangular=pr-pr0,
            independent_two_record_calculation_error=abs(direct-pt),
            finite_unitary_state_normalization_error=float(abs(np.vdot(vt,vt)-1))))
    return dict(round=855,all_checks_passed=True,fresh_test_groups=1,rows=rows,
        scope='Two initially independent commuting ladder/record regions. Total-word compression preserves selected finite jets but has a boundary locality defect; product closure preserves exact locality in this declared tensor model.',
        original_598_graph_or_continuum_local_tensor_factorization_proved=False,
        all_compressions_or_all_graph_theories_disproved=False,
        new_inputs_for_conditional_repair=['finite region tensor factors','finite local factor decomposition of required coefficients','finite tensor-rank initial input embedding'],
        full_goal_completed=False)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args();r=run()
    if args.write:
        assert not TARGET.exists();TARGET.write_text(json.dumps(r,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    else:assert r==json.loads(TARGET.read_text('utf-8'))
    print(json.dumps(r,ensure_ascii=False,indent=2))
