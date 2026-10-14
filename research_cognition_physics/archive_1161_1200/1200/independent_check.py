"""Independent parity/soft-factor checks; imports only frozen 1199 independent stencils.

No import of 1200 author code, no huge Hilbert matrix, and no claim that the
finite lattice has a graviton S matrix. Default only prints deterministic JSON.
"""
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'1199'/'independent_check.py'
spec=importlib.util.spec_from_file_location('frozen_1199_independent',OLD)
old=importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)


def paired_rows(rows, weights, prime):
    terms=defaultdict(F)
    for row, weight in zip(rows,weights):
        for sign in (-1,1):
            terms[tuple(int(x) for x in sign*row % prime)]+=weight
    for key,value in terms.items():
        inverse=tuple((-x)%prime for x in key)
        assert terms[inverse]==value
    return len(terms)


def lorentz(a,b):
    return a[0]*b[0]-sum(a[i]*b[i] for i in (1,2,3))


def run():
    p=5
    eta,v,rd,ro,td,to=old.templates()
    S=old.orbit(eta)
    V=np.concatenate([old.orbit(a) for a in v])
    R=np.concatenate([old.orbit(a) for a in rd+ro])
    T=np.concatenate([old.orbit(a) for a in td+to])
    assert np.all(S@V.T==0)
    assert np.all(R@V.T==0)
    assert np.all(T@S.T==0)
    rank_s,ker_s=old.kernel_basis(S)
    rank_v,ker_v=old.kernel_basis(V)
    n=S.shape[1]
    k=n-rank_s-rank_v
    assert (n,rank_s,rank_v,k)==(384,63,189,132)

    # C is complex-linear (unitary), not complex conjugation.
    omega=np.exp(2j*np.pi/p)
    X=np.roll(np.eye(p,dtype=complex),1,axis=0)
    Z=np.diag(omega**np.arange(p))
    C=np.zeros((p,p),complex)
    for j in range(p):
        C[-j % p,j]=1
    psi=np.arange(1,p+1)+1j*np.arange(p)
    residuals={
      'C_dagger_C':np.linalg.norm(C.conj().T@C-np.eye(p)),
      'C_squared':np.linalg.norm(C@C-np.eye(p)),
      'X_inverse':np.linalg.norm(C@X@C-X.conj().T),
      'Z_inverse':np.linalg.norm(C@Z@C-Z.conj().T),
      'complex_linear':np.linalg.norm(C@(1j*psi)-1j*(C@psi))}
    # Trace(C X^a Z^b)=omega^{-ab/2}; CSS has c dot d=0 mod p.
    inv2=pow(2,-1,p)
    trace_res=0.
    for a in range(p):
      for b in range(p):
        W=np.linalg.matrix_power(X,a)@np.linalg.matrix_power(Z,b)
        trace_res=max(trace_res,float(abs(np.trace(C@W)-omega**((-a*b*inv2)%p))))
    residuals['twisted_trace_formula']=trace_res
    assert max(residuals.values())<1e-13
    anti=np.linalg.norm((C@(1j*psi)).conj()-1j*(C@psi).conj())
    assert anti>1

    # Full original real Hermitian Weyl pairs, not a toy H.
    inversion={
       'S_penalty':paired_rows(S,[F(-p,2)]*len(S),p),
       'V_penalty':paired_rows(V,[F(-p,2)]*len(V),p),
       'R_curvature':paired_rows(R,[F(-p,2)]*192+[F(-p,4)]*192,p),
       'T_curl':paired_rows(T,[F(-p,2)]*192+[F(-p,8)]*384,p)}
    # Each family can carry its independent real coupling; real identity shifts
    # also commute with C. Coefficients do not mix the coupling labels.

    # Quotient ker(V)/row(S): inversion fixes exactly its zero class since p odd.
    # Equivalently all stabilizer terms have Tr(C W)=1, hence Tr(C P)=1.
    dim=p**k
    fixed_classes=1
    dims=((dim+fixed_classes)//2,(dim-fixed_classes)//2)
    assert dims[0]+dims[1]==dim and dims[0]-dims[1]==1 and dims[1]>0

    # Exact normalized code traces for two original, noncentral sine readouts.
    odd={}
    for label,row,kernel in [('R_diag_first',R[0],ker_s),('T_diag_first',T[0],ker_v)]:
        syndrome=row@kernel.T % p
        assert np.any(syndrome)
        assert all(np.any(a*syndrome % p) for a in range(1,p))
        def tr_power(a):
            return F(int(not np.any(a*syndrome % p)))
        trace_sine_squared=(2-tr_power(2)-tr_power(-2))/4
        assert trace_sine_squared==F(1,2)
        odd[label]={'weyl_order_on_code':5,'normalized_trace_sine_squared':str(trace_sine_squared),
                    'nonzero_syndrome_entries':int(np.count_nonzero(syndrome))}
    O=(X-X.conj().T)/(2j)
    assert np.linalg.norm(C@O@C+O)<1e-13
    # Local matrix check of the nonempty odd map; global nonemptiness follows
    # from the preceding exact full-code trace, not from this five-state check.
    ep=(np.eye(p)+C)/2; em=(np.eye(p)-C)/2
    assert np.linalg.norm(ep@O@ep)<1e-13 and np.linalg.norm(em@O@em)<1e-13
    assert np.linalg.norm(em@O@ep)>.1

    # +--- signature, future-pointing HARD limiting momenta. qbar is the soft
    # direction; these are not asserted to conserve five-leg momentum at q>0.
    hard=[(1,0,0,1),(1,0,0,-1),(1,1,0,0),(1,-1,0,0)]
    signs=(-1,-1,1,1)
    qbar=(1,0,1,0)
    assert all(lorentz(a,a)==0 for a in hard+[qbar])
    assert all(sum(s*a[i] for s,a in zip(signs,hard))==0 for i in range(4))
    # e spatial diag(1,0,-1): trace zero, q transverse, deliberately unnormalized.
    assert qbar[1]==qbar[3]==0
    assert -1-(-1)==0
    terms=[F(s*(a[1]**2-a[3]**2),lorentz(a,qbar)) for s,a in zip(signs,hard)]
    residue=sum(terms,F(0))
    assert terms==[F(1)]*4 and residue==4
    return {
      'schema':'round_1200_independent_parity_and_conditional_soft_factor_v1',
      'status':'PASS_INDEPENDENT_FINITE_PARITY_ALGEBRA',
      'dependency':{'path':'../1199/independent_check.py','sha256':hashlib.sha256(OLD.read_bytes()).hexdigest()},
      'independent_implementation':'Frozen independent Laurent stencils; local complex-linear matrix; exact finite-field quotient; exact rational soft contraction; no 1200 author imports',
      'finite_model':{'prime':p,'qudits':n,'scalar_rank':rank_s,'vector_rank':rank_v,'logical_qudits':k,
         'inversion_paired_distinct_terms':inversion},
      'local_matrix_residual_max':float(max(residuals.values())),
      'twisted_trace_exact_formula':'Tr(C X^a Z^b)=omega^(-a*b/2); CSS orthogonality makes every stabilizer contribution 1',
      'code_parity':{'dimension':str(dim),'C_trace':1,'even_dimension':str(dims[0]),'odd_dimension':str(dims[1])},
      'original_odd_readouts':odd,
      'exact_consequences':['[H,C]=0','[P,C]=0','P_minus exp(-itH) P_plus=0 for all real t',
                            'odd total C-parity multi-time correlators vanish in every C-invariant state'],
      'conditional_soft_witness':{'signature':'+---','hard_momenta':hard,'eta':signs,'q_over_omega':qbar,
        'polarization_spatial_diagonal':[1,0,-1],'individual_residues':[str(x) for x in terms],
        'sum_residue':str(residue),'leading_amplitude_coefficient_per_kappa_M4':str(residue/2),
        'condition':'Requires nonzero hard M4, kappa, unbroken C vacuum, C-odd graviton identification and a qualified soft scattering limit; finite p lattice certifies none of these.'},
      'not_certified':['asymptotic_gravitons_or_Lorentz_IR','vacuum_nonbreaking','finite_probability_or_actual_instrument',
                       'all_L_or_N_phases_excluded','C_even_composite_or_broken_vacuum_excluded','complete_GR_generation'],
    }


if __name__=='__main__':
    print(json.dumps(run(),ensure_ascii=False,indent=2))
