"""Independent round 1072 verification. Default execution only prints JSON.
No imports from author or historical project code. No output files are written.
Symmetric subspace: average of six factor permutation matrices.
Boundary model: independently assembled from the explicit round 433 definition.
"""
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np

TOL = 4e-12
counts = {}
residuals = {}
def check(ok, group, message):
    counts[group] = counts.get(group, 0) + 1
    if not bool(ok):
        raise AssertionError(group + ': ' + message)
def close(a, b, group, message, tol=TOL):
    error = float(np.max(np.abs(np.asarray(a) - np.asarray(b))))
    residuals[group] = max(error, residuals.get(group, 0.0))
    check(error <= tol, group, message)

I2 = np.eye(2, dtype=complex)
I4 = np.eye(4, dtype=complex)
I8 = np.eye(8, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
pauli = (X, Y, Z)
def kron3(a,b,c):
    return np.kron(np.kron(a,b),c)
def factor_permutation(p):
    matrix = np.zeros((8,8), dtype=complex)
    for bits in itertools.product((0,1), repeat=3):
        i = 4*bits[0]+2*bits[1]+bits[2]
        y = [bits[p[j]] for j in range(3)]
        matrix[4*y[0]+2*y[1]+y[2], i] = 1
    return matrix
permutations = [factor_permutation(p) for p in itertools.permutations(range(3))]
P_sym = sum(permutations)/6
P_half = I8-P_sym
singlet = np.array([0,1,-1,0], dtype=complex)/math.sqrt(2)
P_s = np.outer(singlet,singlet.conj())
P_t = I4-P_s
Q_s = np.kron(P_s,I2)
Q_t = P_half-Q_s
for p,dimension in ((P_sym,4),(P_half,4),(Q_s,2),(Q_t,2)):
    close(p@p,p,'projectors','idempotence')
    close(p.conj().T,p,'projectors','self-adjointness')
    close(np.trace(p),dimension,'projectors','rank by trace')
close(Q_s@Q_t,0,'projectors','two code summands are orthogonal')

# Four common AB Pauli conjugations, C untouched. This is an operator identity,
# not an empirical inference from the input states tested below.
unitaries = (I2,-1j*X,-1j*Y,-1j*Z)
V = [kron3(u,u,I2) for u in unitaries]
mean_effect = sum(v.conj().T@P_sym@v for v in V)/4
close(mean_effect,(2/3)*np.kron(P_t,I2),
      'universal_adjoint_identity','common Pauli average of the leakage effect')
close(P_half@(sum(v.conj().T@P_sym@v for v in V[1:])/3)@P_half,
      (8/9)*P_half@np.kron(P_t,I2)@P_half,
      'three_pi_code_identity','identity operation has zero leakage in the code')

def marginal_AB(rho):
    return np.trace(rho.reshape(4,2,4,2),axis1=1,axis2=3)
def werner(w):
    return (1-w)*I4/4+w*P_s
def attaining_state(w):
    ps=(1+3*w)/4
    return ps*Q_s/2+(1-ps)*Q_t/2

def rotation(axis,theta):
    axis=np.asarray(axis,dtype=float)
    axis=axis/np.linalg.norm(axis)
    return math.cos(theta/2)*I2-1j*math.sin(theta/2)*sum(a*s for a,s in zip(axis,pauli))

axes=[(1,0,0),(0,1,0),(0,0,1),(1,1,1),(1,-2,3),(0,2,-1)]
angles=[0,math.pi/7,math.pi/2,math.pi,3*math.pi/2,2*math.pi]
# Stronger identity on the original code: every input and every rotation axis.
# A z-axis Clebsch-Gordan calculation plus collective covariance proves the
# continuum statement; this full-matrix check validates its explicit formula.
for axis,theta in itertools.product(axes,angles):
    u=rotation(axis,theta)
    v=kron3(u,u,I2)
    close(P_half@v.conj().T@P_sym@v@P_half,
          (8/9)*math.sin(theta/2)**2*Q_t,
          'all_axis_code_operator_identity','every code input has the same triplet-weight leakage factor')
ws=[Fraction(-1,3),Fraction(0),Fraction(1,4),Fraction(1,3),Fraction(1,2),Fraction(1)]
state_rows=[]
for wf in ws:
    w=float(wf)
    rho=attaining_state(w)
    check(np.linalg.eigvalsh(rho).min()>=-TOL,'preparation','positive state')
    close(np.trace(rho),1,'preparation','trace one')
    close(P_half@rho@P_half,rho,'preparation','starts inside the original code')
    close(marginal_AB(rho),werner(w),'preparation','prescribed Werner marginal')
    leaks=[]
    for v in V:
        moved=v@rho@v.conj().T
        close(marginal_AB(moved),werner(w),'Werner_preservation','same AB state after each revision')
        leaks.append(float(np.trace(P_sym@moved).real))
    close(np.mean(leaks),(1-w)/2,'four_average','exact four-revision average')
    close(np.mean(leaks[1:]),2*(1-w)/3,'three_pi_average','exact three-pi average')
    for axis,theta in itertools.product(axes,angles):
        u=rotation(axis,theta)
        v=kron3(u,u,I2)
        moved=v@rho@v.conj().T
        close(marginal_AB(moved),werner(w),'all_axis_AB_preservation','AB invariant under common U')
        observed=float(np.trace(P_sym@moved).real)
        expected=(1-w)*(1-math.cos(theta))/3
        close(observed,expected,'attainment_formula','all-angle analytic attainment formula')
        check(observed<=2*(1-w)/3+TOL,'attainment_range','bounded by the sharp maximum')
    state_rows.append({'w':str(wf),'singlet_weight':str((1+3*wf)/4),
                       'four_leakages':leaks,'four_average':float(np.mean(leaks)),
                       'three_pi_average':float(np.mean(leaks[1:])),
                       'sharp_minimax':str(2*(1-wf)/3)})

# Generic states need not lie in the code or have Werner marginal: they check
# the stronger, full-8-dimensional adjoint identity independently of preparation.
rng=np.random.default_rng(107209)
for _ in range(12):
    a=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8))
    rho=a@a.conj().T
    rho=rho/np.trace(rho)
    observed=sum(np.trace(P_sym@v@rho@v.conj().T).real for v in V)/4
    expected=(2/3)*np.trace(np.kron(P_t,I2)@rho).real
    close(observed,expected,'generic_state_identity','no code or Werner condition needed')

# A purification supplies a correlated passive reference, rather than silently
# assuming the code input is independent of all other internal systems.
rho=attaining_state(.25)
evals,evecs=np.linalg.eigh(rho)
positive=np.flatnonzero(evals>1e-10)
reference_dim=len(positive)
psi=sum(math.sqrt(float(evals[k]))*np.kron(evecs[:,k],np.eye(reference_dim)[:,j])
        for j,k in enumerate(positive))
pure=np.outer(psi,psi.conj())
close(np.trace(pure.reshape(8,reference_dim,8,reference_dim),axis1=1,axis2=3),
      rho,'passive_reference','purification recovers the code state')
ref_before=np.trace(pure.reshape(8,reference_dim,8,reference_dim),axis1=0,axis2=2)
reference_leaks=[]
for v in V:
    full=np.kron(v,np.eye(reference_dim))
    moved=full@pure@full.conj().T
    close(np.trace(moved.reshape(8,reference_dim,8,reference_dim),axis1=0,axis2=2),
          ref_before,'passive_reference','untouched reference marginal')
    close(np.trace(moved.reshape(4,2*reference_dim,4,2*reference_dim),axis1=1,axis2=3),
          werner(.25),'passive_reference','Werner AB with a correlated reference')
    reference_leaks.append(float(np.trace(np.kron(P_sym,np.eye(reference_dim))@moved).real))
close(np.mean(reference_leaks[1:]),.5,'passive_reference','sharp pi average survives reference extension')

# Original round 433, independently assembled: no new routing or control terms.
pairs=[(0,1),(0,2),(1,2)]
number=np.diag([0,1]).astype(complex)
H=np.zeros((64,64),dtype=complex)
edge_effects=[]
for edge,(i,j) in enumerate(pairs):
    p=list(range(3)); p[i],p[j]=p[j],p[i]
    swap=factor_permutation(p)
    factors=[I2,I2,I2]; factors[edge]=X
    xe=kron3(*factors)
    factors=[I2,I2,I2]; factors[edge]=number
    ne=kron3(*factors)
    H+=np.kron(I8,xe+ne)+np.kron(swap,ne)
    edge_effects.append(np.kron(I8,ne))
Ebar=sum(edge_effects)/3
close(H,H.conj().T,'old433_model','fixed Hamiltonian Hermitian')
t=1/8
evals,evecs=np.linalg.eigh(H)
Ut=(evecs*np.exp(-1j*t*evals))@evecs.conj().T
embedding=np.zeros((64,8),dtype=complex)
for j in range(8): embedding[8*j,j]=1
propagated=Ut@embedding
Eeff=propagated.conj().T@Ebar@propagated
coef_sym=float(np.trace(P_sym@Eeff).real/4)
coef_half=float(np.trace(P_half@Eeff).real/4)
close(Eeff,coef_sym*P_sym+coef_half*P_half,
      'old433_effect_decomposition','same fixed mean-edge readout resolves the two sectors')
check(abs(coef_sym-coef_half)>1e-6,'old433_nonzero_contrast','numerical sector contrast')
rho=attaining_state(.25)
edge_probs=[]
for v in V:
    moved=v@rho@v.conj().T
    full=embedding@moved@embedding.conj().T
    direct=float(np.trace(Ebar@Ut@full@Ut.conj().T).real)
    compressed=float(np.trace(Eeff@moved).real)
    leak=float(np.trace(P_sym@moved).real)
    close(direct,compressed,'old433_direct_readout','full 64-dimensional result agrees')
    close(direct,coef_half+(coef_sym-coef_half)*leak,
          'old433_direct_readout','leakage changes the actual fixed edge readout')
    edge_probs.append(direct)
close(np.mean(edge_probs[1:])-edge_probs[0],(coef_sym-coef_half)/2,
      'old433_direct_readout','sharp average has a mean-edge contrast')

result={
    'status':'PASS','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    'checks_total':sum(counts.values()),'checks_by_group':counts,
    'max_identity_residual':max(residuals.values()),'residuals_by_group':residuals,
    'preparation_and_revision_rows':state_rows,
    'correlated_reference':{'dimension':reference_dim,'four_leakages':reference_leaks},
    'old433_readout':{'time':'1/8','dimension':64,'symmetric_sector_coefficient':coef_sym,
                     'half_spin_sector_coefficient':coef_half,'difference':coef_sym-coef_half,
                     'four_revision_probabilities':edge_probs,
                     'three_pi_average_minus_identity':float(np.mean(edge_probs[1:])-edge_probs[0])},
    'scope':[
        'The full matrix adjoint identity underlies the universal state quantifier; samples are consistency checks.',
        'The analytic compressed operator identity gives the all-axis formula for every admissible extension, not just the displayed preparation or a grid maximum.',
        'Zero initial symmetric-sector weight is required when the final weight is called leakage.',
        'The result concerns common AB revisions on this original ABC fixed code; it does not rule out other encodings.',
        'No post-revision recovery, moving code, or physical space has been classified.',
        'The old433 readout is reused with identical rule and blank edges; no new coupling mechanism is claimed.',
        'The readout sign is numerically checked here; a rigorous finite-time sign needs analytic bounds in the main proof.'
    ]}
print(json.dumps(result,ensure_ascii=False,indent=2))
