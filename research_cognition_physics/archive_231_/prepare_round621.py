"""Freeze 621 primary review; preserve historical receipts and navigation."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        part=re.escape(key)
        if key.isdecimal():part=r'(?<!\d)'+part+r'(?!\d)'
        pats.append(part)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

def write(name,text):
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as stream:stream.write(text)

mapping={'joint_quantum_exchange_noise':'joint_smeared_gauss_noise',
         '619':'620','620':'621','621':'622',
         '3101':'3104','3104':'3106','1169':'1172','1172':'1175',
         '1930':'1937','1937':'1944'}
verify=replace((HERE/'verify_round620.py').read_text('utf8'),mapping)
verify=verify.replace('Verify same-state joint source noise and actual energy exchange.',
                      'Verify finite-window full-Gauss form-source noise and domains.')
verify=verify.replace('==(3,0,0)','==(2,0,0)').replace('run=3,failures=0','run=2,failures=0')
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==12")
write('verify_round621.py',verify)
publish=replace((HERE/'publish_round620.py').read_text('utf8'),
                mapping|{'## 266.':'## 267.','## 171.':'## 172.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第621轮完成：** [原全Gauss热态的有限时间窗联合来源]({p}research_note_621.md)'
         '原固定图形式来源在同一Gauss热态经有限时间窗涂抹后有联合二阶矩及谱截断尾界；'
         '一般形式前提仍不足以保证瞬时噪声。'
         '两组、十二式通过，最新621／3106，1175份编号科学文件、1944份保护证据。'
         '[核验]({p}research_round_621_checks.json)、[条件账]({p}unified_physics_condition_ledger_621.md)。'
         '限固定图与明示来源族；完整因果响应、连续及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（621后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[622共同有限窗量子来源与二阶因果过程]({p}round622_drafts/STATUS.md)，'
       '核噪声、迟致核及接触项能否出自同一保Gauss过程；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('共同来源、混合涨落与交换约束','全Gauss热态、时间窗与共同来源域')
publish=publish.replace('联合量子来源条件不选择维数','有限窗量子来源不构成维数选择')
publish=publish.replace('共同量子来源的交换约束、混合涨落与框架匹配','原全Gauss热态的有限时间窗联合来源')
assert 'hashlib.sha256' in publish
write('publish_round621.py',publish)
write('round621_drafts/research_note_621_draft.md',(HERE/'research_note_621.md').read_text('utf8'))
review='''621 primary-agent review. No independent agent review.
620 is verified progress. Read 590,591,603,612 before assigning this round.
The new theorem uses the actual 603 full finite-graph Gauss Hamiltonian,
stationary thermal state, all energy moments and already verified
Gauss-preserving form sources. It is not inferred from a small matrix.
G=A^(1/2) B A^(1/2) is only a form at entry. Column sums of bounded B give
the smeared second moment using a_m <= a_n+abs(E_m-E_n). This proof
includes arbitrarily small gaps and degeneracies and does not assume G^2.
Real W^(1,1) windows supply finite M0,M1. Triangular windows have unit
integral, transform sinc^2(omega*tau/2), and M1 <= 2/tau. No instantaneous
limit or graph-refinement uniform constant is asserted.
Spectral cutoffs act only on sources while rho remains the same full Gibbs
state. Weighted Hilbert-Schmidt convergence yields the joint Gram matrix.
Symmetry on the finite eigenvector core gives closability, not essential
self-adjointness or a completed physical measurement instrument.
Tail split is n>R/2 or n<=R/2,m>R; the second region has abs(delta)>R/2.
The triangular-window fourth-power Fourier tail and 603 fifth/sixth energy
moments give explicit convergence. General constants are deliberately not
optimized and are not claimed sharp in the finite numerical factor.
Positive transition measure has finite-band linear growth. Detailed balance
uses the same stationary Gibbs weights, including zero-frequency blocks.
Two separately smeared sources do not automatically define a fully
time-ordered retarded kernel with all contact terms.
The rank-two bounded-form counterexample is explicitly abstract, not the
original Gauss model. Its thermal moments are finite but its raw source
variance has a harmonic lower sum; finite windows suppress that tail.
Original 16-dimensional neutral Fock factor validates Fourier normalization,
direct time integration, Gram positivity, bounds and detailed balance only.
Analytic positive-tail bounds exclude floating point rounding error and
are not advertised as interval arithmetic certification.
No completed continuum UV construction, dynamic geometry, internal clock,
Einstein equation or unified physical theory is claimed. No image checks.
'''
for name in ('research_note_621.md','joint_smeared_gauss_noise.py',
             'joint_smeared_gauss_noise_results.json','unified_physics_condition_ledger_621.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round621_drafts/final_review.txt',review)
print('621 verification/publication scripts and primary review prepared')
