"""Freeze 622 and prepare append-only navigation after science verification."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def replace(text,mapping):
    patterns=[]
    for key in sorted(mapping,key=len,reverse=True):
        pat=re.escape(key)
        if key.isdecimal():pat=r'(?<!\d)'+pat+r'(?!\d)'
        patterns.append(pat)
    return re.sub('|'.join(patterns),lambda m:mapping[m.group()],text)

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

mapping={'joint_smeared_gauss_noise':'joint_causal_source_functional',
         '620':'621','621':'622','622':'623',
         '3104':'3106','3106':'3109','1172':'1175','1175':'1178',
         '1937':'1944','1944':'1951'}
verify=replace((HERE/'verify_round621.py').read_text('utf8'),mapping)
verify=verify.replace('Verify finite-window full-Gauss form-source noise and domains.',
                     'Verify common quadratic CTP sources and causal coefficient limits.')
verify=verify.replace('==(2,0,0)','==(3,0,0)').replace('run=2,failures=0','run=3,failures=0')
verify=verify.replace("text['display_formulas']==12","text['display_formulas']==14")
write('verify_round622.py',verify)
publish=replace((HERE/'publish_round621.py').read_text('utf8'),
                mapping|{'## 267.':'## 268.','## 172.':'## 173.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第622轮完成：** [同一量子过程的有限窗噪声、迟致核与接触项]({p}research_note_622.md)'
         '原Gauss形式族的规整二阶CTP系数共同收敛；原费米双历史复现噪声、迟致与接触项。'
         '零脉冲噪声仍可有高能虚跃迁响应，排除仅凭噪声签收截断。'
         '三组、十四式通过，最新622／3109，1178份编号科学文件、1951份保护证据。'
         '[核验]({p}research_round_622_checks.json)、[条件账]({p}unified_physics_condition_ledger_622.md)。'
         '全模型过程求导交换、连续及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（622后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[623原完整Hamiltonian的共同算符域与实际参数过程]({p}round623_drafts/STATUS.md)，'
       '核原具体动能、势与来源的图范数和演化，保持同一物理对象；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('全Gauss热态、时间窗与共同来源域','共同二阶过程、时间排序与虚跃迁')
publish=publish.replace('有限窗量子来源不构成维数选择','二阶因果匹配不构成维数选择')
publish=publish.replace('原全Gauss热态的有限时间窗联合来源','同一量子过程的有限窗噪声、迟致核与接触项')
assert 'hashlib.sha256' in publish
write('publish_round622.py',publish)
write('round622_drafts/research_note_622_draft.md',(HERE/'research_note_622.md').read_text('utf8'))
review='''622 primary-agent review. No independent agent review.
Previous goal turn 620-621 is verified progress. Current work read README,
direction,state,621 evidence,622 entry and 591/602/603/586 before proceeding.
No new current-date Python jobs were present at the initial process check.
Mature CTP formula checked against Hu-Verdaguer section4; gravity action
and continuum state are not inferred from use of that formula.
Same full Gauss Gibbs state is used at every cutoff. Only the source
perturbation is spectrally projected; spectator blocks remain in the trace.
Common C2 form family and both derivative bounds are explicit conditions.
Dyson coefficients use plus Hamiltonian source, chi=i theta<[G,G]>, and
the state response has minus chi. Contact sign and one-half factors were
checked by independent actual nonlinear two-branch Fock evolution.
Half-axis ordered kernel has endpoint q(0), so |omega K| is bounded but
need not go to zero. Absolute spectral convergence follows from form
row/column bounds and thermal A,A2 moments, not from naive multiplication
of a singular distribution by theta. Both ordered terms are controlled.
Proof is convergence of derivatives of regulated processes; no uniform
Taylor remainder or interchange with full-process differentiation claimed.
Noise, retardation and the original Hessian contact are from one state.
The exact neutral original factor preserves complex Dirac/Majorana masses
and lapse-source mixed contact. Numerical time and amplitude convergence
are checks, not interval certificates or a full Gauss spectrum computation.
During initial development a diagnostic magnitude assertion expected the
retarded pulse integral >1e-3. Its actual independently converged value
is 0.0007693489923, so the nonzero diagnostic was corrected to >1e-5;
no result, model parameter or mathematical identity was adjusted.
Counterexample varies a norm-one form source across a fixed H and Gibbs
state, not a fixed source across cutoffs. Fourier orthogonality makes its
specified pulse noise exactly zero, while the ordered integral is nonzero.
This excludes a uniform class tail bound based only on those data; it does
not contradict convergence for each fixed source or show original Gauss
sources fail. Nonzero second-order phase is not an all-order claim.
Original operator regularity and actual parameter process remain next.
No cognitive implementation, image checks or GR completion claimed.
'''
for n in ('research_note_622.md','joint_causal_source_functional.py',
          'joint_causal_source_functional_results.json','unified_physics_condition_ledger_622.md'):
    review+=n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest()+'\n'
write('round622_drafts/final_review.txt',review)
print('622 verification/publication scripts and primary review prepared')
