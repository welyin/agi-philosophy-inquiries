"""Freeze 636: original quotient joint labels, boundary entropy and sources."""
import hashlib,re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)
def write(name,text):
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)
mapping={'joint_entropy_geometric_matching':'joint_quotient_boundary_entropy',
         '634':'635','635':'636','636':'637',
         '3143':'3146','3146':'3149','1214':'1217','1217':'1220',
         '2036':'2043','2043':'2050'}
verify=replace((HERE/'verify_round635.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original mass spectra, replica UV entropy and target-gradient term.',
                      'Verify original quotient labels, regional entropy and electric-source witness.')
verify=verify.replace("['display_formulas']==12","['display_formulas']==14")
write('verify_round636.py',verify)
publish=replace((HERE/'publish_round635.py').read_text('utf8'),
                mapping|{'## 281.':'## 282.','## 186.':'## 187.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第636轮完成：** [原商群的联合边界记录、区域熵与电几何来源]({p}research_note_636.md)'
         '原Z₆边界标签独立化产生13/18禁标签权；'
         '原非对角电来源同均值同协方差的两个合法态仍有不同区域熵。'
         '三组、十四式通过，最新636／3149，1220份编号科学文件、2050份保护证据。'
         '[核验]({p}research_round_636_checks.json)、[条件账]({p}unified_physics_condition_ledger_636.md)。'
         '限原固定图的明示态及电来源；全热态熵、连续匹配和GR仍开放。')
order=('**当前执行顺序（636后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[637原完整热态、区域能源与熵]({p}round637_drafts/STATUS.md)，'
       '复用603完整H与617切分，核共同态的区域熵条件；目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原物种和曲目标的共同几何熵','原商群边界记录与电来源的共同限制')
publish=publish.replace('几何熵匹配不选择维数或Newton常数','边界熵约束不等于连续几何已涌现')
publish=publish.replace('原物种、曲目标梯度与区域熵的共同几何系数','原商群的联合边界记录、区域熵与电几何来源')
write('publish_round636.py',publish)
write('round636_drafts/research_note_636_draft.md',(HERE/'research_note_636.md').read_text('utf8'))
review="""636 primary-agent review. No independent agent review.
Previous goal turn only confirmed priority; classified no progress for research.
Current turn adds actual original-object science, not another status restatement.
README, direction, state, 635 report/results and 636 entry reread.
CIM process inspection was denied; Get-Process fallback returned no Python rows.
No process was killed and no protected history or frozen entry was changed.
Relevant 360/362/365, 589/598/603 and 617 history was checked.
Donnelly 1109.0036 section II and Eq32 supply inherited decomposition.
Donnelly-Wall Maxwell result is not extended to interacting non-Abelian fields.
Original quotient group and integer charge convention retained.
Original six modules checked through existing representation code and actual
four-link character covariance, not just the congruence formula.
Loop states multiplied by gauge-invariant smooth compact scalar packets and
the original CAR vacuum are legitimate states in the full Gauss space.
No claim that vacuum or selected loops form an invariant sector of full H.
All original interactions retained; no pure-electric Gibbs substitution.
Region is a vertex bipartition with two cut edges; scalar packets and Fock
vacuum are pure regional factors. The two labels are correlated, not duplicated.
Extended ordinary trace vs invariant algebra unit-block trace is explicit.
The representation entropy is not called independently accessible bits.
Non-diagonal K is background geometry, not the dynamical scalar F.
Gauss at the bottom-left vertex gives negative cross-Casimir; all other active
outgoing edges counted once. All other graph edges are constant wavefunctions.
Analytic geometric derivative identity holds for arbitrary fixed parameters.
Two positive rational distributions have equal first and second electric
Casimir moments, different third moment and different genuine reduced entropy.
Full stress, magnetic contributions and subsequent evolution may differ.
Initial assertion demanded center entropy difference >.01; actual .00489585
was retained and requirement corrected to >.001. No physics parameters changed.
Three groups now pass; source derivative error 5.38e-11, exact mean/cov
differences zero, entropy difference .417151515830.
No continuum area coefficient, physical cutoff, dimension or GR derived.
Next 637 targets full original thermal region entropy and required energy
control; cognition design remains postponed. Goal is unchanged.
"""
for name in ('research_note_636.md','joint_quotient_boundary_entropy.py',
             'joint_quotient_boundary_entropy_results.json','unified_physics_condition_ledger_636.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round636_drafts/final_review.txt',review)
print('636 science frozen; verifier and publisher prepared')
