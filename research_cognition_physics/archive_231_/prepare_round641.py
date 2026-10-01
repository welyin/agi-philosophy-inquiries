"""Freeze actual nonlinear records on the shared spatial Gaussian process."""
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
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

mapping={'joint_spatial_block_reference':'joint_blocked_record_history',
         '639':'640','640':'641','641':'642',
         '3157':'3160','3160':'3162','1229':'1232','1232':'1235',
         '2074':'2081','2081':'2089'}
verify=replace((HERE/'verify_round640.py').read_text('utf8'),mapping)
verify=verify.replace('Verify actual spatial blocking, reference, temporal kernel and geometry.',
                      'Verify original nonlinear records and common spatial geometry response.')
verify=verify.replace("['display_formulas']==14","['display_formulas']==12")
verify=verify.replace('(3,0,0)','(2,0,0)').replace('fresh_tests=dict(run=3,','fresh_tests=dict(run=2,')
start=verify.index('    names=');end=verify.index('    preserved=',start)
verify=verify[:start]+"""    names=('unified_physics_condition_ledger_641.md','round641_drafts/research_note_641_draft.md',
           'round641_drafts/final_review.txt','round642_drafts/STATUS.md',
           'round641_drafts/first_results_before_tail_bound.json')
"""+verify[end:]
write('verify_round641.py',verify)
publish=replace((HERE/'publish_round640.py').read_text('utf8'),
                mapping|{'## 286.':'## 287.','## 191.':'## 192.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第641轮完成：** [同一空间合并过程的非线性记录与几何响应]({p}research_note_641.md)'
         '原二次场的同一时间关联与对易核确定正弦二值读取的两时记录及几何得分；'
         '等时热态替代和删除对易核均给实际记录差异。'
         '两组、十二式通过，最新641／3162，1235份编号科学文件、2089份保护证据。'
         '[核验]({p}research_round_641_checks.json)、[条件账]({p}unified_physics_condition_ledger_641.md)。'
         '字段线性化仍为输入；完整Gauss、连续与GR未完成。')
order=('**当前执行顺序（641后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[642完整非线性物质与Gauss的共同尺度接口]({p}round642_drafts/STATUS.md)，'
       '回查原商群、曲目标和非对角电动能，不继续Gaussian记录精度优化；目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原空间合并的共同状态、时间核与来源','原空间过程的非线性记录与几何响应')
publish=publish.replace('二次空间合并不等于完整连续物理','二次记录接口不替代完整Gauss尺度匹配')
publish=publish.replace('一次真实空间合并的共同热态、时间核与几何来源','同一空间合并过程的非线性记录与几何响应')
write('publish_round641.py',publish)
write('round641_drafts/research_note_641_draft.md',(HERE/'research_note_641.md').read_text('utf8'))
review="""641 primary-agent review; no independent agent review.
640 completed, verified, published and postchecked before work on641.
Original L+/- square-root sine function is retained on the LINEARIZED original
singlet pair average. Gaussian states/processes refer to the declared original
neutral tree-level quadratic branch, not the nonlinear interacting Gauss model.
Both original neutral masses and radial coefficients are used on an actual
N=8 three-dimensional periodic lattice with all lattice momenta present.
Haar alias covariance and independent unblocked Fourier sum agree.
Same coarse thermal Gaussian state is reconstructed at every momentum before
the mean-force dynamical comparison. No separate state fit hides a mismatch.
Weyl sandwich phase sign follows [S0,St]=iD; two Fourier coefficients multiply
without an extra conjugation because each original multiplication is L_r.
Non-Gaussian postmeasurement states are handled by the exact ordered product;
no Gibbs reset or Gaussian replacement instrument is inserted.
Analytic strip1.2 gives uniform coefficient tail and per-event probability
bound; finite FFT alias is negligible under the same strip estimate.
Direct t=0 identity, probability normalization, same first-record marginals,
and K16/K24 comparisons are checked. Floating point is not interval arithmetic.
Common geometric covariance derivatives are evaluated numerically, then
propagated through the exact analytic chain rule; distinct-step full history
differences verify convergence. No analytic derivative interval claimed.
Scores are record statistics, not stress tensors. Scalar normalization needed
for640 physical source cannot be reconstructed from normalized probabilities.
First result file is preserved before adding analytic tail and independent
fine-grid checks; previous frozen evidence remains unchanged.
Two groups passed. No continuous spacetime, full gauge spatial map, autonomous
implementation or GR derived. Next642 returns to nonlinear/Gauss matching.
"""
for name in ('research_note_641.md','joint_blocked_record_history.py',
             'joint_blocked_record_history_results.json','unified_physics_condition_ledger_641.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round641_drafts/final_review.txt',review)
print('641 science frozen; verifier and publisher prepared')
