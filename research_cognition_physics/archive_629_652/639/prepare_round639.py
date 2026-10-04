"""Freeze 639: original regional read, physical support and shared geometry."""
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

mapping={'joint_record_relative_entropy':'joint_regional_read_scale',
         '637':'638','638':'639','639':'640',
         '3151':'3154','3154':'3157','1223':'1226','1226':'1229',
         '2059':'2067','2067':'2074'}
verify=replace((HERE/'verify_round638.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original record reference support and common relative entropy limits.',
                      'Verify original regional read, fixed support and geometry variation.')
verify=verify.replace("['display_formulas']==16","['display_formulas']==15")
start=verify.index('    names=');end=verify.index('    preserved=',start)
verify=verify[:start]+"""    names=('unified_physics_condition_ledger_639.md','round639_drafts/research_note_639_draft.md',
           'round639_drafts/final_review.txt','round640_drafts/STATUS.md')
"""+verify[end:]
write('verify_round639.py',verify)

publish=replace((HERE/'publish_round638.py').read_text('utf8'),
                mapping|{'## 284.':'## 285.','## 189.':'## 190.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第639轮完成：** [原区域读量、物理体积与尺度统一的相对信息预算]({p}research_note_639.md)'
         '原完整H的区域注能范数由有效物理体积给出，固定支撑可保统一信息上界；'
         '非均匀几何变分还须保读取权重的变化。'
         '三组、十五式通过，最新639／3157，1229份编号科学文件、2074份保护证据。'
         '[核验]({p}research_round_639_checks.json)、[条件账]({p}unified_physics_condition_ledger_639.md)。'
         '集体读取为候选操作；因果实现、跨图态极限及GR仍开放。')
order=('**当前执行顺序（639后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[640实际尺度映射、共同参考与联合观测]({p}round640_drafts/STATUS.md)，'
       '核原群、参考、记录及物质／几何来源能否共同输送；不继续优化剖面或装置，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原记录、参考和区域相对熵的共同近似','原区域读量与物理体积的统一预算')
publish=publish.replace('谱相对熵极限不替代空间连续匹配','区域信息上界不替代跨图连续状态')
publish=publish.replace('原记录、共同参考与区域相对熵的联合近似','原区域读量、物理体积与尺度统一的相对信息预算')
write('publish_round639.py',publish)
write('round639_drafts/research_note_639_draft.md',(HERE/'research_note_639.md').read_text('utf8'))

review="""639 primary-agent review. No independent agent review.
Read current README, direction, state, 638 results and 639 candidate note/code.
Original full nonlinear scalar target, Gauss, gauge and CAR interactions retained.
Regional weighted singlet instrument is a new conditional operation family:
no instantaneous causal implementation or cognition axiom is inferred.
Dirichlet product identity retains the original inverse target metric.
Exact norm uses continuous multiplication supremum at the target origin;
Gauss-invariant finite-width packets approach it. No uniformly bounded energy
sequence or Gibbs saturation is asserted.
Fixed positive effective physical volume and bounded beta give a uniform
information UPPER bound only; no lower signal bound or cross-graph state limit.
The original 3D periodic graph, geometry and chosen profile are explicit inputs.
The numerical support test uses cos(x)>0, without an artificial positive cutoff.
Smooth fixed support plus explicit node count establishes a conservative bound
for all N divisible by 8; four sampled configurations do not replace that proof.
The full Gibbs conclusions are analytic. Configuration injections are not
thermal expectations or calculations of actual Gibbs relative entropy.
Geometry changes weights and the instrument: total injection derivative
includes both original fixed-instrument source and operation variation.
Uniform conformal direction recovers the previous round593 identity.
The scalar geometry-energy shift leaves normalized Gibbs information unchanged;
this inherited limitation is not counted as a separate new theorem.
Next640 concerns one actual spatial scale map for reference, regions, record
and matter/geometry sources. Instrument engineering remains deferred.
Three numerical groups, fifteen formulas and links are checked by verifier.
No continuum spacetime, area law, Newton value or GR derivation is claimed.
"""
for name in ('research_note_639.md','joint_regional_read_scale.py',
             'joint_regional_read_scale_results.json','unified_physics_condition_ledger_639.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round639_drafts/final_review.txt',review)
print('639 science frozen; verifier and publisher prepared')
