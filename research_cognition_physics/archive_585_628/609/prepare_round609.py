"""Prepare 609 verification/publication without changing prior evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def replace(text,mapping):
    return re.sub("|".join(re.escape(x) for x in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)
def write(name,text):
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)
mapping={'joint_projection_local_composition':'joint_fock_covariant_completion',
         '607':'608','608':'609','609':'610','3066':'3069','3069':'3072',
         '1133':'1136','1136':'1139','1842':'1849','1849':'1856'}
verify=replace((HERE/'verify_round608.py').read_text('utf8'),mapping)
verify=verify.replace("display_formulas']==12","display_formulas']==14")
write('verify_round609.py',verify)
publish=replace((HERE/'publish_round608.py').read_text('utf8'),mapping|{'254':'255','159':'160'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第609轮完成：** [一粒子协变补全、Fock组合与共同质量来源]({p}research_note_609.md)'
         '同一投影的dΓ联络及正项保原允许形式并兼容独立组合；'
         '质量须同步保正常pp／qq及仅pp配对，补模式成对产生有明确泄漏反例。'
         '三组、十四式通过，最新609／3072，1139份编号科学文件、1856份保护证据。'
         '[核验]({p}research_round_609_checks.json)、[条件账]({p}unified_physics_condition_ledger_609.md)。'
         '固定图Gauss／热态结论有条件；原物种嵌入、传播、测度及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（609后，优先于下方历史安排）：** 继续用共同条件连接各部门，认知实现后置。'
       '接[610真实谱投影、链路导数与共同配置域]({p}round610_drafts/STATUS.md)，'
       '先核一致谱隙的局部系数界及热态适用域，不把系数衰减当无界传播证明，目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('独立组合检验约束全局投影补全','一粒子补全、Fock组合与质量共同约束')
publish=publish.replace('共同过程与局域支持的条件性整合','一粒子组合相容仍须实际空间传播')
publish=publish.replace('投影补全的独立组合、局域性与共同过程','一粒子协变补全、Fock组合与共同质量来源')
write('publish_round609.py',publish)
write('round609_drafts/research_note_609_draft.md',(HERE/'research_note_609.md').read_text('utf8'))
review="""609 primary-agent final review; no independent agent review.
Checked p'=off-diagonal, C=[p,p'], covariant invariance, exterior lift and positive dGamma Gram term.
Allowed-sector form equals the original differential form. The full extension differs from 607 and is not claimed uniquely physical.
Complement particle-number averaging of an even quadratic mass retains normal pp+qq and pair pp only. Retaining qq pair creation is explicitly tested and fails.
The actual 598 32-mode mass has not been embedded into the 606 Euclidean spinor projector; four-mode mass diagnostic is not advertised as that missing embedding.
Independent one-particle direct sums give even Fock tensor sums, removing the 608 remote-sector condition for independent systems, not proving interacting region factorization.
Closedness and heat comparison require bounded projection derivatives, positive metric and fixed graph. No volume-uniform spectral or unbounded-gauge LR assertion.
Gauge covariance is conditional on the original gauge representation and domain. It is not continuum chiral anomaly cancellation.
The scalar instrument identity uses unchanged scalar kinetic directions and a link-only projector. Geometry differentiation requires gamma-independent p.
The grid discretizes the new covariant form; it is not exact pinching of the old finite-difference grid. Sources differentiate kinetic/geometric terms, not the gamma-independent mass.
The actual Wilson zero-momentum projector is independently extracted from the 64-dimensional kernel. Full histories and arbitrary references have an analytic intertwining proof.
"""
for name in ('research_note_609.md','joint_fock_covariant_completion.py',
             'joint_fock_covariant_completion_results.json','unified_physics_condition_ledger_609.md'):
    review+=name+' sha256='+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round609_drafts/final_review.txt',review)
print('609 review and publication helpers prepared')
