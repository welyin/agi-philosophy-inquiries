"""Prepare current review and publication; preserve earlier science."""
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
mapping={'joint_fock_covariant_completion':'joint_gapped_link_locality',
         '608':'609','609':'610','610':'611','3069':'3072','3072':'3075',
         '1136':'1139','1139':'1142','1849':'1856','1856':'1863'}
verify=replace((HERE/'verify_round609.py').read_text('utf8'),mapping)
verify=verify.replace("display_formulas']==14","display_formulas']==12")
write('verify_round610.py',verify)
publish=replace((HERE/'publish_round609.py').read_text('utf8'),mapping|{'255':'256','160':'161'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第610轮完成：** [共同谱隙、局部链路扰动与Fock来源系数]({p}research_note_610.md)'
         '同一谱隙与局部有限秩导数同时控制投影、联络、正能源及Fock尾部；'
         '并核远端改场依赖，界不乘总费米占据数。'
         '三组、十二式通过，最新610／3075，1142份编号科学文件、1863份保护证据。'
         '[核验]({p}research_round_610_checks.json)、[条件账]({p}unified_physics_condition_ledger_610.md)。'
         '硬域非原全热态，系数衰减非无界传播；真实质量、测度及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（610后，优先于下方历史安排）：** 继续整合同一物理对象，认知实现后置。'
       '接[611原物种质量与共同投影的真实字典]({p}round611_drafts/STATUS.md)，'
       '核原32模式与候选旋量投影的质量、规范和来源条件，不继续单独优化衰减率，目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('一粒子补全、Fock组合与质量共同约束','谱隙、链路导数与共同局部系数')
publish=publish.replace('一粒子组合相容仍须实际空间传播','局部系数不替代无界电动能传播')
publish=publish.replace('一粒子协变补全、Fock组合与共同质量来源','共同谱隙、局部链路扰动与Fock来源系数')
write('publish_round610.py',publish)
write('round610_drafts/research_note_610_draft.md',(HERE/'research_note_610.md').read_text('utf8'))
review="""610 primary-agent final review; no independent agent review.
Checked common finite-range, row-block norm, spectral gap and local derivative trace-norm assumptions.
Positive spectral rectangle length is 2(M+g), distance at least g/2; weighted perturbation is at most g/4 for both weight signs.
Double-sided weighted derivative trace norm uses W D=D W=D, not a volume-dependent sum of matrix elements.
Connection and positive-term weighted trace bounds give dGamma norms on all finite Fock occupancies via singular-value/CAR decomposition.
Spatial-index truncation alone is not gauge configuration locality. Two gapped configurations with remote kernel change are compared explicitly, including the p-change term in the connection.
No interpolation gap or arbitrary far-field freezing is assumed. No strict local cylinder operator on the whole gauge configuration space is claimed.
Whole interaction summability requires a growth/coefficient contract. Coefficient bounds do not establish LR for unbounded electric derivatives.
The actual Wilson perturbation is one link; reconstruction of the old full-holonomy kernel and gauge covariance are checked.
The fixed four Euclidean dimensions are input. The old Hamiltonian spatial-domain dictionary remains unproved.
605 full-support and finite-penalty results are inherited without rerunning them, not strengthened into failure of every rough configuration.
Constants are deliberately loose uniform sufficient bounds; small matrices are not used to fit an asymptotic decay theorem.
"""
for name in ('research_note_610.md','joint_gapped_link_locality.py',
             'joint_gapped_link_locality_results.json','unified_physics_condition_ledger_610.md'):
    review+=name+' sha256='+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round610_drafts/final_review.txt',review)
print('610 review and publication helpers prepared')
