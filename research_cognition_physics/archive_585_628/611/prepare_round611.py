"""Prepare primary review, reproducibility verification and publication."""
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
mapping={'joint_gapped_link_locality':'joint_original_mass_spinor_bridge',
         '609':'610','610':'611','611':'612','3072':'3075','3075':'3078',
         '1139':'1142','1142':'1145','1856':'1863','1863':'1870'}
verify=replace((HERE/'verify_round610.py').read_text('utf8'),mapping)
write('verify_round611.py',verify)
publish=replace((HERE/'publish_round610.py').read_text('utf8'),mapping|{'## 256.':'## 257.','## 161.':'## 162.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第611轮完成：** [原物种质量、旋量投影与共同来源字典]({p}research_note_611.md)'
         '指定旋量字典使原32模式Dirac／Majorana与BdG质量在静止片保持；'
         '朴素I₄嵌入删掉全部Dirac项，变化投影则共同改变质量块与来源。'
         '三组、十二式通过，最新611／3078，1145份编号科学文件、1870份保护证据。'
         '[核验]({p}research_round_611_checks.json)、[条件账]({p}unified_physics_condition_ledger_611.md)。'
         '非物理极点质量或完整实时手征重建；一般背景、测度及GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（611后，优先于下方历史安排）：** 继续核同一对象的完整过程，认知实现后置。'
       '接[612Euclidean关联与正Hamiltonian共同重建]({p}round612_drafts/STATUS.md)，'
       '对接成熟反射正性与自由传播，保留相互作用和尺度范围，不把静止谱相等视为统一完成，目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('谱隙、链路导数与共同局部系数','原物种质量、旋量及来源共同字典')
publish=publish.replace('局部系数不替代无界电动能传播','静止质量不替代完整实时间过程')
publish=publish.replace('共同谱隙、局部链路扰动与Fock来源系数','原物种质量、旋量投影与共同来源字典')
write('publish_round611.py',publish)
write('round611_drafts/research_note_611_draft.md',(HERE/'research_note_611.md').read_text('utf8'))
review="""611 primary-agent final review; no independent agent review.
Uses original 598 mass_matrices(), charges, diagnostic Y and scalar F dependence; no replacement arbitrary flavour data.
Explicit chirality frames and beta=gamma4 are declared inputs. Identity spin lift kills LR masses; beta lift preserves original normal and pairing coefficients in the rest slice.
Original 32-mode, 64-dimensional BdG spectrum is checked. The 64-dimensional ambient one-particle space is not counted as new physical species.
Complement Majorana pair creation is explicitly removed by the 609 rule; allowed original pair stays unchanged.
Internal gauge dictionary and original full-group mass covariance are checked; local-in-spacetime reconstruction is not inferred from a zero-momentum slice.
Both spatial and selected temporal Wilson holonomy fibres are independently extracted for all original charges. Distinct cosine factors concern compressed mass blocks, not physical pole masses or Lorentz tests.
Link derivative includes moving-frame terms. Scalar sources retain the same factors and original F. The covariant connection in this slice is zero in the allowed frame.
Co-rotating ambient mass restores the old block only by changing the selected model unless the entire process and boundaries are transformed; a global covariant local frame is not supplied.
Luscher's Euclidean psi/bar-psi projectors differ. The analytic beta-conjugation defect excludes only a naive pointwise Hilbert-adjoint dictionary, not reflection-positive reconstruction.
Static compatibility is not full Hamiltonian dynamics, chiral measure, continuum, gravity, or the whole goal.
"""
for name in ('research_note_611.md','joint_original_mass_spinor_bridge.py',
             'joint_original_mass_spinor_bridge_results.json','unified_physics_condition_ledger_611.md'):
    review+=name+' sha256='+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round611_drafts/final_review.txt',review)
print('611 review and publication helpers prepared')
