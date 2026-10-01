"""Prepare review and publication helpers; never modify earlier science files."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
def replace(text,mapping):
    pattern="|".join(re.escape(x) for x in sorted(mapping,key=len,reverse=True))
    return re.sub(pattern,lambda m:mapping[m.group()],text)
def write(name,text):
    path=HERE/name
    path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)
mapping={'joint_projected_process_completion':'joint_projection_local_composition',
    '606':'607','607':'608','608':'609','3063':'3066','3066':'3069',
    '1130':'1133','1133':'1136','1835':'1842','1842':'1849'}
verify=replace((HERE/'verify_round607.py').read_text('utf8'),mapping)
verify=verify.replace("display_formulas']==14","display_formulas']==12")
write('verify_round608.py',verify)
publish=replace((HERE/'publish_round607.py').read_text('utf8'),mapping|{'253':'254','158':'159'})
start=publish.index('summary=')
end=publish.index('planned={}',start)
header="""summary=('**第608轮完成：** [投影补全的独立组合、局域性与共同过程]({p}research_note_608.md)'
         '整体P／Q规则一般不保独立组合，全空间可出现距离无关影响；'
         '实际在位约束的逐项补全保局部支持，并保原P内全部历史、参考和来源。'
         '三组、十二式通过，最新608／3069，1136份编号科学文件、1849份保护证据。'
         '[核验]({p}research_round_608_checks.json)、[条件账]({p}unified_physics_condition_ledger_608.md)。'
         '反例操作出P；实际GW、无界传播及完整物理连接仍开放，无新增独立代理审查。')
order=('**当前执行顺序（608后，优先于下方历史安排）：** 继续由共同过程与组合约束合并条件，认知实现后置。'
       '接[609一粒子协变补全与Fock共同来源]({p}round609_drafts/STATUS.md)，'
       '返回真实投影及原质量接口，不把在位乘积分解输入GW核，目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('投影补全的历史、Gauss与来源共同闭合','独立组合检验约束全局投影补全')
publish=publish.replace('修改动力学后的交织不替代局域性','共同过程与局域支持的条件性整合')
publish=publish.replace('保完整历史与共同来源的投影动力学补全','投影补全的独立组合、局域性与共同过程')
write('publish_round608.py',publish)
note=(HERE/'research_note_608.md').read_text('utf8')
write('round608_drafts/research_note_608_draft.md',note)
review="""608 primary-agent final review; no independent agent review.
Checked independent tensor-sum defect analytically, including complement-sector off-diagonal terms.
The finite-time receiver counterexample is exact for any N. Its remote encoding leaves P; it is not a physical-sector signalling proof.
Onsite-factor assumption is explicit and is not attributed to GW projectors or Gauss physical regions.
Local conditional expectations preserve support and contract each bounded interaction norm. No unbounded gauge LR theorem is claimed.
Both completions reduce the common P sector and have identical restriction; each specified Kraus map preserves it. Full history and arbitrary reference are preserved there only.
Geometry sources require parameter-independent projectors. Instrument injection is checked as an actual commuting square, including nonzero injection.
Ambient Q dynamics and Gibbs states change; numerical checks explicitly detect this. The inherited P-sector dynamics are not claimed newly derived from cognition.
The next step returns to the actual one-particle/Fock projector rather than replacing it with the qutrit diagnostic.
"""
for name in ('research_note_608.md','joint_projection_local_composition.py',
             'joint_projection_local_composition_results.json','unified_physics_condition_ledger_608.md'):
    review+=name+' sha256='+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round608_drafts/final_review.txt',review)
print('608 review and publication helpers prepared')
