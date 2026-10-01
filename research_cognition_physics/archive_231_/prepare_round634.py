"""Freeze 634: original record charges, noise and source boundary."""
import hashlib
import re
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

mapping={'joint_continuum_record_sources':'joint_record_ward_boundary',
         '632':'633','633':'634','634':'635',
         '3136':'3140','3140':'3143','1208':'1211','1211':'1214',
         '2022':'2029','2029':'2036'}
verify=replace((HERE/'verify_round633.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original continuous record states, sources and causal counterexample.',
                      'Verify original record charge noise and distributional source boundary.')
verify=verify.replace("['display_formulas']==18","['display_formulas']==14")
verify=verify.replace("==(4,0,0)","==(3,0,0)").replace('fresh_tests=dict(run=4,','fresh_tests=dict(run=3,')
write('verify_round634.py',verify)
publish=replace((HERE/'publish_round633.py').read_text('utf8'),
                mapping|{'## 279.':'## 280.','## 184.':'## 185.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第634轮完成：** [记录条件化、守恒来源与几何约束的共同边界]({p}research_note_634.md)'
         '原连续记录的概率共同固定分支能源／动量及交叉噪声；'
         '无补偿地拼接前后来源产生非零表面散度，临时接触项不能抵消净荷。'
         '三组、十四式通过，最新634／3143，1214份编号科学文件、2036份保护证据。'
         '[核验]({p}research_round_634_checks.json)、[条件账]({p}unified_physics_condition_ledger_634.md)。'
         '仅排除指定接法；真实总过程、图映射及GR仍开放。')
order=('**当前执行顺序（634后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[635区域熵、原物种与共同几何有效项]({p}round635_drafts/STATUS.md)，'
       '先回查旧面积熵与边界结果，核同一处方的真实连接；不继续记录硬件设计，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原连续记录与物质几何来源的共同连接','原记录四动量与几何守恒的共同边界')
publish=publish.replace('记录条件的整合不选择空间维数','记录来源匹配不完成空间或GR推导')
publish=publish.replace('连续记录、态正则性与同一物质的几何来源','记录条件化、守恒来源与几何约束的共同边界')
write('publish_round634.py',publish)
write('round634_drafts/research_note_634_draft.md',(HERE/'research_note_634.md').read_text('utf8'))
review="""634 primary-agent review. No independent agent review.
633 was completed, verified and published before proceeding to this unit.
600/601 full-matter constraint matching, 618 boundary conditions, 620 exchange
noise, 624 conditional sources, and 593 earlier resource obstruction reviewed.
General finite-graph energy compensation and Bianchi identities are inherited,
not restated as new scientific discoveries.
FV 1810.06512 section 3.3 equations 3.31/3.32 and theorem 3.4 were checked;
conditioning may change retrodictions and is not a physical stress jump rule.
Hu-Ver daguer 0802.0658 section 3 was checked for conserved stress/noise sources.
The original 633 compact packet is phase modulated, preserving its support
and original sterile Dirac/Majorana mixing. No new field or device is added.
The asymmetric p is evaluated with original BdG projector, not imposed.
Both energy and momentum annihilate the vacuum; conditional charge formulas
therefore have no vacuum-pair cross terms. Momentum holes require k -> -k.
Their positive charge means are independently constrained by the same p.
The between-record covariance is a lower positive contribution to total
quantum covariance, not the entire stress-noise kernel or a classical QFT.
Source splicing is an explicitly tested additional ansatz, not an automatic
consequence of conditioning. Body-wise conservation does not cancel its delta.
The nonzero integrated energy AND momentum jump rules out compact-time
compensators with no asymptotic flux. Persistent resource changes are allowed.
The free fields live on a fixed flat background. Geometry is tested only
against the declared linearized Einstein/Bianchi necessary condition.
No claim of full nonlinear metric no-go, full GR derivation or model failure.
Changing scalar or geometry dynamics requires the other sources from 600/620.
Charge compensation per outcome compares the SAME conditional initial and
final data. It must not substitute the unconditional initial apparatus mean.
The joint +/- covariance block additionally assumes a sharp common total
initial charge and charge-compatible terminal records. These are explicit.
No compensating stress field, autonomous apparatus, metric solution, or
physical finite-time implementation is claimed from charge bookkeeping.
Quadrature normalization adjustment is explicit. The inherited compact-packet
tail estimate is modified for momentum shift; grid agreement is not rigorous
total error, and the tail coefficient is not interval-certified.
Three groups passed first execution; 14 display formulas in the note.
Further apparatus work is deferred per user. Next 635 turns to unintegrated
region entropy/area and the original matter/geometric action, after history review.
All older evidence and frozen STATUS files preserved; broad goal unchanged.
"""
for name in ('research_note_634.md','joint_record_ward_boundary.py',
             'joint_record_ward_boundary_results.json','unified_physics_condition_ledger_634.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round634_drafts/final_review.txt',review)
print('634 science frozen; verifier and publisher prepared')

