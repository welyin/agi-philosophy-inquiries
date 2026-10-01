"""Freeze 627 and prepare verified history-preserving publication."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    patterns=[]
    for key in sorted(mapping,key=len,reverse=True):
        pattern=re.escape(key)
        if key.isdecimal():pattern=r'(?<!\d)'+pattern+r'(?!\d)'
        patterns.append(pattern)
    return re.sub('|'.join(patterns),lambda m:mapping[m.group()],text)


def write(name,text):
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_equilibrium_memory_matching':'joint_anomaly_scale_partition',
         '625':'626','626':'627','627':'628',
         '3118':'3121','3121':'3123','1187':'1190','1190':'1193',
         '1973':'1980','1980':'1987'}
verify=replace((HERE/'verify_round626.py').read_text('utf8'),mapping)
verify=verify.replace('Verify shared thermal memory and source-only locality obstruction.',
                      'Verify anomaly constraints on original matter and scale partitions.')
verify=verify.replace("text['display_formulas']==14","text['display_formulas']==12")
verify=verify.replace("==(3,0,0)","==(2,0,0)").replace('fresh_tests=dict(run=3,','fresh_tests=dict(run=2,')
write('verify_round627.py',verify)
publish=replace((HERE/'publish_round626.py').read_text('utf8'),
                mapping|{'## 272.':'## 273.','## 177.':'## 178.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第627轮完成：** [物质尺度划分、反常相位与共同参考域]({p}research_note_627.md)'
         '原一代完整模块的无补偿删减受反常联立限制，64子集仅4个通过所列必要检查；'
         '原全热态也不支持夸克质量的统一正点态隙。'
         '消去部门须保其测度相位、适用态域及来源。'
         '两组、十二式通过，最新627／3123，1193份编号科学文件、1987份保护证据。'
         '[核验]({p}research_round_627_checks.json)、[条件账]({p}unified_physics_condition_ledger_627.md)。'
         '未完成全局反常、连续手征过程或GR；不排除有补偿项的有效理论。')
order=('**当前执行顺序（627后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[628完整物质、手征测度与共同区域过程]({p}round628_drafts/STATUS.md)，'
       '核原一代、配置域、态、拼接及真实过程的共同接口；不重复阈值扫描，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('共同参考态与局部作用的记忆限制','原物质划分、测度相位与质量域')
publish=publish.replace('时间来源记忆不选择空间维数','物质子集约束不选择空间维数')
publish=publish.replace('共同参考态对局部有效作用的记忆限制','物质尺度划分、反常相位与共同参考域')
write('publish_round627.py',publish)
write('round627_drafts/research_note_627_draft.md',(HERE/'research_note_627.md').read_text('utf8'))
review='''627 primary-agent review. No independent agent review.
626 was verified, published and its seven new hashes and seven historical
navigation snapshots checked before pursuing this next matter interface.
Historical searches found the full anomaly table, weighted-trace distinction,
seesaw/threshold work and faithful thermal support theorem; none is relabelled
as a new discovery. Existing 627 entry explicitly permits this joint interface.
The new exact classification uses six original one-generation Weyl modules,
each retained zero or one time. No added mirror, compensation phase, or
different low-energy gauge target is hidden in the classification.
The 16 internal Weyl components are not doubled for spin or Nambu counting.
Physical right-handed fields are correctly conjugated in the anomaly table.
Trace normalization matches original Q=6Y, cubic color fundamental=1 and
quadratic fundamental index 2T=1. Integer component traces and all 64 subsets
agree with the independent three-equation proof. Only four subsets pass.
This is necessary anomaly arithmetic, not a full construction, a proof of
Standard Model uniqueness, or all global anomalies of the Z6 quotient.
Usual SU2 sign is checked on the spin covering-group sector only.
Finite CAR gauge covariance is not called a continuum chiral anomaly result.
Preskill's actual paper was checked at sections 6 and 7. It explicitly allows
massive gauge effective theories with compensating terms under stated
conditions. The report does not rule out every anomalous light spectrum,
nor assert a handwritten Wess-Zumino/global sign term suffices for a theory.
Needed heavy-sector variation is identified; complete phase, regulator,
determinant line and global implementation remain open.
Original finite mass data are unchanged. Quark BdG blocks are tested against
analytic masses using full original scalar F and potential along h=t h0.
The full configuration domain has no uniform positive pointwise quark mass
gap. This does not assert the full Hamiltonian has no gap or compact resolvent.
The nonzero thermal weight of low-mass open sets follows by applying the
existing 605 faithfulness theorem to a new explicit invariant region; its
probability is not computed. Positive probability does not rule out finite
error low-energy or background-restricted effective theories.
Original vacuum min-quark < max-lepton excludes only the stated single
threshold sorting, with no Yukawa refit or claim of experimental pole masses.
Optional sterile module satisfies these anomaly tests but still requires the
earlier mass/threshold/seesaw analysis. No rederivation of that work.
Two scientific groups passed first run. Full source and result checks are
performed by the verifier. No images, independent agents or new tasks.
Next entry preserves all unified goals and returns to actual chiral measure,
regions and process mapping rather than more subset or threshold scans.
'''
for name in ('research_note_627.md','joint_anomaly_scale_partition.py',
             'joint_anomaly_scale_partition_results.json','unified_physics_condition_ledger_627.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round627_drafts/final_review.txt',review)
print('627 draft, review, verifier and publisher prepared')
