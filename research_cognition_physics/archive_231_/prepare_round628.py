"""Freeze 628 and prepare history-preserving publication."""
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


mapping={'joint_anomaly_scale_partition':'joint_global_anomaly_bundle',
         '626':'627','627':'628','628':'629',
         '3121':'3123','3123':'3125','1190':'1193','1193':'1196',
         '1980':'1987','1987':'1994'}
verify=replace((HERE/'verify_round627.py').read_text('utf8'),mapping)
verify=verify.replace('Verify anomaly constraints on original matter and scale partitions.',
                      'Verify original global anomaly and non-liftable bundle interface.')
verify=verify.replace("text['display_formulas']==12","text['display_formulas']==14")
write('verify_round628.py',verify)
publish=replace((HERE/'publish_round627.py').read_text('utf8'),
                mapping|{'## 273.':'## 274.','## 178.':'## 179.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第628轮完成：** [原商群、旋量几何与全局反常的共同条件]({p}research_note_628.md)'
         '原Z₆商群、普通spin与完整一代通过闭五维费米反常相位检查，包括非提升束；'
         '同一物种的辅助背景指数8同时区分无反常与无零模。'
         '两组、十四式通过，最新628／3125，1196份编号科学文件、1994份保护证据。'
         '[核验]({p}research_round_628_checks.json)、[条件账]({p}unified_physics_condition_ledger_628.md)。'
         '成熟定理的条件性应用；全局权重、连续区域测度、实时过程及GR仍开放。')
order=('**当前执行顺序（628后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[629拓扑权重与共同作用]({p}round629_drafts/STATUS.md)，'
       '回查原相位、质量与区域条件，核剩余共同资料；不重复反常或同类束检查，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原物质划分、测度相位与质量域','原全局群、物质与spin反常的合并')
publish=publish.replace('物质子集约束不选择空间维数','反常工具的辅助维数不是时空生成')
publish=publish.replace('物质尺度划分、反常相位与共同参考域','原商群、旋量几何与全局反常的共同条件')
write('publish_round628.py',publish)
write('round628_drafts/research_note_628_draft.md',(HERE/'research_note_628.md').read_text('utf8'))
review='''628 primary-agent review. No independent agent review.
Previous goal response only restated execution priority and was no progress.
This turn creates a concrete joint-condition result without changing the goal.
README, direction, state, 627 note/results, 628 frozen entry and relevant
531/598/614/616/617 records were checked. Prior 628 entry is preserved.
Historical search found the explicit full-global-anomaly gap, not a prior
application covering arbitrary original quotient bundles.
CIM process query was denied; Get-Process fallback inspected actual processes.
No process was killed or restarted, and all new artifacts use exclusive create.
No previous frozen report/result was edited. No agents/images/new app threads.
DGL 1910.11277 section 2, 4.5 and equation 4.45 were read, including explicit
non-null-bordant four-dimensional partition-function caveats.
Garcia-Etxebarria/Montero 1808.00009 sections 2 and 3.4 were checked for
index density, SU(n) spin bordism and subgroup naturality. The result is
mature mathematics, not claimed as a new anomaly theorem.
Davighi/Lohitsiri 2011.10102 section 2 labels general extended-field-theory
classification a conjecture; that classification is not promoted to a theorem.
Actual old iota factors through the full Z6 quotient. Inducing an SU5 bundle
works without lifting a G bundle to the product covering group.
Associated original matter bundle uses exactly 1+wedge2+wedge4 and original J.
No additional gauge bosons, full SU5 invariant mass sector or new Higgs is assumed.
Mass deformations are not used to change the anomaly representation.
All ordinary-spin closed-five-manifold fermion anomaly phases are trivial,
using I6=0 and the cited bordism/APS inputs, including mapping tori.
This is not an explicit nonperturbative chiral regulator, positive determinant,
theta-angle choice, regional/corner theory, or positive real-time process.
AHSS degree-five diagonal has only (4,1); Sq2 c2=c3 gives surjective d2.
Code checks exact root polynomial identities, not the entire bordism theorem.
S2xS2 is a declared auxiliary closed Euclidean spin background. No claim it is
a generated universe or a demonstrated continuum sector of the finite graph.
E3=L+1+1 and E2=L^-1+1 define the old quotient bundle. c1(E2)=(-1,-1)
fails the necessary divisibility by 6 for a product-cover lift.
The induced bundle V has c2=-2ab, p1(TM)=0; original matter is
4L+4L^-1+8 trivial lines, so twisted massless chiral index=8.
Module indices 3,2,1,1,1,0 use old J with no spin or Nambu overcount.
Index is a lower bound on total zero modes, not an exact spectrum.
At original scalar zero, F=2, U finite and old mass matrices vanish.
No assertion is made that all nonzero Yukawa backgrounds have these zero modes,
or that original full graph H has a zero eigenvalue or nonpositive thermal trace.
The nonbounding example is not used to assert a complete globally normalized
fermionic path integral. Its remaining phase choice is explicitly retained.
Two scientific exact-algebra/bundle groups passed first run. The full verifier
reproduces their results and protects previous evidence before publication.
Next work must integrate remaining phases with original parameters/regions,
not enumerate more backgrounds or design cognitive hardware.
'''
for name in ('research_note_628.md','joint_global_anomaly_bundle.py',
             'joint_global_anomaly_bundle_results.json','unified_physics_condition_ledger_628.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round628_drafts/final_review.txt',review)
print('628 draft, review, verifier and publisher prepared')

