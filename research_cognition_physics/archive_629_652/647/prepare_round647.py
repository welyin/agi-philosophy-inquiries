"""Freeze original relational readout and joint location sources."""
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
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_overlap_continuum_sources':'joint_relational_readout_source',
         '645':'646','646':'647','647':'648',
         '3174':'3177','3177':'3179','1247':'1250','1250':'1253',
         '2119':'2127','2127':'2135',
         'cross_branch_entry_audit.md':'relational_quantum_entry.md'}
verify=replace((HERE/'verify_round646.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original overlap continuum window and coincident energy source boundary.',
                      'Verify original material reference, readout and joint location sources.')
verify=verify.replace("==(3,0,0)","==(2,0,0)").replace('fresh_tests=dict(run=3,','fresh_tests=dict(run=2,')
write('verify_round647.py',verify)
publish=replace((HERE/'publish_round646.py').read_text('utf8'),
                mapping|{'## 292.':'## 293.','## 197.':'## 198.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第647轮完成：** [原物质参考、关系读取与几何定位来源]({p}research_note_647.md)'
         '原读口与原物质参考共同定位后，完整来源须保读量、体积和定位三项；'
         '冻结定位有严格纯坐标反例，原初片复算验证三项抵消并连接关系体积。'
         '两组、十六式通过，最新647／3179，1253份编号科学文件、2135份保护证据。'
         '[核验]({p}research_round_647_checks.json)、[全条件账]({p}unified_physics_condition_ledger_647.md)。'
         '限经典关系接口；量子instrument、连续映射及动态量子几何仍开放。')
order=('**当前执行顺序（647后，优先于下方历史安排）：** 先尽量整合已有条件，认知系统设计后置。'
       '接[648关系区域、共同量子来源与实际操作]({p}round648_drafts/STATUS.md)，'
       '优先检查同一过程、同一区域和同一来源的连接；不另开装置或参考优化，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('同一原谱的有效连续窗口与来源边界','原物质参考与关系定位的共同来源')
publish=publish.replace('分离时间连续不替代重合来源匹配','读取、体积与定位必须共同输送')
publish=publish.replace('原传播核的共同连续窗口与重合来源边界','原物质参考、关系读取与几何定位来源')
write('publish_round647.py',publish)
write('round647_drafts/research_note_647_draft.md',(HERE/'research_note_647.md').read_text('utf8'))
review="""647 primary-agent review; no independent agent review.
646 was completed, verified, published and postchecked. Its2127 evidence files
and7 navigation snapshots were preserved. Entry audit distinguishes the prior
progress from the still-open full quantum and dynamic geometry joint mapping.
The classical branch is the ORIGINAL573 constrained Einstein/Yang-Mills/sigma
initial data, not a new independent clock field or a quantized Gauss GR model.
548 rank obstruction and573 local reference existence are inherited, not new.
X=(h,s,A,B) uses Lorentz metric contractions, so both matter and metric variations
are retained. Fixed relational labels include the inverse-map derivative.
Original624 bounded readout e_r(s) becomes a label function at X^1=s. This does
not remove quantum fluctuations or invalidate the original fixed-graph CP map.
The region functional has matter, volume and location derivatives. Pure Lie
variation is a total divergence with transported support/branch or closed slice.
Positive normalized classical averages are not claimed as quantum effects.
At zeros of f the undivided integral variation is used, not a bounded log claim.
The strict frozen-location counterexample fixes its diagnostic vector field at
the baseline solution. An open full-rank patch and monotone e(s) make the integral
strict. Diagnostic chi selects a gauge variation, not a new physical field.
Relational volume variation includes determinant and inverse-map terms. Single
inverse expressions apply only locally;573 global degeneracies remain explicit.
Numerics reuse573 data. The y-shift depends on z, has unit Jacobian and preserves
the torus. Full pulled-back inverse spatial metric includes yz entries. Original
normal velocities transform as scalars; momentum densities have unit Jacobian.
Fourier shift exactly evaluates the collocation polynomial; it is not an interval
certificate for the continuum solution. The finite profiles are diagnostic only.
Numerics test the spatial identity, not full4D evolution or BV quantization.
The generic relational formulas are established tools; the contribution is the
original readout/reference restriction and joint source on the same old fields.
Development: guessed finite wrong-response threshold1e-5 failed (actual7.03e-6).
No result file was written on that failure. A diagnostic optimized-Python run
was used only to inspect amplitudes, not counted as verification. Replaced the
amplitude guess by separation from measured numerical residuals, without changing
data, parameters or physics. Final ordinary Python with all assertions passed.
Two check groups,16 formulas; independent difference steps and three grids agree.
No image checks, new cognition axioms, apparatus design or goal change.
Next addresses joint region/quantum-source interfaces, not profile optimization.
"""
for name in ('research_note_647.md','joint_relational_readout_source.py',
             'joint_relational_readout_source_results.json','unified_physics_condition_ledger_647.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round647_drafts/final_review.txt',review)
print('647 science frozen; verifier and publisher prepared')
