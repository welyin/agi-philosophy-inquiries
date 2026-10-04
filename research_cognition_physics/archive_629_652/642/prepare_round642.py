"""Freeze full quotient-matter tree dictionary and actual source obstruction."""
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


mapping={'joint_blocked_record_history':'joint_tree_matter_source',
         '640':'641','641':'642','642':'643',
         '3160':'3162','3162':'3165','1232':'1235','1235':'1238',
         '2081':'2089','2089':'2097'}
verify=replace((HERE/'verify_round641.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original nonlinear records and common spatial geometry response.',
                      'Verify full quotient tree transport and same-state source obstruction.')
verify=verify.replace("['display_formulas']==12","['display_formulas']==16")
verify=verify.replace('(2,0,0)','(3,0,0)').replace('fresh_tests=dict(run=2,','fresh_tests=dict(run=3,')
start=verify.index('    names=');end=verify.index('    preserved=',start)
verify=verify[:start]+"""    names=('unified_physics_condition_ledger_642.md','round642_drafts/research_note_642_draft.md',
           'round642_drafts/final_review.txt','round643_drafts/STATUS.md',
           'round642_drafts/nonlinear_scale_entry.md')
"""+verify[end:]
write('verify_round642.py',verify)

publish=replace((HERE/'publish_round641.py').read_text('utf8'),
                mapping|{'## 287.':'## 288.','## 192.':'## 193.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第642轮完成：** [完整规范物质的共同表示、空间消元与几何来源]({p}research_note_642.md)'
         '原完整商群、曲目标及CAR共用树表示；合法Gauss态在同一内部CAR消元后相同，'
         '完整能源和几何来源仍不同，限定了状态单独闭合的空间合并。'
         '三组、十六式通过，最新642／3165，1238份编号科学文件、2097份保护证据。'
         '[核验]({p}research_round_642_checks.json)、[条件账]({p}unified_physics_condition_ledger_642.md)。'
         '见证非Gibbs／低能态；连续统一与GR仍未完成。')
order=('**当前执行顺序（642后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[643共同模型与剩余跨分支接口]({p}round643_drafts/STATUS.md)，'
       '回填全条件账，核固定图、连续手征、参考与动态几何的共同资料；不继续树控制器设计，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原空间过程的非线性记录与几何响应','完整物质的空间表示与共同来源')
publish=publish.replace('二次记录接口不替代完整Gauss尺度匹配','完整Gauss表示不等于物质来源可删除')
publish=publish.replace('同一空间合并过程的非线性记录与几何响应','完整规范物质的共同表示、空间消元与几何来源')
write('publish_round642.py',publish)
write('round642_drafts/research_note_642_draft.md',(HERE/'research_note_642.md').read_text('utf8'))

review="""642 primary-agent review; no independent agent review.
641 was completed, published and postchecked before642 entry.
The user's condition-integration priority is preserved; no new goal or design task.
Full original quotient G, nonlinear H5 density and32-CAR matter are used.
Tree gauge fixing is established mathematics, not a new theorem attributed here.
Burbano/Bauer reference is pure SU2 with no fermions; no unproved matter result
is imported from that source. The original-model extensions are derived locally.
The root global Gauss constraint survives. Haar changes are triangular, target
rotations preserve F and density, and finite Fock second quantization is unitary.
The tree dictionary is gamma-independent. Original closed forms and geometric
derivatives are transported together. This is representation, not compression.
Tree electric derivatives act on descendant matter AND loop endpoints. Original
non-diagonal same-source K cross terms are retained. Explicit original derivatives
of the witness coefficients fix momentum signs without guessed fermion generators.
The32-mode filled sea has determinant1: charges sum to zero with both spins.
Holes transform in the dual of an existing original module, not a new species.
The particle-hole wavefunction is Gauss-equivariant at every vertex, smooth in
compact group variables, and multiplied by common invariant compact scalar packets.
It is normal and finite-form-energy, not a delta state, full-H eigenstate or Gibbs
state. It need not be low energy. No low-energy coarse-graining no-go is claimed.
Rooted internal singlets have the same retained bosonic/loop/root-CAR state.
General partial traces of root-singlet states need not be supported in the
retained singlet subspace. The two explicit outputs here ARE singlets; the proof
uses only these inputs, not an invalid assertion of a physical channel for all.
Pointwise norm1 and scalar independence give equal bosonic energy means. Original
onsite h has zero module-diagonal blocks and trace; anomalous terms and hopping
have zero means by fixed global/local particle-number selection respectively.
These facts hold throughout the gamma family, so the exact FULL-H mean gap and
FULL geometric source gap equal the electric differences. No H_electric-only
truncation is confused with full-H dynamics. Root charge and occupancies agree.
Original Yukawa changes module occupation; a frozen label is not a full process.
No exponential Fock matrix is computed. Exact selection rules and representation
identities support the proof;32-mode matrices check original coefficients.
Numerical derivatives, positive K and finite nonlinear tree covariance pass.
The result is limited to a specified internal-CAR erasure and declared states.
Continuous limits, full Gibbs compression, quantum geometry and GR remain open.
Next entry returns to the whole-condition ledger, not repeated tree examples.
"""
for name in ('research_note_642.md','joint_tree_matter_source.py',
             'joint_tree_matter_source_results.json','unified_physics_condition_ledger_642.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round642_drafts/final_review.txt',review)
print('642 science frozen; verifier and publisher prepared')
