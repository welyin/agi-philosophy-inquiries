"""Freeze657 and prepare verification/publication without overwriting history."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    parts=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        parts.append(p)
    return re.sub('|'.join(parts),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as stream:stream.write(text)


mapping={'joint_spatial_auxiliary_geometry':'joint_auxiliary_reflection_gluing',
    '655':'656','656':'657','657':'658','3206':'3210','3210':'3214',
    '1277':'1280','1280':'1283','2218':'2229','2229':'2239',
    'actual_spatial_measure_entry.md':'positive_gluing_entry.md',
    'spatial_auxiliary_probe':'reflection_kernel_probe'}
verify=replace((HERE/'verify_round656.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original spatial auxiliary measure, shared propagation and response.',
                      'Verify original auxiliary reflection Gram identity and S9 integration.')
verify=verify.replace("'round657_drafts/reflection_kernel_probe_results.json',\n           'round657_drafts/initial_planar_results.json'",
                      "'round657_drafts/reflection_kernel_probe_results.json'")
verify=verify.replace("text['display_formulas']==20","text['display_formulas']==21")
write('verify_round657.py',verify)
publish=replace((HERE/'publish_round656.py').read_text('utf8'),mapping|{'## 302.':'## 303.','## 207.':'## 208.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第657轮完成：** [原空间辅助测度的反射内积与完整球面积分]({p}research_note_657.md)'
         '原自由偶数时间盒的带相位辅助Pfaffian给同一BCS反射Gram核；全S⁹积分为显式平方范数。'
         '两时间片严格归一，任意偶数长度未归一正性成立。'
         '四组、二十一式通过，最新657／3214，1283份编号科学文件、2239份保护证据。'
         '[核验]({p}research_round_657_checks.json)、[条件账]({p}unified_physics_condition_ledger_657.md)。'
         '一般长度严格归一、共同单步过程、物理CAR与量子GR仍开放。')
order=('**当前执行顺序（657后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[658严格归一与跨时间共同过程]({p}round658_drafts/STATUS.md)，'
       '核原平均向量支撑及不同时间盒的实际映射；不把每盒正内积当同一演化，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原空间测度与共同传播','原空间测度的反射内积')
publish=publish.replace('空间测度、时间耦合与来源须共同签收','原测度、反射内积与完整球面平均须共用对象')
publish=publish.replace('原空间手征测度与共同传播投影子的连接','原空间辅助测度的反射内积与完整球面积分')
write('publish_round657.py',publish)
write('postcheck_round657.py',replace((HERE/'postcheck_round656.py').read_text('utf8'),mapping))
write('round657_drafts/research_note_657_draft.md',(HERE/'research_note_657.md').read_text('utf8'))
review="""657 primary-agent review; no independent agent review.
Previous goal turn was progress:656 completed/published/postchecked and657
actual signed reflection-kernel probe saved. Current turn read root/research
navigation,656 note/ledger/results,657 entry and reverified all2229 evidence.
No live Python process found, no goal change, task, automation or image check.
Historical search and612/646/653/656 reviewed. Free reflection positivity and
BCS overlap machinery not counted as new; original auxiliary mapping is new.
Primary Kikukawa-Usui1005.3751v3 sectionIII equations30/31 opened and mapped:
C=Gamma-dagger P(theta x,y)=-gamma4 Dov(theta x,y), Dov action not inverse.
Both finite AP spectral terms positive after reflection. Rectangular spatial
boxes only change finite Fourier sum. m0=1 exceptional zero spatial momentum
handled by m0 approaching1 in finite gapped AP box, not discarded modes.
S=Rtheta tensor gamma5 gamma4 complements P. After folding lower spin, block
P is [Delta,C;C,I-Delta]; P^2 and C positivity give canonical square root.
Gamma-transpose B Gamma=B preserves original scalar E reflection contract.
Original M(E) unitary skew, constant Pfaffian from connected S9 and nonzero
determinant. No continuity across zeros used for a sign claim.
BCS overlap identity retains sign through complementary Pfaffian expansion;
Robledo0901.3213 checked as mature source, not asserted new. Full Fock six-mode
enumeration independently verifies transposes, complex rotations and endpoints.
WDelta defined in eigenmodes of Delta TRANSPOSE, with vacuum d^(1/4) and
occupied (1-d)^(1/4); endpoint formula continuous, no division by zero.
Original signed Pfaffian, not modulus, equals Gram. A positive kernel need not
be pointwise positive. Sphere average uses exact original full S9 moments and
finite pure-creation polynomial through degree64 per site, not planar measure.
Full integral is squared norm. General semidefinite W may annihilate the mean:
strict general-time normalization NOT inferred from epsilon positive limit.
Two-time Delta strictly between0/1; vacuum component gives actual3-space-site
Z>=2^-160. This is loose bound, not a numerical partition value.
Global original G only, no local interacting Gauss claim. Same W transports
kappa source; kappa not metric stress. Auxiliary64CAR per half-site distinct
from original physical32CAR. Same positive kernel for each box does not prove
common one-step transfer; at kappa0 two-slice weight is k squared, not T^2.
Numerics: original full S9 pairings2/4/6 time slices, projector blocks through8;
square roots near endpoints have O(sqrt(machine epsilon)) error recorded,
not an artificial physical truncation. Abstract six-mode Fock and exact
degree4 sphere quadrature only test sign/moment formula, not full huge Fock.
4 groups,21 equations. Next658 strict normalization and actual across-time
mapping. Cognitive architecture design deferred; full unified goal active.
"""
for name in ('research_note_657.md','joint_auxiliary_reflection_gluing.py',
             'joint_auxiliary_reflection_gluing_results.json','unified_physics_condition_ledger_657.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round657_drafts/final_review.txt',review)
print('657 frozen; verification and publication prepared')
