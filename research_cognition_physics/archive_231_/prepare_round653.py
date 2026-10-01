"""Freeze 653 actual chiral auxiliary character and temporal transfer."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        pat=re.escape(key)
        if key.isdecimal():pat=r'(?<!\d)'+pat+r'(?!\d)'
        pats.append(pat)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_quantum_reference_forms':'joint_chiral_character_state',
    '651':'652','652':'653','653':'654','3191':'3194','3194':'3198',
    '1265':'1268','1268':'1271','2175':'2186','2186':'2196',
    'common_quantum_reference_entry.md':'common_chiral_state_entry.md',
    'curved_reference_velocity_probe':'holonomy_charge_probe'}
verify=replace((HERE/'verify_round652.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original quantum reference forms and common energy representation.',
                      'Verify original subgroup auxiliary character and actual temporal transfer.')
start=verify.index('    names=(');end=verify.index('\n    preserved=',start)
verify=verify[:start]+"""    names=('unified_physics_condition_ledger_653.md','round653_drafts/research_note_653_draft.md',
           'round653_drafts/final_review.txt','round654_drafts/STATUS.md',
           'round653_drafts/common_chiral_state_entry.md','round653_drafts/holonomy_charge_probe.py',
           'round653_drafts/holonomy_charge_probe_results.json')"""+verify[end:]
verify=verify.replace("'round653_drafts/integration_before_design.md',",'')
verify=verify.replace("==(3,0,0)","==(4,0,0)").replace('fresh_tests=dict(run=3,','fresh_tests=dict(run=4,')
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==20")
write('verify_round653.py',verify)
publish=replace((HERE/'publish_round652.py').read_text('utf8'),mapping|{'## 298.':'## 299.','## 203.':'## 204.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第653轮完成：** [原手征辅助测度、规范状态与时间拼接的共同表示]({p}research_note_653.md)'
         '原32CAR直接迹接法有电荷谱障碍；原整个商群的辅助积分有正表示，'
         '实际空间平凡时间链给同一正转移、规范投影与来源。'
         '四组、二十式通过，最新653／3198，1271份编号科学文件、2196份保护证据。'
         '[核验]({p}research_round_653_checks.json)、[条件账]({p}unified_physics_condition_ledger_653.md)。'
         '限零标量时间分支；原全图动力学、连续手征和引力仍开放。')
order=('**当前执行顺序（653后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[654原空间传播、质量与正转移]({p}round654_drafts/STATUS.md)，'
       '检验原全图对象与本次手征状态／演化的实际连接；不优化辅助维数，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原量子参考与同一能量表示','原辅助测度与规范时间拼接')
publish=publish.replace('参考量子表示须与原过程共同相容','手征状态须与测度及时间拼接共同相容')
publish=publish.replace('原物质量子参考、共同能量与区域表示的联合连接','原手征辅助测度、规范状态与时间拼接的共同表示')
write('publish_round653.py',publish)
write('postcheck_round653.py',replace((HERE/'postcheck_round652.py').read_text('utf8'),mapping))
write('round653_drafts/research_note_653_draft.md',(HERE/'research_note_653.md').read_text('utf8'))
review="""653 primary-agent review; no independent agent review.
Previous goal turn is progress:652 completed, published, postchecked;653 exact
charge support probe and entry saved. Current turn revalidated navigation,
latest652 note/results/ledger,2186 protected evidence, and653 pending drafts.
No live Python process or ongoing handle found. No task/automation/goal edit.
Read actual614 representation,615 Pfaffian and measure,616 thermal boundary,
628/629 anomaly/phase interfaces and646 free-window limit. Searched existing
notes for auxiliary characters/complete homogeneous/495/D=32440320; no prior
completed interface found. Ordinary character theory is not a new theorem.
Primary Kikukawa3.2/3.6/5.2, Etingof Weyl character notes, and a primary paper's
Funk-Hecke appendix checked. The old candidate's general locality problem is
not claimed solved. No screenshots/images.
Original32 CAR conditional trace has Q support[-36,36] for ALL theta-independent
positive states.615 full factor has exact rational coefficient at60. Obstruction
is to same-Q conditional trace, NOT to full Gauss trace or all chiral models.
General original subgroup flat one-cell Pfaffian reduces to |C_alpha E|^16;
Dirichlet(1^5) integral gives h8(a)/495. detV=1 excludes V=-I5, so M>0 in this
whole one-cell family. This is not a general configuration nonzero theorem.
Integer Laurent polynomial, Weyl alternating identity and all495 subgroup
irreps certified exactly. All multiplicities positive and Z6 conditions hold.
Independent spherical-harmonic weighted character identity gives a structural
positive representation, not only sampled Gram positivity.
NEW actual temporal-chain result: m0=1, spatial links identity, one cell in each
spatial direction, zero scalar. Covariant AP shift S is unitary; X=-Sdag P+ -S P-.
Actual chiral frames u/v use same S. Pairing reduces to block average of T(E)
and S* T(E) Sdag. This yields K(E_x,O(g_x)E_next)=((1+dot)/2)^8.
O(g) defined via R(g)* T(E) R(g)dag is a representation; order agrees with
shift holonomy g0...gN-1. AP signs square out of auxiliary pairing only.
All N gluing follows actual kernel integration, not a fitted transfer matrix.
Funk-Hecke eigenvalues checked by exact Gegenbauer/Beta moments. Positive
support=sum harmonics l0..8 has35750 dimensions; traceT=1. Full L2 has a kernel,
so logT is defined only on support. No complete inserted-operator reconstruction.
Physical determinant retains 2^(-32N); detR=1 removes direction conjugation in
the square. H_* includes 32log2/a_tau identity, no lost vacuum normalization.
Actual 1,2,3-time-slice matrices independently verify projection, Pfaffian and
physical determinant; no repeated one-point interpolation presented as proof.
Same original G projection exists on EXTENDED Fock x support and gives a
positive state. This space/H_* is not original598 full matter H or603 Gibbs.
Auxiliary dimension is not a count of physical particle species. No new
physical gauge group, cognitive axiom, hardware or optimized implementation.
Source matching is original holonomy response only; not GR stress or full
spatial/finite-mass response. Next restores those true joint conditions.
All four groups passed stated thresholds.20 numbered formulas, exact domain,
normalization, branch assumptions and negative/positive conclusions audited.
"""
for name in ('research_note_653.md','joint_chiral_character_state.py',
             'joint_chiral_character_state_results.json','unified_physics_condition_ledger_653.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round653_drafts/final_review.txt',review)
print('653 frozen; verification/publication prepared')
