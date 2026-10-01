"""Freeze full Gauss regional replicas and source matching."""
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

mapping={'joint_gauss_fermion_influence':'joint_region_replica_source',
         '642':'643','643':'644','644':'645',
         '3165':'3168','3168':'3171','1238':'1241','1241':'1244',
         '2097':'2104','2104':'2111'}
verify=replace((HERE/'verify_round643.py').read_text('utf8'),mapping)
verify=verify.replace('Verify full Gauss heat influence and original32-CAR conditional weights.',
                      'Verify original full Gauss replica identities and64-CAR diagnostics.')
write('verify_round644.py',verify)
publish=replace((HERE/'publish_round643.py').read_text('utf8'),
                mapping|{'## 289.':'## 290.','## 194.':'## 195.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第644轮完成：** [同一Gauss区域态的费米复制、二阶熵与几何来源]({p}research_note_644.md)'
         '原完整区域态、规范边界及CAR符号共用复制热历史，完整来源同时给二阶熵变分；'
         '两节点64模式的条件计算三路相符，原合法Gauss态另核符号。'
         '三组、十四式通过，最新644／3171，1244份编号科学文件、2111份保护证据。'
         '[核验]({p}research_round_644_checks.json)、[条件账]({p}unified_physics_condition_ledger_644.md)。'
         '数值非全Gauss热积分；面积律、连续匹配和GR未完成。')
order=('**当前执行顺序（644后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[645区域态与连续几何熵的共同尺度]({p}round645_drafts/STATUS.md)，'
       '核参考、区域、复制参数和尺度的真实映射；不继续小型热因子扫描，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('完整Gauss物质的共同热历史与来源','完整Gauss区域态的二阶熵与共同来源')
publish=publish.replace('完整热桥表示不替代连续几何推导','区域复制一致不等于面积熵或GR成立')
publish=publish.replace('完整Gauss物质、共同热参考与来源的有序历史表示','同一Gauss区域态的费米复制、二阶熵与几何来源')
write('publish_round644.py',publish)
write('round644_drafts/research_note_644_draft.md',(HERE/'research_note_644.md').read_text('utf8'))
review="""644 primary-agent review; no independent agent review.
643 was verified, published and postchecked before this independent next unit.
617's same gamma-independent J is retained. K=J P exp(-beta H) Jdag is zero
outside the matched physical space; no redundant thermal states are counted.
The Hilbert swap identity is established mathematics. Its conversion to CAR
second quantization is derived with the actual fermion parity convention.
r_A sends copy1 to copy2 and copy2 to minus copy1; det r_A=1. Its lift produces
ordinary purity for even reduced states, including number-nonconserving states.
Unsigned single-particle exchange gives parity-weighted purity, NOT ordinary
Hilbert swap. Fixed CAR ordering and graded product of even copy states justify
restriction to A; no arbitrary Jordan-Wigner sign is inserted.
Recent preprint2605.21632 was inspected: its displayed Koszul-only exchange,
when read using our ordinary Hilbert trace, fails the occupied one-mode check.
Its formulas are not imported as a theorem here; conventions are resolved by
the self-contained parity proof and explicit Fock test, without a claim to audit
the entire preprint. Casini/Huerta and Klich are used within stated tool scope.
Both full Gauss projections remain in the replica weight. They do not generally
commute with regional exchange. Same bosonic endpoint rewiring, original CAR
closure and J multiplication are used; there is no independent regional bath.
Absolute bridge integrability follows from643 scalar domination and product
measure Cauchy-Schwarz after swapping region coordinates. Haar cut volumes are1.
Trace-norm Duhamel uses623 operator-relative bounds, stronger than only603 form
means, splitting the integral at beta/2. Region and J remain fixed under gamma.
Normalization derivative has the correct -2 beta <G> sign in S2'. Full G keeps
kinetic, potential, CAR and non-diagonal gauge terms. No Brownian density score.
Numerics freeze original two scalar backgrounds and identity gauge link, and
explicitly CHOOSE allowed hopping t=.17. All64 original CAR modes, complex Y,
Majorana and nonzero inter-node hopping remain; no claim that t was derived.
Exact16/256-dimensional Fock factorization with multiplicities sums to32 modes
per node; independent128-dimensional Nambu covariance and256 replica determinant
agree. The state is CONDITIONAL, not the full Gauss bosonic thermal integral.
642 physical particle-hole states separately verify odd regional parity with
original boundary Schmidt rank d_R. They are normal finite-energy Gauss states,
not full-H Gibbs/eigenstates. Old state construction itself not counted as new.
Replica2 is not alpha1 geometric entropy. No area law, Newton coupling, continuum
chiral limit, dynamic geometry or GR theorem claimed. Next step retains that gap.
Three groups pass; goal unchanged; cognition implementation remains deferred.
"""
for name in ('research_note_644.md','joint_region_replica_source.py',
             'joint_region_replica_source_results.json','unified_physics_condition_ledger_644.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round644_drafts/final_review.txt',review)
print('644 science frozen; verifier and publisher prepared')
