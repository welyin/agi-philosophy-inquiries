"""Freeze656 and prepare bounded verification/publication; exclusive writes."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    patterns=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        patterns.append(p)
    return re.sub('|'.join(patterns),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_full_graph_transfer_sources':'joint_spatial_auxiliary_geometry',
    '654':'655','655':'656','656':'657','3202':'3206','3206':'3210',
    '1274':'1277','1277':'1280','2206':'2218','2218':'2229',
    'normalization_geometry_entry.md':'actual_spatial_measure_entry.md',
    'lapse_extension_probe':'spatial_auxiliary_probe'}
verify=replace((HERE/'verify_round655.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original full graph transfer, lapse sources and geometric normalization.',
                      'Verify original spatial auxiliary measure, shared propagation and response.')
verify=verify.replace("'round656_drafts/round655_formula15_erratum.md','round656_drafts/development_diagnostic.md'",
                      "'round656_drafts/initial_planar_results.json'")
verify=verify.replace("text['display_formulas']==22","text['display_formulas']==20")
write('verify_round656.py',verify)
publish=replace((HERE/'publish_round655.py').read_text('utf8'),mapping|{'## 301.':'## 302.','## 206.':'## 207.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第656轮完成：** [原空间手征测度与共同传播投影子的连接]({p}research_note_656.md)'
         '原自由P给全S⁹辅助权重模长、平面精确正性及九方向共同Laplacian。'
         '恢复空间同时改写时间耦合；权重与其来源不能另拟。'
         '四组、二十式通过，最新656／3210，1280份编号科学文件、2229份保护证据。'
         '[核验]({p}research_round_656_checks.json)、[条件账]({p}unified_physics_condition_ledger_656.md)。'
         '限自由被积式；全S⁹积分、正物理过程、连续与量子引力仍开放。')
order=('**当前执行顺序（656后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[657原空间积分与正时间拼接]({p}round657_drafts/STATUS.md)，'
       '核实际测度能否共用正过程与原态；不把平面正性等同完整量子重建，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原全图过程与几何归一','原空间测度与共同传播')
publish=publish.replace('全图演化、时间来源与体积归一须共同签收','空间测度、时间耦合与来源须共同签收')
publish=publish.replace('原完整图的质量转移、共同时间来源与几何归一边界','原空间手征测度与共同传播投影子的连接')
write('publish_round656.py',publish)
write('postcheck_round656.py',replace((HERE/'postcheck_round655.py').read_text('utf8'),mapping))
write('round656_drafts/research_note_656_draft.md',(HERE/'research_note_656.md').read_text('utf8'))
review="""656 primary-agent review; no independent agent review.
Previous conversational goal turn was no scientific progress: it checked the
latest ledger and restated priority. This turn revalidated655 and resumed the
available656 calculation. All2218 protected prior evidence hashes intact.
Root/research READMEs/direction/state, latest655, ledger, entry/results read.
No live Python process found; no goal change, tasks, automation or image check.
Historical search and615/653/655 read. Mature Slater determinant and overlap
locality not claimed as new. Kikukawa1710.11618v3 section4.2 checked: global
free-field S9 positivity argument uses realness and zero-locus assumptions;
we do not use it as an unconditional theorem. No images inspected.
Original606/615 B,T and negative spectral projector retained. Free real shifts
give B-dagger P-star B=P. A0 unitary follows, compression formula then exact.
K(E) unitary for every S9, not only the two-plane. Full S9 modulus identity
det(I-L)^(1/4), leakage bound and pair-distance trace use original matrices.
Modulus is never substituted for signed measure. General phase unresolved.
Planar signed Pfaffian selected by Laurent polynomial identity, not unsafe
continuity across determinant zeros. Sign and zero example checked directly.
Equality criterion uses finite-dimensional invariant/reducing subspace and
injectivity of K(E). Constant maximum does not assert ordered phase/vacuum.
Nine-tangent Hessian uses sphere second derivative; diagonal fixed by P^2=P.
Six-site Fourier weights, spectrum and kappa derivatives checked separately
against original Pfaffian, matrix Hessian and spectral divided differences.
Counterexample excludes only normalized pure equal-time spatial multiplier;
does not exclude other transfer representations or all integrated measures.
Full4D m0=1 gap1/norm7 from exact symbol identity; binomial tail is volume-
uniform but loose. Does not imply hard locality, Lorentz causality or emergent
dimension. Kappa remains a Wilson parameter, not automatically metric stress.
Euclidean Slater auxiliary state not identified with original physical32CAR.
Initial planar-only results preserved; completed new identity added to final
results with exclusive output creation. No failed numerical tests concealed.
4 meaningful test groups,20 equations; checks are finite original matrix
instances, while stated general free identities have analytic proofs.
Next657 asks actual spatial measure positive gluing/S9 integration, not
parameter optimization or cognitive system design. Full goal stays active.
"""
for name in ('research_note_656.md','joint_spatial_auxiliary_geometry.py',
             'joint_spatial_auxiliary_geometry_results.json','unified_physics_condition_ledger_656.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round656_drafts/final_review.txt',review)
print('656 draft frozen; verifier and publication scripts prepared')
