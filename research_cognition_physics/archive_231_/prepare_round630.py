"""Freeze substantive round 630; preserve all historical evidence."""
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

mapping={'joint_topological_mass_phases':'joint_continuum_source_spectrum',
         '628':'629','629':'630','630':'631',
         '3125':'3127','3127':'3130','1196':'1199','1199':'1202',
         '1994':'2001','2001':'2008'}
verify=replace((HERE/'verify_round629.py').read_text('utf8'),mapping)
verify=verify.replace('Verify joint original mass phases and quotient-group topological weights.',
                      'Verify original full matter cuts, UV matching and joint metric sources.')
verify=verify.replace("==(2,0,0)","==(3,0,0)").replace("run=2,failures=0","run=3,failures=0")
verify=verify.replace("['display_formulas']==14","['display_formulas']==19")
write('verify_round630.py',verify)
publish=replace((HERE/'publish_round629.py').read_text('utf8'),
                mapping|{'## 275.':'## 276.','## 180.':'## 181.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第630轮完成：** [原物质谱、真实时间响应与几何来源的共同连接]({p}research_note_630.md)'
         '原整代径向质量的正实时谱复现两个旧UV系数，并共同固定真空噪声与共形几何交叉来源；'
         '最低成对阈值限制同一局部展开窗口。'
         '三组、十九式通过，最新630／3130，1202份编号科学文件、2008份保护证据。'
         '[核验]({p}research_round_630_checks.json)、[全账]({p}unified_physics_condition_ledger_630.md)。'
         '限明确连续真空分支；一般张量应力、图映射、实际态和GR仍开放。')
order=('**当前执行顺序（630后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[631完整物质与张量几何响应]({p}round631_drafts/STATUS.md)，'
       '核其它度规通道能否与同一物质、态及旧几何作用共用响应；不继续脉冲精度优化，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原质量、拓扑权重与参数冗余','原物质谱与几何来源的共同连接')
publish=publish.replace('相位参数商不是空间维数','共形响应不是一般度规或空间维数')
publish=publish.replace('原质量相位、拓扑权重与共同量子测度的压缩','原物质谱、真实时间响应与几何来源的共同连接')
write('publish_round630.py',publish)
write('round630_drafts/research_note_630_draft.md',(HERE/'research_note_630.md').read_text('utf8'))
review="""630 primary-agent review. No independent agent review.
Previous conversational turn only reconfirmed ordering and was not a scientific
progress turn. This round supplies a new original-model spectral interface.
629 receipt, frozen entry, reports and prior evidence are preserved.
Historical checks excluded repeating generic Ward relations, counterterms,
scalar heat-kernel derivative orders, anomalies and finite-source cross noise.
Djouadi hep-ph/0503172 sec 2.1.1 eq 2.6 was read to cross-check the beta^3
threshold. No outdated experimental data were imported. Normalization is
derived independently by canonical Dirac projectors and phase-space delta.
599/Vassilevich eq 4.28 supplies the already-established heat-kernel input.
All sixteen original Weyl masses use half weight; 7 Dirac plus 2 Majorana.
Original 64-dimensional BdG positive eigenvalues match twice each singular
value. No additional spin or Nambu multiplicity is counted.
Same 599 ray and original complex Yukawas; q=.27 is diagnostic, not vacuum fit.
The continuum spacetime, canonical kinetic operator and vacuum are inputs.
Free fermion bilinear correlations are exact in that chosen background;
external source effective action is evaluated to quadratic source order
(one fermion loop), not the full interacting SM or graph Gauss process.
All masses scale on the radial ray, so one constant Takagi basis suffices.
Spectral convention Pi_R=+i theta<[O,O]> is stated; for +eta O the response
is -Pi eta, and Gamma_E''=-Pi(iQ) for linear mass source.
Direct sea-energy curvature and dispersive zero-frequency integral agree.
Frequency cutoff quadratics are scheme dependent. Only logarithms are matched
to dimensional regularization (1/epsilon corresponds to 2 log cutoff).
Both source Hessian factors and Euclidean-frequency sign are included.
Twice-subtracted response converges; c0,c2 and any allowed finite higher local
terms remain explicit. No spectral discontinuity can be changed by a real
local polynomial alone. Explicit local matter fields remain a valid option.
L4 and L6 are analytic beta-integrals. Bounds hold below the smallest positive
pair threshold. Small-mass and small-frequency limits cannot be interchanged.
Band-limited and compact-time sources are not identified.
Gaussian work is a diagnostic per volume and per epsilon squared, with an
exponentially controlled integration tail; no apparatus or exact compact
time support is asserted. This vacuum has no thermal zero-frequency occupancy
noise; old 602/626 thermal memory is not discarded or contradicted.
Conformal fermion rescaling gives e^sigma(1+eta) mass. Only nonlocal
connected quadratic kernels have rank one with vector (b,1).
Local trace anomaly, renormalization and tadpole/coordinate contacts remain.
The compensated direction does not leave full scalar/gravity action invariant.
All-source tensor stress, other scalar directions, real gauge interactions,
state matching and graph-to-continuum limits remain open.
27-condition ledger is a synthesis and adds no scientific round by itself.
Next 631 returns to full geometric channels, not source-precision design.
Three test groups passed on the first execution. Report has nineteen formulas.
"""
for name in ('research_note_630.md','joint_continuum_source_spectrum.py',
             'joint_continuum_source_spectrum_results.json','unified_physics_condition_ledger_630.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round630_drafts/final_review.txt',review)
print('630 verification and publication prepared; science frozen')

