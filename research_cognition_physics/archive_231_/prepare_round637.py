"""Freeze 637: complete-model regional entropy via original energy control."""
import hashlib,re
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
mapping={'joint_quotient_boundary_entropy':'joint_full_thermal_region_entropy',
         '635':'636','636':'637','637':'638',
         '3146':'3149','3149':'3151','1217':'1220','1220':'1223',
         '2043':'2050','2050':'2057'}
verify=replace((HERE/'verify_round636.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original quotient labels, regional entropy and electric-source witness.',
                      'Verify complete-model regional entropy and original quotient heat bounds.')
verify=verify.replace("(3,0,0)","(2,0,0)").replace("fresh_tests=dict(run=3,","fresh_tests=dict(run=2,")
write('verify_round637.py',verify)
publish=replace((HERE/'publish_round636.py').read_text('utf8'),
                mapping|{'## 282.':'## 283.','## 187.':'## 188.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第637轮完成：** [原完整Gauss热态、区域熵与同尺度近似的共同能源条件]({p}research_note_637.md)'
         '原完整H能源下界经同一切分控制区域热迹与有限熵；'
         '原有界注能记录及保能源的同图全态近似共同覆盖。'
         '两组、十四式通过，最新637／3151，1223份编号科学文件、2057份保护证据。'
         '[核验]({p}research_round_637_checks.json)、[条件账]({p}unified_physics_condition_ledger_637.md)。'
         '限固定图与正背景；空间连续、面积律、Newton匹配及GR仍开放。')
order=('**当前执行顺序（637后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[638跨尺度区域态、规范边界与几何熵]({p}round638_drafts/STATUS.md)，'
       '联查旧细化障碍、原区域态及635局部几何项，核共同有限量；目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原商群边界记录与电来源的共同限制','原完整能源与区域熵的共同成立')
publish=publish.replace('边界熵约束不等于连续几何已涌现','固定图有限区域熵不等于连续面积律')
publish=publish.replace('原商群的联合边界记录、区域熵与电几何来源','原完整Gauss热态、区域熵与同尺度近似的共同能源条件')
write('publish_round637.py',publish)
write('round637_drafts/research_note_637_draft.md',(HERE/'research_note_637.md').read_text('utf8'))
review="""637 primary-agent review. No independent agent review.
636 completed, verified, published and postchecked before beginning 637.
One 636 postpublication helper initially used Path('.') without resolving,
so its parent lookup failed after publication; corrected read-only helper
verified all 2050 protected files, 7 snapshots and 7 published navigation files.
603 full thermal proof, 598 full CAR bound, 617 cut and 636 were reread.
History search found no prior full-H regional entropy proof to repeat.
Winter 1507.07775 Gibbs hypothesis and Lemma15 were read as inherited tools.
The new model-specific work is full original coercivity plus regional
comparison construction, exact half-Casimir weights and quotient heat trace.
Full non-diagonal geometry, confining five-scalar potential and complete
finite CAR are retained. Comparison K is NOT a replacement physical H.
H+Cstar >= K0 is a quadratic-form statement; no invalid exponential
operator monotonicity is used. Entropy bound follows finite truncations.
J*(KA+KB)J=K0 counts every original edge once, including cuts.
Positive KA,KB give regional energy bounds for arbitrary finite-energy
normal original Gauss states, not just selected loop states.
Scalar heat trace uses inherited exact H5 kernel and original U integral.
Quotient heat trace uses squared representation dimensions and the actual
Z6 congruence, independently checked by the six center characters.
Rectangle tails use monotone polynomial Gaussian ratios and a union bound.
Analytic truncation tails do not include floating-point rounding error.
Internal Gauss projection commutes with reference K and cannot enlarge trace.
The regional reference Gibbs state is not the original physical Gibbs state.
Full actual Gibbs energy is finite by 603; it is not numerically computed here.
Original bounded-injection L+/- instrument preserves finite average energy;
conditional bound includes 1/p and never assigns a state to p=0.
Uniform-energy trace continuity uses the SAME finite graph and SAME KA.
Actual full-H Gibbs spectral projections preserve Gauss and have uniform
energy, so their regional entropy convergence is analytic, not fictitious
numerical evaluation of the enormous original Gibbs density matrix.
Numerics check original slab counts, actual coupling/geometry constants,
quotient trace, and entropy bounds on the 636 allowed loop fixtures only.
Two groups passed first execution; fourteen formulas checked textually.
No continuum-uniform estimate, area law, physical Newton value, autonomous
instrument or dynamic GR is claimed. Next638 is cross-scale joint matching,
not a continuation of device design or pure-electric model substitution.
"""
for name in ('research_note_637.md','joint_full_thermal_region_entropy.py',
             'joint_full_thermal_region_entropy_results.json','unified_physics_condition_ledger_637.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round637_drafts/final_review.txt',review)
print('637 science frozen; verifier and publisher prepared')
