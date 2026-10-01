"""Freeze658 actual sphere-average support and prepare sequential publication."""
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


mapping={'joint_auxiliary_reflection_gluing':'joint_auxiliary_normalization_support',
    '656':'657','657':'658','658':'659','3210':'3214','3214':'3217',
    '1280':'1283','1283':'1286','2229':'2239','2239':'2246'}
verify=replace((HERE/'verify_round657.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original auxiliary reflection Gram identity and S9 integration.',
                      'Verify exact zero-sector support and strict original auxiliary normalization.')
a=verify.index('    import importlib.util');b=verify.index('    result=core.read(model.TARGET)',a)
verify=verify[:a]+verify[b:]
a=verify.index('    names=(');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_658.md','round658_drafts/research_note_658_draft.md',
           'round658_drafts/final_review.txt','round659_drafts/STATUS.md')
"""+verify[b:]
verify=verify.replace("text['display_formulas']==21","text['display_formulas']==16")
verify=verify.replace('(4,0,0)','(3,0,0)').replace('run=4','run=3')
verify=verify.replace("'round659_drafts/STATUS.md',\n                 'round658_drafts/positive_gluing_entry.md'", "'round659_drafts/STATUS.md'")
write('verify_round658.py',verify)
publish=replace((HERE/'publish_round657.py').read_text('utf8'),mapping|{'## 303.':'## 304.','## 208.':'## 209.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
header="""summary=('**第658轮完成：** [原球面平均的零模支撑与严格归一]({p}research_note_658.md)'
         '原半区零模只在空间零动量，原S⁹平均在强制占据部门有精确正开链系数。'
         '全部有限自由空间／偶数时间辅助积分严格正，657反射内积可实际归一。'
         '三组、十六式通过，最新658／3217，1286份编号科学文件、2246份保护证据。'
         '[核验]({p}research_round_658_checks.json)、[条件账]({p}unified_physics_condition_ledger_658.md)。'
         '共同单步过程、物理CAR、一般规范场与连续引力仍开放。')
order=('**当前执行顺序（658后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[659跨时间共同演化]({p}round659_drafts/STATUS.md)，'
       '核实际不同时间盒及原过程映射；不优化归一下界或装置，目标不改。')
"""
publish=publish[:a]+header+publish[b:]
publish=publish.replace('原空间测度的反射内积','原球面平均的零模支撑')
publish=publish.replace('原测度、反射内积与完整球面平均须共用对象','原平均、端点支撑与归一反射内积共同成立')
publish=publish.replace('原空间辅助测度的反射内积与完整球面积分','原球面平均的零模支撑与严格归一')
write('publish_round658.py',publish)
write('postcheck_round658.py',replace((HERE/'postcheck_round657.py').read_text('utf8'),mapping))
write('round658_drafts/research_note_658_draft.md',(HERE/'research_note_658.md').read_text('utf8'))
review="""658 primary-agent review; no independent agent review.
Continues completed657 in the same goal turn, after its full verifier,
publication and postcheck of2239 protected evidence and8556 navigation links.
No goal change, new task, automation, image check or concurrent writer.
Reuses653 exact open/closed temporal kernels and657 free reflection spectrum.
New gap is nonvanishing of original sphere mean under original WDelta support.
For nonzero spatial momentum b>0, both free spectral matrix eigenvalues are
strictly positive above cut threshold: squared difference=2b coshE-A>0.
Finite nonzero temporal exponential polynomial cannot vanish on interval;
cross-reflection C(k)>0, hence0<Delta(k)<1. Original k=0 is the only endpoint.
Explicit uniform spatial Wplus-Wminus internal temporal-edge frame spans all
2(L-1) null spin modes; complementary plus frame eigen1 and4 boundary modes
eigen1/2. Analytic frame, not numerical threshold, determines null space.
Restricted original M(E) is block diagonal with T((mean E_j+mean E_j+1)/2).
Spatial means NOT renormalized. Signed Pfaffian constant reference phase
times product of sixteenth powers; identity inherits original Clifford data.
Sphere mean coefficient I0>0 by nonnegative integrand and positive-measure
neighborhood. This is one component of full original mean, not replacing
the whole graph with spatial averages. Fock modes for Delta TRANSPOSE mean
null frame is Zstar, hence restricted pairing is Z-transpose M Z.
All zero modes occupied, every other mode empty is allowed support and has
nonzero mean component; squared weight=(pdet Delta_spin)^8 I0^2. No other
orthogonal Fock components can cancel norm. All finite even times now strict.
Lift Delta+Z Zdagger computes pseudodeterminant only; no mass added to model.
S9 cap density128/(35pi)*(1-u^2)^(7/2); interval[1/2,3/4] gives explicit
positive lower bound. Cap is proof subset, not threshold or changed measure.
No claim of uniform volume/time/continuum lower bound or actual Z numerical
value. One-space-site exact lambda0 open chain and old closed spectrum check
power exponents with Fractions; mature spectrum not counted again as new.
3 test groups,16 equations. Full finite free reflection functional normalized;
common one-step dynamics, physical32CAR and interacting gauge/GR still open.
Next659 actual cross-time mapping; no optimization of the cap or resources.
"""
for name in ('research_note_658.md','joint_auxiliary_normalization_support.py',
             'joint_auxiliary_normalization_support_results.json','unified_physics_condition_ledger_658.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round658_drafts/final_review.txt',review)
print('658 frozen; verification/publication prepared')
