"""Freeze 620 primary review and prepare guarded append-only publication."""
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
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)

mapping={'joint_chiral_boundary_gluing':'joint_quantum_exchange_noise',
         '618':'619','619':'620','620':'621',
         '3098':'3101','3101':'3104','1166':'1169','1169':'1172',
         '1922':'1930','1930':'1937'}
verify=replace((HERE/'verify_round619.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original Nambu boundary ranks, transmission, and frame sources.',
                      'Verify same-state joint source noise and actual energy exchange.')
verify=verify.replace(",'round620_drafts/boundary_current_entry.md'",'')
assert 'boundary_current_entry' not in verify
write('verify_round620.py',verify)
publish=replace((HERE/'publish_round619.py').read_text('utf8'),
                mapping|{'## 265.':'## 266.','## 170.':'## 171.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第620轮完成：** [共同量子来源的交换约束、混合涨落与框架匹配]({p}research_note_620.md)'
         '原32模式同态来源变换与实际能量交换共同固定交叉关联；'
         '删掉混合噪声会破坏框架方差和交换身份。'
         '三组、十四式通过，最新620／3104，1172份编号科学文件、1937份保护证据。'
         '[核验]({p}research_round_620_checks.json)、[条件账]({p}unified_physics_condition_ledger_620.md)。'
         '限条件质量部门及明示Ward前提；全Gauss连续来源、动态几何与GR仍开放，无新增独立代理审查。')
order=('**当前执行顺序（620后，优先于下方历史安排）：** 继续先整合共同条件，认知设计后置。'
       '接[621全Gauss过程、涂抹来源与联合响应域]({p}round621_drafts/STATUS.md)，'
       '复用有限能源／零频旧界，核实际联合来源域与连续接口；目标不变。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('手征区域、法向流与共同来源','共同来源、混合涨落与交换约束')
publish=publish.replace('连续旋量边界要求不推出维数','联合量子来源条件不选择维数')
publish=publish.replace('原手征物种、边界电流与共同区域拼接','共同量子来源的交换约束、混合涨落与框架匹配')
assert 'hashlib.sha256' in publish
write('publish_round620.py',publish)
write('round620_drafts/research_note_620_draft.md',(HERE/'research_note_620.md').read_text('utf8'))
review='''620 primary-agent review. No independent agent review.
Previous turn 619 publication was verified progress; 620 adds actual
same-model joint noise data rather than another CTP overview.
Reused 590 finite-energy noise obstruction, 600 mean Ward/source reduction,
602 static/isolated mismatch and Wick half weight, 603 thermal response
domain, and 618-619 frame/boundary conventions. No historical files edited.
The Lorentz lower-metric source convention is explicit. The generating
functional yields a mean Ward identity, not an operator identity by itself.
Noise claims assume a valid renormalized operator Ward identity and test
function domains. Nonzero gauge curvature and boundary flux cannot be
discarded; the displayed continuum branch has zero external gauge curvature.
Gaussian stochastic representation is restricted to symmetric two-point
data, with no claim about higher cumulants, commutators, or full processes.
Same fixed CAR mass family is pulled back using n_E=sqrt(F)n_J. The s
source mixes scalar and lapse sources; the noise transforms by the same
Jacobian. This is not a full continuum frame equivalence or a Hessian rule.
At fixed Higgs sqrt(F)B has only the original singlet Majorana derivative;
the code retains all 32 original modes, complex Yukawas and pairing.
The state is identical during the pullback. During the driven diagnostic
one initial Gibbs state evolves unitarily without instantaneous resetting.
d(U^dagger B U)/dt=sdot U^dagger B' U follows directly, including noncommuting
masses. The integrated operator is not a two-projective-measurement work law.
W is independently midpoint-integrated; DeltaE uses endpoints. Observed
second-order convergence is not a rigorous interval error certificate.
Symmetric Wick covariance has factor 1/2 for the doubled Nambu matrix.
Independent 16-dimensional neutral Fock dynamics verifies means and all
tested covariances. The full 64-matrix computation retains charged modes.
Dropping mixed covariance violates the exact exchange identity when its
variance is nonzero. Positivity alone would not exclude the independent
matrix; compatibility with the actual operator process is what excludes it.
Prescribed backgrounds are diagnostic inputs, not self-generated internal
resources or an already closed universe. No continuum Gauss matching,
Einstein equation, quantum gravity, or complete physical prediction claimed.
Notes, results, code and condition ledger primary-reviewed; no image checks.
'''
for name in ('research_note_620.md','joint_quantum_exchange_noise.py',
             'joint_quantum_exchange_noise_results.json','unified_physics_condition_ledger_620.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round620_drafts/final_review.txt',review)
print('620 verification/publication scripts and primary review prepared')
