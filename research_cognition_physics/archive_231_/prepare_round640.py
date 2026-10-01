"""Freeze and prepare publication of the original spatial-block audit."""
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
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

mapping={'joint_regional_read_scale':'joint_spatial_block_reference',
         '638':'639','639':'640','640':'641',
         '3154':'3157','3157':'3160','1226':'1229','1229':'1232',
         '2067':'2074','2074':'2081'}
verify=replace((HERE/'verify_round639.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original regional read, fixed support and geometry variation.',
                      'Verify actual spatial blocking, reference, temporal kernel and geometry.')
verify=verify.replace("['display_formulas']==15","['display_formulas']==14")
write('verify_round640.py',verify)
publish=replace((HERE/'publish_round639.py').read_text('utf8'),
                mapping|{'## 285.':'## 286.','## 190.':'## 191.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第640轮完成：** [一次真实空间合并的共同热态、时间核与几何来源]({p}research_note_640.md)'
         '原中性二次模相邻合并后，精确热态Hamiltonian不能代替二频过程；'
         '同一空间消元须保记忆核及几何归一化来源。'
         '三组、十四式通过，最新640／3160，1232份编号科学文件、2081份保护证据。'
         '[核验]({p}research_round_640_checks.json)、[条件账]({p}unified_physics_condition_ledger_640.md)。'
         '限树级二次分支；完整Gauss、跨图连续及GR仍开放。')
order=('**当前执行顺序（640后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[641原空间消元后的共同区域过程]({p}round641_drafts/STATUS.md)，'
       '核原记录、共同参考及来源能否共用有记忆过程；保留尺度正则性条件，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原区域读量与物理体积的统一预算','原空间合并的共同状态、时间核与来源')
publish=publish.replace('区域信息上界不替代跨图连续状态','二次空间合并不等于完整连续物理')
publish=publish.replace('原区域读量、物理体积与尺度统一的相对信息预算','一次真实空间合并的共同热态、时间核与几何来源')
write('publish_round640.py',publish)
write('round640_drafts/research_note_640_draft.md',(HERE/'research_note_640.md').read_text('utf8'))
review="""640 primary-agent review; no independent agent review.
Previous goal turn was progress: 639 published, 2074 protected hashes and
seven navigation files postverified. Active goal remains unchanged.
Read root/research README, direction, state, latest entry/results and related
608/612/613/617/625/567/558/574/576/599 history. No live Python worker seen.
Initial stale file guesses were absent; actual original mass implementation
was located and read. No existing file or process was overwritten/stopped.
Yeo-Shim 2412.02074 SectionII and Morinelli et al 2010.11121 Sections3.2/5.1
checked from primary text/PDF. Failed HTML paths not treated as evidence.
Numerics use both original neutral masses, original curved target metric and
potential Hessian, and original geodesic nearest-neighbor quadratic expansion.
No added pointer. hbar=1 is explicit; no original hbar=.7 numerical equality
claimed. Tree-level canonical quadratic quantization is a declared branch:
not full interacting Gauss Gibbs, nonlinear ordering correction or controlled
finite-temperature interaction approximation.
Real Haar average and its actual discarded partner form an invariant block
of the original N=8 spatial Laplacian. Embedding is transverse-zero in given3D.
Partial trace is mathematical coarse observation, not free physical erasure.
Mean-force coefficients exactly reconstruct both covariances; temperature
dependence and mixed limiting state rule out one beta-independent generator.
Two-pole positive moment determinant independently excludes one-mode quadratic
all-time reproduction even at fixed beta. Schur kernel retains both poles.
Geometric source uses fixed canonical coordinates and original psi=exp(g);
partition normalization is fixed from the true two-mode thermal trace.
Retaining its derivative restores source; no arbitrary vacuum subtraction.
Frequency determinant identity is not an unregularized path-integral product;
kinetic normalization and measure would also need tracking in that representation.
Three groups passed first run. Note writing patch first failed atomically,
then complete files were created; no frozen history was altered.
No full gauge spatial map, nonlinear multi-time record, continuum or GR claimed.
Next641 verifies one common retained process and scale regularity, not hardware.
"""
for name in ('research_note_640.md','joint_spatial_block_reference.py',
             'joint_spatial_block_reference_results.json','unified_physics_condition_ledger_640.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round640_drafts/final_review.txt',review)
print('640 science frozen; verifier and publisher prepared')
