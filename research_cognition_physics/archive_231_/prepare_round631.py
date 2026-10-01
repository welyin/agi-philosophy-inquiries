"""Freeze round 631 and prepare its evidence-preserving publication."""
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

mapping={'joint_continuum_source_spectrum':'joint_tensor_stress_spectrum',
         '629':'630','630':'631','631':'632',
         '3127':'3130','3130':'3133','1199':'1202','1202':'1205',
         '2001':'2008','2008':'2015'}
verify=replace((HERE/'verify_round630.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original full matter cuts, UV matching and joint metric sources.',
                      'Verify original tensor stress cut and common curvature response.')
verify=verify.replace("['display_formulas']==19","['display_formulas']==14")
write('verify_round631.py',verify)
publish=replace((HERE/'publish_round630.py').read_text('utf8'),
                mapping|{'## 276.':'## 277.','## 181.':'## 182.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第631轮完成：** [原整代物质的张量应力谱与共同曲率系数]({p}research_note_631.md)'
         '原同一连续真空的应力切割具有迹与五个无迹通道；'
         '张量正谱复现旧费米Weyl平方运行系数，质量和几何来源共用同一物质计数。'
         '三组、十四式通过，最新631／3133，1205份编号科学文件、2015份保护证据。'
         '[核验]({p}research_round_631_checks.json)、[条件账]({p}unified_physics_condition_ledger_631.md)。'
         '限平直非接触二点；背景自洽、接触项、图映射及GR仍开放。')
order=('**当前执行顺序（631后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[632共同态、背景方程与接触项]({p}round632_drafts/STATUS.md)，'
       '先核旧真空与同阶变分，再接入实时张量来源；不把诊断背景当自洽几何，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原物质谱与几何来源的共同连接','原物质张量谱与曲率系数的共同连接')
publish=publish.replace('共形响应不是一般度规或空间维数','应力自旋2不是引力子或空间生成')
publish=publish.replace('原物质谱、真实时间响应与几何来源的共同连接','原整代物质的张量应力谱与共同曲率系数')
write('publish_round631.py',publish)
write('round631_drafts/research_note_631_draft.md',(HERE/'research_note_631.md').read_text('utf8'))
review="""631 primary-agent review. No independent agent review.
630 was scientifically verified, published, and checked after publication.
Its science, full condition ledger, next entry, and navigation snapshots stay
frozen. No running process is claimed as a wait without a live handle.
Historical search found old curvature running and general source identities,
not this original full-mass tensor cut / common metric normalization test.
Mature Lorentz spectral decomposition was checked against Karateev 2012.08538
sections 3 and 5. Normalization here is independently derived by phase space,
not copied across different spectral conventions.
Hu/Verdaguer 0802.0658 section 3 was read for positive distributional stress
noise, state dependence and the input semiclassical background equation.
It is not invoked as a derivation of Einstein dynamics.
Same 630 one-generation 16-Weyl masses and half weights. Same continuum
canonical free-fermion vacuum and fixed background are inputs.
Original gauge/interaction/record structure has not silently been replaced
by a claim that this free sector is the entire unified theory.
Dirac stress vertex e_ij T_ij=(e k).alpha fixes source normalization.
Direct projector Gram integrals over 72 angular nodes reproduce the seven
mass/trace/shear channels, including cross entries and the null direction.
The angular integrand is polynomial of degree at most four; quadrature is
exact in principle at the chosen angular order, up to roundoff.
Trace transition equals -mass transition, consistent with H having no
positive-negative transition and the (-+++) convention.
rho_tensor=rho2 P2 + rho_mass P0/3 satisfies the cut Ward identity.
This excludes local contacts/tadpoles and is not a nonlinear Ward proof.
Five traceless channels are a timelike off-shell little-group representation,
not five massless graviton polarizations, new particles or a graviton pole.
In the all-zero-mass check the relative mass source vanishes but physical
absolute scalar perturbations need separate differentiation; no contradiction
with the finite Yukawa matrix at zero scalar background is asserted.
Stress spin2 has positive omega^4 tail, requiring three subtractions.
Finite source polynomials remain unconstrained here; when covariantly
completed their coefficients and contact terms cannot be chosen independently.
Heat-kernel fractions (-7,-8,5)/360 become -W^2/20+11E4/360.
The bare fermion log has negative Weyl-square coefficient; the counterterm
and old 553 beta convention are positive. Pi and Euclidean Hessian differ
by a minus sign, as declared in 630.
h_ij=2 epsilon e_ij gives integral W^2=2 integral epsilon''^2;
therefore Hessian factor four matches the spectral z^4 logarithm.
c_f=.4 is for this one generation only. Old three-generation scalar/vector
total c=293/120 is not copied into this calculation.
No lower metric polynomial is claimed fixed by cuts without contacts.
The sixth dispersive moment and bound are analytic beta integrals.
The result applies below the common strictly positive threshold; it is not
a bound on all finite local terms or general curved backgrounds.
Gaussian stochastic variables can reproduce two-point noise only. No full
Gaussian quantum process, actual apparatus or self-consistent geometry is
claimed. Next 632 must address common state/background/contact conditions.
Three test groups passed on first execution. Report has fourteen formulas.
"""
for name in ('research_note_631.md','joint_tensor_stress_spectrum.py',
             'joint_tensor_stress_spectrum_results.json','unified_physics_condition_ledger_631.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round631_drafts/final_review.txt',review)
print('631 verification and publication prepared; science frozen')

