"""Freeze 629 and prepare history-preserving publication."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    patterns=[]
    for key in sorted(mapping,key=len,reverse=True):
        pattern=re.escape(key)
        if key.isdecimal():pattern=r'(?<!\d)'+pattern+r'(?!\d)'
        patterns.append(pattern)
    return re.sub('|'.join(patterns),lambda m:mapping[m.group()],text)


def write(name,text):
    path=HERE/name;path.parent.mkdir(exist_ok=True)
    with path.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_global_anomaly_bundle':'joint_topological_mass_phases',
         '627':'628','628':'629','629':'630',
         '3123':'3125','3125':'3127','1193':'1196','1196':'1199',
         '1987':'1994','1994':'2001'}
verify=replace((HERE/'verify_round628.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original global anomaly and non-liftable bundle interface.',
                      'Verify joint original mass phases and quotient-group topological weights.')
write('verify_round629.py',verify)
publish=replace((HERE/'publish_round628.py').read_text('utf8'),
                mapping|{'## 274.':'## 275.','## 179.':'## 180.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第629轮完成：** [原质量相位、拓扑权重与共同量子测度的压缩]({p}research_note_629.md)'
         '原群的四个拓扑角与原单代五个非零质量相位共用指数Jacobian；'
         '对常量模块重定相取商后剩三个角不变量，单独擦掉质量相位会改变它们。'
         '两组、十四式通过，最新629／3127，1199份编号科学文件、2001份保护证据。'
         '[核验]({p}research_round_629_checks.json)、[条件账]({p}unified_physics_condition_ledger_629.md)。'
         '仅压缩明示参数族；角值、连续区域过程、尺度及GR仍开放。')
order=('**当前执行顺序（629后，优先于下方历史安排）：** 继续整合条件，认知设计后置。'
       '接[630完整过程与连续物理的共同接口]({p}round630_drafts/STATUS.md)，'
       '回填623—629并核同一尺度的过程、测度及几何来源；不继续相位优化，目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('原全局群、物质与spin反常的合并','原质量、拓扑权重与参数冗余')
publish=publish.replace('反常工具的辅助维数不是时空生成','相位参数商不是空间维数')
publish=publish.replace('原商群、旋量几何与全局反常的共同条件','原质量相位、拓扑权重与共同量子测度的压缩')
write('publish_round629.py',publish)
write('round629_drafts/research_note_629_draft.md',(HERE/'research_note_629.md').read_text('utf8'))
review='''629 primary-agent review. No independent agent review.
628 was verified and published before 629. Its reports/results/next entry stay frozen.
Historical search distinguishes old Euler terms, phase-sensitive records and
the 615 holonomy parameter from the present global topological theta weights.
No earlier joint calculation of old five mass phases and four theta angles was found.
User priority remains condition integration, not cognitive implementation.
DGL Omega4^Spin(BG)=Z^4 and Garcia-Etxebarria/Montero nonbounding phase discussion
were read, together with Fileviez Perez/Patel's original mass/weak-angle paper.
The latter omits hypercharge in its simplified illustration; this calculation
retains the original quotient bundle and full gravitational index.
Tong's global-group discussion motivates using integral characteristic data,
not independent naive SU3/SU2/U1 instanton numbers on nonliftable bundles.
u=a^2/2, v=c2(E3), w=c2(E2), t=p1/48 are integer on the declared smooth closed
spin background class. Spin divisibility and K3 are mature theorem inputs.
Four geometric witnesses give characteristic matrix diag(1,1,1,-1).
Together with torsion-free rank-four bordism this proves the coordinate basis;
an arbitrary integer test matrix alone is not a geometric existence proof.
Original Q,uc,dc,L,ec,nuc bundles give the exact six-by-four index matrix A.
Weights are obtained using the frozen original J and polynomial Chern roots;
spin and Nambu multiplicities are not included a second time.
The exact matrix sums match 628's full-generation Chern character/index.
Original holomorphic annihilation coefficients are extracted from old
particle-hole BdG creation blocks and conjugated, not guessed from RH labels:
(-Yu*,-Yd*,-Ye*,-Ynu*,-Ys), after removing positive scalar factors.
Passive field convention is psi=exp(i alpha) psi_prime; couplings shift by C alpha
and the declared index Jacobian shifts theta by -A^T alpha. All signs consistent.
Six constant module rephasings act with rank six on the nine phase coordinates.
The three invariant integer rows have an identity theta subblock, so their
kernel is connected. The compact connected rank-six image is exactly that
kernel; this proves the full stated torus quotient, not only a tangent count.
Baryon vector (1,-1,-1,0,0,0) leaves every original coupling invariant while
shifting only theta_w by 3 beta in this integer characteristic basis.
Majorana lepton-number violation does not break this baryon redefinition.
Extra Higgs phase is already in the same image: A^T q=0 and
C(-q)=(3,-3,-3,3,0); this exact check was added and saved results reproduced.
No additional physical gauge symmetry or conservation law is asserted.
Finite CAR rephasing is a mass/spectrum dictionary check. Continuum anomalous
Jacobian is an explicit mature input, not derived from a finite matrix.
SU5 c2 direction (-2,1,1) is a restriction in a fixed mass convention, not
forced by the use of SU5 in the preceding anomaly proof. No claim is made
that it is an invariant constraint without simultaneously specifying masses.
The three phase invariants classify only this one-generation nonzero-coupling
nine-angle family under constant rephasings. No CKM/all-QFT/CP-duality count.
No phase values or experimental observables are predicted. Theta=0 is only
a coordinate example. Nonzero determinant amplitudes are not assumed on 628
zero-mode backgrounds; insertions/allowed sectors still require a real process.
Closed fixed-sector variation of topological numbers vanishes; boundary data
and full quantum sector sums remain open, not replaced by local-source checks.
Two groups passed first run; exact Higgs-scope check subsequently also passed,
with identical saved results. No images/new agents/new application tasks.
Next work returns to the shared quantum-process/continuum/geometry interface;
no more phase optimization or bundle enumeration is scheduled.
'''
for name in ('research_note_629.md','joint_topological_mass_phases.py',
             'joint_topological_mass_phases_results.json','unified_physics_condition_ledger_629.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round629_drafts/final_review.txt',review)
print('629 draft, review, verifier and publisher prepared')

