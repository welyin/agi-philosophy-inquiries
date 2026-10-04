"""Prepare675; preserve history, original scope and active goal."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_674.md').read_text('utf8')
ledger=ledger.replace('# 联合条件总账：674互补配对与完整标量泛函实性','# 联合条件总账：675实际标量边缘、辅助平均与共同正性')
ledger=ledger.replace('接[673全账](unified_physics_condition_ledger_673.md)，回填[674报告](research_note_674.md)。[结果](joint_gauss_reflection_reality_results.json)、[核验](research_round_674_checks.json)。',
    '接[674全账](unified_physics_condition_ledger_674.md)，回填[675报告](research_note_675.md)。[结果](joint_physical_boundary_average_results.json)、[核验](research_round_675_checks.json)。')
ledger=ledger.replace('674完整标量核Hermitian与总权重实性|',
    '674完整标量核Hermitian与总权重实性；675实际标量边缘、精确辅助多项式及大τ必要条件|')
ledger=ledger.replace('674原标量／质量来源反射已证；原物理态',
    '674原标量／质量来源反射已证；675实际标量边缘保留全辅助积分；原物理态')
ledger=ledger.replace('其与643原Gauss有序过程仍须匹配',
    '675固定λ极限不是已识别的原物理时间；其与643原Gauss有序过程仍须匹配')
ledger=ledger.replace('674标量核Hermitian已证；群选择',
    '674标量核Hermitian已证；675精确辅助平均后仍保留完整双Haar；群选择')
ledger=ledger.replace('任意费米来源反射、动态正性、指定归一与连续局域性仍缺',
    '675非平坦投影上的全S⁹有限多项式与秩选择规则已接；任意费米来源反射、动态正性、指定归一与连续局域性仍缺')
ledger=ledger.replace('674总权重及反射配对实质量来源为实；',
    '674总权重及反射配对实质量来源为实；675复用591旧基态给当前核受控极限，未指定实际宇宙参考；')
at=ledger.index('## 本轮合并与下一项')
ledger=ledger[:at]+'''- 675保留全部辅助S⁹后定义实际标量边缘。固定E、部分子式或有限样本的符号不能替代这一对象的正性。
- 657球面积分多项式直接拉回非平坦u，得到原质量因子乘精确辅助最高次系数；r_u非4倍数的纯物理来源边缘为零。没有宣布所有这些谱秩均由真实Wilson背景实现。
- 591唯一正Gauss基态与谱隙直接复用；675新增的是它到当前双边界标量核的显式指数误差界。固定λ的大τ数学族尚未与原H_F时间等同。

## 本轮合并与下一项

675连接C01、C04、C14、C16、C17、C19：同一非平坦规范背景中的原质量、物理来源、球面平均及双Gauss边界有明确共同对象。精确消去E的有限代数式保留相位；可把当前候选正性约束转成完整群平均B的可认证必要条件，而非仅看固定辅助背景。

两组验证包含完整512维带相位分解，以及原128维非平坦手征空间中的4／8模式精确矩交叉核验。后者只是代码接口检查，不替代完整最高次收缩。先导抽样仍未确定B符号，未声称正性、非零归一或原H_F同一性。

四个分支仍未统一。不给固定图谱隙新增细化一致性，不把基态参考制备、微观时间与空间维数自动补齐；原连续映射、量子约束、普适动态几何仍缺。

接[676](round676_drafts/STATUS.md)：回到643原有序CAR过程与659边界字典，核时间、实际观测及共同来源的具体身份，复用653—655／657／667已完成限制。精确球面式留作接口，避免以张量枚举或采样优化代替统一目标。旧空间合同及整体研究目标不变。
'''
write('unified_physics_condition_ledger_675.md',ledger)
mapping={'joint_gauss_reflection_reality':'joint_physical_boundary_average',
         '673':'674','674':'675','675':'676','3253':'3255','3255':'3257',
         '1331':'1334','1334':'1337','2419':'2433','2433':'2444'}
verify=replace((HERE/'verify_round674.py').read_text('utf8'),mapping)
verify=verify.replace('Verify scalar reflection via complementary Pfaffians and full Gauss reality.',
    'Verify physical marginal, original nonflat auxiliary contraction and ground-kernel scope.')
verify=verify.replace('import joint_gauss_boundary_functional as prior_model','import joint_gauss_reflection_reality as prior_model')
verify=verify.replace('unnormalized_source_probe','ground_haar_probe')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_675.md','round675_drafts/research_note_675_draft.md',
           'round675_drafts/final_review.txt','round676_drafts/STATUS.md',
           'round675_drafts/ground_haar_entry.md','round675_drafts/ground_haar_probe.py',
           'round675_drafts/ground_haar_probe_results.json','round675_drafts/entry_checks.json')
"""+verify[b:]
verify=verify.replace("assert text['display_formulas']==16","assert text['display_formulas']==14")
write('verify_round675.py',verify)
publish=replace((HERE/'publish_round674.py').read_text('utf8'),mapping|{'## 320.':'## 321.','## 225.':'## 226.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+'''summary=('**第675轮完成：** [物理标量边缘、完整辅助平均与共同正性判据]({p}research_note_675.md)'
         '657精确球面多项式接入原非平坦投影与质量；实际标量边缘保留完整双Gauss平均。'
         '591旧基态给此候选的受控大τ必要正性条件；固定λ尚非原物理时间。'
         '完整积分符号、RP及原过程身份仍开放。'
         '两组、十四式通过，最新675／3257，1337份编号科学文件、2444份保护证据。'
         '[核验]({p}research_round_675_checks.json)、[全条件账]({p}unified_physics_condition_ledger_675.md)。'
         '连续及量子GR未完成。')
order=('**当前执行顺序（675后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[676完整平均后的共同物理过程]({p}round676_drafts/STATUS.md)，'
       '回查原有序过程、实际观测与时间／来源字典，不继续扩张同类抽样；目标不改。')
'''+publish[b:]
publish=publish.replace('互补手征配对、完整规范边界与标量泛函的实性','物理标量边缘、完整辅助平均与共同正性判据')
publish=publish.replace('原配对、完整边界与标量反射的共同合同','实际标量边缘与旧基态的共同正性限制')
publish=publish.replace('完整标量泛函的实性与正性的区别','非平坦辅助收缩与实际共同时间的边界')
write('publish_round675.py',publish)
write('postcheck_round675.py',replace((HERE/'postcheck_round674.py').read_text('utf8'),mapping))
write('round675_drafts/research_note_675_draft.md',(HERE/'research_note_675.md').read_text('utf8'))
review='''675 primary-agent proof/code/scope review; no independent agent.
Previous turn completed674, protected evidence and executed675 entry: progress.
Current README/direction/state,674 note/result,675 pilot/result/receipt read.
No running Python found. Goal unchanged. No task/automation/image/install.
382-386,425,522-523 remain inherited.591 ground-state existence not new.
657 sphere formula reused; nonlinear local polynomial pulled back to ux*b.
Each site still has only64 underlying Grassmann linear forms, hence k<=16.
Coefficient formula valid for all nonzero Wilson ranks; odd or2mod4 ru
give zero physical-source edge after auxiliary averaging. No claim that all
algebraic ranks occur for actual group or topological backgrounds.
Physical factor independent of E and inverse-free, all16 original channels.
Compression norm<=1, rectangular L norm<=1, NW bound1+|lambda|*normP;
Pf power(rv+r)/2 valid when even, odd sector vanishes.
Physical scalar marginal integrates both E fields; fixed-E sign not a no-go.
591 gap applied after two heat smoothings; A_i finite, gauge invariant.
Remainder norm exp[-(tau-2sigma)Delta], products give stated error.
M finite uniformly over full compact gauge and S9, locally in q.
673 Haar-null Wilson zeros and domination give continuity at each q.
Strict negative finite Gram transfers to invariant wavepackets by continuity.
B negative only falsifies fixed-lambda all-tau candidate; positive B not RP.
No actual original HF time, low-temperature state, infinite-volume gap claim.
Full finite factor checked; exact eight/four minors use actual spectral u,
not invented reduced physical theory. Whole128-variable coefficient unevaluated.
Independent cubature separates degree4 and2+2; odd-site monomials average0.
Old pilot mean/sign unresolved, large contribution concentration preserved.
General Spin9 sufficient result not applied to original full gauge group.
2 groups,14 equations. Next676 returns to original process/time/source mapping,
not another sphere Monte Carlo batch or cognitive-control design.
'''
for name in ('research_note_675.md','joint_physical_boundary_average.py','joint_physical_boundary_average_results.json','unified_physics_condition_ledger_675.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round675_drafts/final_review.txt',review)
print('675 preparation completed')

