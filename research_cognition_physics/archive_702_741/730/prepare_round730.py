"""Prepare730 continuous-reference bridge and preserve all prior evidence."""
import hashlib
import json
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent

def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)

def remap(text,mapping):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m.group()],text)

ledger=(HERE/'unified_physics_condition_ledger_729.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：730原动态连续参考、记录与同一来源\n'+rest
ledger=ledger.replace(
 '接[728全账](unified_physics_condition_ledger_728.md)，回填[729报告](research_note_729.md)。[结果](joint_curved_periodic_reference_results.json)、[核验](research_round_729_checks.json)。',
 '接[729全账](unified_physics_condition_ledger_729.md)，回填[730报告](research_note_730.md)。[结果](joint_dynamic_continuum_reference_results.json)、[核验](research_round_730_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**730当前增量：** 在573原给定光滑动态背景上，完整矩阵Dirac／Majorana及规范连接可接辅助过去的Hadamard构造；原目标邻域和非零时间jet完全保留。原光滑sterile记录的后态及相对来源共同正则，无需未来统一瞬时基态隙。背景协变态族不是全量子Gauss态，辅助过去不是内部制备；原图映射、绝对基准与自洽引力仍开放。\n\n'+marker)
updates={
 'C01':'730原完整矩阵物质的自对偶CAR与纯Hadamard参考族接通，限给定经典动态背景',
 'C03':'730原光滑sterile记录在变化质量／曲背景保同一短距结构与相对来源；因果装置仍缺',
 'C04':'730真实Cauchy运输保原质量／规范／非零几何jet，不以每刻负谱投影替代',
 'C19':'730可容许连续参考不要求未来统一瞬时隙，辅助过去选择不唯一且非物理制备',
 'C20':'730连续构造不等于原有限Gauss态的细化极限，604及699边界保持',
 'C22':'730相同减除给光滑记录来源差；变化参考的准备导数须与直接几何源共同保留'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.replace('|指定连续物质|630—635／645—646谱、来源和熵；663旋量字典；664联络匹配|',
 '|指定连续物质|630—635／645—646谱、来源和熵；663旋量字典；664联络匹配；730原给定动态背景的Hadamard参考族、自对偶CAR、光滑记录及同背景相对来源|')
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 730原动态连续参考的共同实现

- 原573短时发展中的度量、外曲率、全五标量及规范电场对应速度保留。数值只核原一阶jet，没有数值解整个未来。
- 原完整Gamma和质量允许固定Dirac模块识别；变化内部连接／质量为正确伴随的矩阵零阶项。Nambu自对偶性明确保留。
- 辅助插值只在目标区外，紧T³、正度量凸组合、正F球和固定spin结构给合法背景。它不是原Einstein方程的整个宇宙过去。
- 从原常量质量过去输送得到纯Hadamard态族，未来无需始终有瞬时谱隙；不同辅助处方不选择唯一宇宙态。
- 原光滑sterile仪器保短距奇性，全部质量及应力的同背景状态差有限；绝对应力和有限反项没有被固定。
- 背景规范协变不等于全量子Gauss投影。已知634来源跳接反例仍保持；真实因果装置和内生准备未完成。
- 真实参考随背景变化时保准备响应，与704合同相同；原完整局部矩阵诊断验证非零贡献，不另拟合压力。

## 本轮合并与下一项

C01／C03／C04／C19／C22在原给定连续背景的线性物质范围内共同成立；C20原图映射、全Gauss和动态反馈仍未由此关闭。

接[731](round731_drafts/STATUS.md)：回查原初始Gauss／Einstein约束与同态来源的具体反作用连接；复用325、572—574、578—583、630—635及649—651。停止瞬时谱、一般Hadamard和辅助过去参数优化。先明确相对源及绝对基准、选择性记录及闭合补偿，不把同背景状态差自动当自洽解。旧空间及统一目标保持。
'''
write('unified_physics_condition_ledger_730.md',ledger)
write('round730_drafts/research_note_730_draft.md',(HERE/'research_note_730.md').read_text('utf8'))
write('round731_drafts/STATUS.md','''# 第731轮入口：同一量子来源与原初始约束的反作用

接[730](../research_note_730.md)、[全账](../unified_physics_condition_ledger_730.md)。原给定动态背景已有完整线性费米参考族、光滑记录和相对来源；尚未证明量子来源加入后保原Gauss／Einstein约束。

1. 先回查325、572—574、578—583、630—635、649—651、727—730；现有Hadamard定理、参考导数、来源跳接及Bianchi提醒不再另计轮次。
2. 把原CMC初始约束、规范Gauss、标量／规范动量和同一量子来源的变化列在同一方程中。先检验核空间与总荷相容，不能只改正能量密度。
3. 572原两径向反流的Gram只有秩2，曾明确z方向不在其范围。634记录已有非零z动量；若此旧接法再次出现，应直接合并旧证据，禁止重新扫描或把它扩大成所有反馈不可能。
4. 同背景Hadamard态差可避免重复UV减除，但绝对参考的应力、有限反项及初始几何仍須共同给出。两套分别自洽的参考不能无映射相加。
5. 背景规范协变不是全Gauss物理态，记录条件化不是自动物理坍缩。可检验替代包括原物质其余动量／规范电流的补偿、允许几何自由资料的约束校正，以及完整有限时间守恒操作。
6. 先整合并验证真实共同接口，暂不转入具体认知装置设计。旧空间、604、649／699限定结论和统一目标保持。
''')
write('round730_drafts/literature_scope_audit.json',json.dumps(dict(
 sources=[
 dict(url='https://arxiv.org/html/2108.11630',use='Definition2.3, matrix formal adjoint, deformation and Theorem7.2; finite-multiplicity original module and self-dual reality checked explicitly.',
      caution='Fix original spin structure; do not infer uniqueness from orientation. Prove global hyperbolicity of our interpolation using compact T3; no arbitrary noncompact convex-gluing claim.'),
 dict(url='https://arxiv.org/html/1210.4031',use='Proposition3.2 external gauge/scalar Yukawa deformation and local Wick context.',
      caution='Scalar Yukawa theorem not directly claimed for all chiral Majorana couplings; no importing its constant-mass conserved stress result to general dynamic sources.'),
 dict(url='https://arxiv.org/html/0911.1304',use='Inherited Proposition4.12 context for Hadamard operation stability; actual finite-polynomial smooth-difference argument given.')],
 own_connection='573 dynamical initial collar, original full matrix matter/reality, state transport and same sterile record/source; no future uniform ground gap required.',
 excluded='Full quantum Gauss, lattice continuum, physical detector preparation, unique cosmological state or semiclassical Einstein solution.'
),ensure_ascii=False,indent=2)+'\n')
write('round730_drafts/scope_and_dedup_review.md','''# 730范围与去重审查

主代理审查，无独立代理；三组检查，无图像检验。

- 573原经典短时发展直接复用；辅助过去在目标邻域之外，原K和所有场速度不冻结。
- 原未来精确解未在代码中数值重算。一阶jet曲线只核字典，解析拼接用原实际光滑发展。
- 一般Hadamard存在是成熟定理；新连接是全部原矩阵物质、Nambu现实条件、非零jet及原记录来源。
- 原平凡束的有限Clifford多重性明确；不以Zahn标量Yukawa文献直接覆盖整个手征理论，不抄录定向便有唯一spin结构的错误说法。
- 紧T3保证所用平滑正度量插值全局双曲；该证明不推广任意非紧空间。
- 辅助过去不是物理制备，也不决定唯一态。未来不需瞬时基态隙；禁止继续用谱扫描替代动态参考。
- Cauchy输送态、Hadamard短距、同背景来源差及绝对应力四者分开验收。
- 有限±k矩阵包含全部原物种，但不是空间连续演化或Hadamard数值证明。差分只核实际源和准备导数。
- 633分布instrument的因果实现限制、634跳接来源、704准备响应、旧空间及699范围直接继承。
- 下一项应核原初始约束的全来源相容，不重做已知一般能量账或真空反项。
''')
files=['research_note_730.md','joint_dynamic_continuum_reference.py','joint_dynamic_continuum_reference_results.json',
       'unified_physics_condition_ledger_730.md']
write('round730_drafts/final_review.txt','730 primary review; no independent agent review.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+'\n'+
      'Three checks pass. Continuum proof uses stated mature tools and original matrix mapping; numerical evidence only validates initial jets and finite local source bookkeeping. Full goal remains open.\n')
v=(HERE/'verify_round729.py').read_text('utf8')
v=remap(v,{
 'joint_curved_periodic_reference':'joint_dynamic_continuum_reference',
 'joint_ground_gauss_lift':'joint_curved_periodic_reference',
 'round730':'round731','round729':'round730','round728':'round729',
 '_729':'_730','_728':'_729','range(584,729)':'range(584,730)',
 '3200':'3214','3214':'3228','==16':'==18',
 'round=729':'round=730','round=728':'round=729',
 '3392':'3395','1499':'1502','nonflat_reference_entry':'varying_reference_entry'})
v=v.replace('Reproduce729','Reproduce730')
write('verify_round730.py',v)
post=(HERE/'postcheck_round729.py').read_text('utf8')
post=remap(post,{'range(584,730)':'range(584,731)','round729':'round730','_729':'_730',
                '3214':'3228','第729':'第730','round=729':'round=730'})
post=post.replace('Check729','Check730')
write('postcheck_round730.py',post)
pub=(HERE/'publish_round729.py').read_text('utf8')
pub=remap(pub,{'730':'731','729':'730','728':'729','3392':'3395','3389':'3392',
              '1499':'1502','3214':'3228','## 375.':'## 376.','## 280.':'## 281.',
              'joint_curved_periodic_reference':'joint_dynamic_continuum_reference',
              '原曲几何、完整周期物质参考与同一来源':'原动态背景上的连续参考、真实记录与同一来源',
              '旧空间合同保持，有限曲参考不替代共同连续过程':'旧空间合同保持，背景连续参考不替代自洽引力'})
lines=pub.splitlines()
for i,line in enumerate(lines):
    if line.startswith('summary='):
        lines[i]="summary='**第730轮完成：** [原动态背景上的连续参考、真实记录与同一来源]({p}research_note_730.md)保原非零几何／物质时间jet，完整Dirac／Majorana可接Hadamard参考族及同一光滑记录来源差；无需未来统一瞬时谱隙。三组、十八式通过，最新730／3395，1502份编号科学文件、3228份保护证据。[核验]({p}research_round_730_checks.json)、[全条件账]({p}unified_physics_condition_ledger_730.md)。原图映射、全Gauss及自洽引力仍开放。'"
    if line.startswith('order='):
        lines[i]="order='**当前执行顺序（730后，优先于下方历史安排）：** 接[731同一量子来源与原初始约束]({p}round731_drafts/STATUS.md)，回查325、572—574、578—583、630—635及649—651，核全部来源、核空间和反作用相容。停止谱隙及辅助过去参数优化；旧空间、604、649／699及统一目标保持。'"
write('publish_round730.py','\n'.join(lines)+'\n')

