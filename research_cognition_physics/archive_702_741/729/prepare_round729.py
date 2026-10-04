"""Prepare729 artifacts; preserve frozen history and separate the continuum branches."""
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

ledger=(HERE/'unified_physics_condition_ledger_728.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：729原曲几何、完整周期物质参考与同一来源\n'+rest
ledger=ledger.replace(
 '接[727全账](unified_physics_condition_ledger_727.md)，回填[728报告](research_note_728.md)。[结果](joint_ground_gauss_lift_results.json)、[核验](research_round_728_checks.json)。',
 '接[728全账](unified_physics_condition_ledger_728.md)，回填[729报告](research_note_729.md)。[结果](joint_curved_periodic_reference_results.json)、[核验](research_round_729_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**729当前增量：** 原共形几何的规范半密度动能为psi^-1 D_A psi^-1；原非平坦链接、全部32模式与完整4³／6³周期图共同给孤立参考数值证据及代表矩阵残差下界。其局部几何来源与同一基态能源导数一致。当前有限图参考可接728稳定子提升；未认证精确连续输入，不是所有网格统一隙、手征调节器或自洽Einstein解。停止继续放大中心差分网格。\n\n'+marker)
updates={
 'C12':'729同一给定共形几何进入原费米动能与局部几何变分；不是生成该几何',
 'C19':'729原完整周期图的孤立参考与728 Gauss提升相容，限声明的有限代表背景',
 'C20':'729两网格非零谱及残差下界不外推统一连续；604倍增、699指定候选反例保持',
 'C22':'729同一局部源等于参考能源导数；e0及其来源仍须加入几何反馈，未完成自洽'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 729原几何、全部物种及真实来源

- 原空间旋量连接与规范CAR测度相容，得到C K0 C；此映射属于继承的三几何／spin输入。
- 574实际弱／圆链接及598完整质量保留；颜色分块是原对称性，Nambu半权不删物种。
- 4³／6³代表矩阵的最小谱约.03767／.03790，保守残差下界约.006167／.006095；未计输入椭圆和连续误差。
- 原x链接的非零弱混合保证574稳定子机制在两个网格适用；728字符分析和局部Gauss带构造直接继承。
- 实际局部共形源为-r与K的半反对易式；同一基态的能量差分与来源相符。
- 不能把冻结初值当成其后超静态时空，不能把有限正Hamiltonian当成已通过604／699的物理手征候选。
- 578—579、613的减除一致性及630—635的响应／局部反项已存在，不作为后续新发现重复编号。

## 本轮合并与下一项

本轮合并C12／C19／C22的有限原背景接口；C20物理连续、真实过程和完整几何反馈继续开放。

接[730](round730_drafts/STATUS.md)，审计原连续物质、可容许参考及同一几何来源的共同连接，复用633／663／667和630—635。停止更多中心差分谱隙扫描；不能把原非零几何动量删掉后宣称同一宇宙。旧空间、649／699范围和统一目标保持。
'''
write('unified_physics_condition_ledger_729.md',ledger)
write('round729_drafts/research_note_729_draft.md',(HERE/'research_note_729.md').read_text('utf8'))
write('round730_drafts/STATUS.md','''# 第730轮入口：共同连续物质参考与真实几何过程

接[729](../research_note_729.md)、[全账](../unified_physics_condition_ledger_729.md)。原有限周期背景的规范链接、曲几何动能、完整物质参考及其来源已经共同检查；尚不是物理连续或自洽引力结果。

1. 回查633的原连续自由参考与Hadamard／实际记录，663的原几何自旋主部，667的局部参考与传播字典；先检查能否共用原时空、物种和状态。
2. 578—579／613的真空减除一致性，630—635的连续响应、局部反项与熵来源直接复用；不把一般重整化提醒另立一轮。
3. 必须区分573非零几何动量的初值、729冻结的瞬时Hamiltonian以及另选的静态／超静态时空。不能通过偷偷改背景得到参考后称原自洽解。
4. 若采用成熟曲时空Dirac参考定理，逐项核原规范连接、Dirac／Majorana混合、初始态短距和实际局域记录合同；标明外加的全局双曲、光滑、spin及重整化输入。
5. 停止N更大网格、有限逆范数优化及局部谱扫描。604中心差分仍有倍增，699只否定已声明固定手征拼接候选，700／701共同极限条件保持。
6. 正式新轮需产生真实共同连接或限定反例；若仅重述旧结论，保入口去重记录，不记轮次。旧空间382—386、425、522—523、649及统一目标保持。
''')
write('round729_drafts/literature_scope_audit.json',json.dumps(dict(
 external_sources_newly_checked=[dict(title='Matthias Fischmann, On Conformal Powers of the Dirac Operator on Spin Manifolds',
 url='https://arxiv.org/html/1311.4182',scope='Mature conformal covariance and spin-connection conventions; not a new reconstruction theorem.')],
 inherited='573 positive geometry;574 original nonflat source and stabilizer;598 complete matter;604 finite Weyl edges and doubling obstruction;663 spin geometry;727 source;728 Gauss lift.',
 own_mapping='Original conformal half-density kinetic, complete periodic represented spectra with arithmetic-model inverse residual bounds, and same-ground local geometric source.',
 exclusions='No exact continuum coefficient certificate, uniform gap, one-generation chiral regulator, autonomous preparation, adiabatic theorem or self-consistent Einstein backreaction.'
),ensure_ascii=False,indent=2)+'\n')
write('round729_drafts/scope_and_dedup_review.md','''# 729范围及去重审查

主代理审查，无独立代理。三组实际检查通过，不进行图像检验。

- Dirac共形协变和Hellmann–Feynman为成熟工具，贡献是原几何／规范／质量及真实来源的共同映射。
- 原64维Gamma中的旋量连接与三维半密度直接逐项检查；不同于619四维框架权。
- 颜色分块保原全部32模式，轻子Majorana和Nambu半权完整。未以少物种取得参考。
- 逆残差下界明确限定代表binary64矩阵及标准误差模型；不是全程区间认证，未含原psi连续误差。
- 两个周期图实际参考是有限结论。604、699和原正有限H之间不串换证明。
- 574稳定子证明确认在N=4、6的原x链接仍非退化，再复用728基态字符。
- 同一来源固定phi、holonomy和坐标图；没有把约束解总导数当局部metric变分。
- 573原几何动量不因冻结Hamiltonian而消失；接超静态参考必须声明新分支。
- 已回查578—579、613、630—635；真空减除与压力身份不重复作新轮次。
- 旧空间382—386、425、522—523按667复用。停止更多有限网格谱隙优化，转共同连续参考的实际前提。
''')
files=['research_note_729.md','joint_curved_periodic_reference.py','joint_curved_periodic_reference_results.json',
       'unified_physics_condition_ledger_729.md']
write('round729_drafts/final_review.txt','729 primary review; no independent agent review.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+'\n'+
      'Three groups passed; conformal identity and original source jointly retained. Represented finite matrices distinguished from exact continuum data. Goal remains open.\n')
v=(HERE/'verify_round728.py').read_text('utf8')
v=remap(v,{
 'joint_ground_gauss_lift':'joint_curved_periodic_reference',
 'joint_matter_ground_source':'joint_ground_gauss_lift',
 'round729':'round730','round728':'round729','round727':'round728',
 '_728':'_729','_727':'_728','range(584,728)':'range(584,729)',
 '3186':'3200','3200':'3214','==14':'==16',
 'round=728':'round=729','round=727':'round=728',
 '3389':'3392','1496':'1499','macroscopic_source_entry':'nonflat_reference_entry'})
v=v.replace('Reproduce726','Reproduce729')
write('verify_round729.py',v)
post=(HERE/'postcheck_round728.py').read_text('utf8')
post=remap(post,{'range(584,729)':'range(584,730)','round728':'round729','_728':'_729',
                '3200':'3214','第728':'第729','round=728':'round=729'})
post=post.replace('Check726','Check729')
write('postcheck_round729.py',post)
pub=(HERE/'publish_round728.py').read_text('utf8')
pub=remap(pub,{'729':'730','728':'729','727':'728','3389':'3392','3386':'3389',
              '1496':'1499','3200':'3214','## 374.':'## 375.','## 279.':'## 280.',
              'joint_ground_gauss_lift':'joint_curved_periodic_reference',
              '完整物质基态的Gauss提升与原来源分支':'原曲几何、完整周期物质参考与同一来源',
              '完整物质基态的Gauss提升与来源分支':'原曲几何、完整周期物质参考与同一来源',
              '旧空间合同保持，平坦带隙不替代原非平坦来源':'旧空间合同保持，有限曲参考不替代共同连续过程'})
lines=pub.splitlines()
for i,line in enumerate(lines):
    if line.startswith('summary='):
        lines[i]="summary='**第729轮完成：** [原曲几何、完整周期物质参考与同一来源]({p}research_note_729.md)原共形旋量动能与完整周期物质参考接通；4³／6³实际非平坦链接均有孤立谱证据及代表矩阵残差界，局部几何源与同态能源差分一致。三组、十六式通过，最新729／3392，1499份编号科学文件、3214份保护证据。[核验]({p}research_round_729_checks.json)、[全条件账]({p}unified_physics_condition_ledger_729.md)。未完成物理连续、手征调节器或自洽Einstein反馈。'"
    if line.startswith('order='):
        lines[i]="order='**当前执行顺序（729后，优先于下方历史安排）：** 接[730共同连续参考与真实几何过程]({p}round730_drafts/STATUS.md)，复用633／663／667及630—635，核原物质、允许参考和来源能否共用。停止更多中心差分谱隙扫描；不抹掉原非零几何动量。旧空间、649／699范围及统一目标保持。'"
write('publish_round729.py','\n'.join(lines)+'\n')

