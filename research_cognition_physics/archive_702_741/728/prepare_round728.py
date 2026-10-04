"""Preserve728 Gauss-band results and original nonflat source distinction."""
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


ledger=(HERE/'unified_physics_condition_ledger_727.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：728完整物质基态的Gauss提升与来源分支\n'+rest
ledger=ledger.replace('接[726全账](unified_physics_condition_ledger_726.md)，回填[727报告](research_note_727.md)。[结果](joint_matter_ground_source_results.json)、[核验](research_round_727_checks.json)。',
                      '接[727全账](unified_physics_condition_ledger_727.md)，回填[728报告](research_note_728.md)。[结果](joint_ground_gauss_lift_results.json)、[核验](research_round_728_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**728当前增量：** 原574非平坦来源的完整Fock稳定子字符障碍关闭，但实际谱隙待验；727平坦内部链路分支另有严格有理谱隙下界和完整一代电荷抵消，局部Gauss基态参考存在。两分支不混为原Einstein解；真实记录后态合法而不自动回基态。下一项恢复原非平坦链路共同核验。\n\n'+marker)
updates={
 'C02':'728条件基态管状截面的Gauss代表元独立性明确；不由谱投影协变直接推状态合法',
 'C15':'728原完整一代夸克轻子电磁字符抵消；使用已有表示，不重新推出规范群',
 'C19':'728原非平坦源的基态字符无阻碍；平坦链路分支给有理正隙与真实局部Gauss带',
 'C20':'728特定平坦分支质量隙有图无关下界，不等于参考导数或连续手征极限统一',
 'C22':'728源分支和动量变更显式列账；平坦参考不能继承旧非平坦Einstein初值'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,val in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+val;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 728参考态的Gauss合法性与来源分支

- 634记录条件方差、624—625／704来源运输直接复用，不以再次拆方差虚增研究。
- 574原稳定子与等变管直接复用。新增是598完整Fock基态线字符与商群中心的具体核验。
- 原非平坦源稳定子在物理表示中仅留连续颜色SU3的一维字符，因此任何单重基态线无额外字符障碍；其全图谱隙另验。
- 平坦内部链路分支的原方向边与均匀比较质量反对易；原非均匀径向质量为节点直和小扰动。有理外界给全图谱隙大于.013。
- 原完整一代t=0时夸克／轻子电荷相消，无隙关闭的路径保持零电磁字符；Gauss截面由真实Fock投影构造。
- 平坦参考可与新的零相位Gauss包组合，但明确改变初始源，不继承573／651的Einstein解；原非平坦条件分支保其原相位。
- 原sterile完整两结果后态仍为合法有限能源Gauss态，不保证留在基态带、重新热化或原约束来源不变。

## 本轮合并与下一项

C02／C15／C19共同参考的合法性取得正面连接；特定分支中谱隙与字符均不再是任意输入，但原非平坦全图谱隙、实际动态及几何自洽仍缺。

接[729](round729_drafts/STATUS.md)，恢复574原非平坦弱链路与几何来源，核同一物质参考；停止一般占据和时间窗优化，旧空间、649／699及统一目标保持。
'''
write('unified_physics_condition_ledger_728.md',ledger)
write('round728_drafts/research_note_728_draft.md',(HERE/'research_note_728.md').read_text('utf8'))
write('round729_drafts/STATUS.md','''# 第729轮入口：恢复原非平坦来源中的共同物质参考

接[728](../research_note_728.md)、[全账](../unified_physics_condition_ledger_728.md)。原574稳定子在完整Fock中无基态字符障碍；728平坦链路分支另有显式正隙，但不能直接代回非平坦原Einstein来源。

1. 先回查569—574原字段、弱／圆链路采样、正几何、原p及Gauss；后续604的Weyl边是已声明的额外动力分支。
2. 恢复同一原非平坦链路，审计实际物质谱与参考。不能用单位链路、随意平移质量或改Y获得谱隙后声称原源相同。
3. 区分有限图实际谱证据、可认证隙下界、所有网格的统一界。特定图无隙或变化不推出整个统一计划失败。
4. 基态若存在，原完整能源与几何来源依旧包含e及带间噪声；不能丢弃原非零动量、记录后态或真空反作用。
5. 634条件方差、591／602谱响应、623—625／704实际记录与来源直接复用。旧空间、649、699及统一目标保持。
''')
write('round728_drafts/literature_scope_audit.json',json.dumps(dict(
 external_sources_newly_imported=[],
 inherited='574 compact-group tube, exact nonflat stabilizer and Gauss phase;598 original Fock representation;604 declared Weyl edges;634 conditional record sources;727 isolated band and source accounting.',
 own_mapping='Fock stabilizer character, complete-generation electromagnetic cancellation and rational graph-independent gap in the declared flat-link branch; actual local Gauss band and record states.',
 elementary_tools='Finite Hermitian perturbation, compact group character and local spectral projector construction are derived explicitly, not claimed as new general theorems.',
 excluded='Repeating634 variance decomposition; using a flat-link gap for the nonflat Einstein source; treating projection covariance as a full Fock phase proof.'
),ensure_ascii=False,indent=2)+'\n')
write('round728_drafts/scope_and_dedup_review.md','''# 728范围与去重

主代理审查，无新增独立代理。上一727及728入口已取得实质进展。

- 634已经给记录间／内方差；624—625／704已有响应运输，本轮转向727实际缺失的Gauss截面条件。
- 574已有原非平坦来源稳定子，本轮只把其剩余中心作用核到598完整费米表示，未复写一轮稳定子。
- 均匀比较质量与原三方向边的反对易是604身份的图上应用；新增为原非均匀字段的有理外界及无谱流全路径。
- 原平坦链接分支、原非平坦Einstein来源分开列。后者的实际谱隙没有得到证明；前者的新零相位态不保原p或原约束解。
- 唯一Fock真空与简并单粒子谱不混淆；完整物种电磁字符相消，不把单粒子投影协变当字符平凡。
- Gauss截面用Fock投影及稳定子平凡构造；有限矩阵测试仅校准谱、荷和原表示，不声称显式存储2^128态。
- 实际sterile仪器保Gauss与有限能源，但会出基态带。原能源中心e、噪声及制备缺口保留。
- 有理gap界对图大小不变只属于声明的质量窗口／平坦内部链路分支，不是共同连续或SM手征修复。
''')
files=['research_note_728.md','joint_ground_gauss_lift.py','joint_ground_gauss_lift_results.json',
       'unified_physics_condition_ledger_728.md']
write('round728_drafts/final_review.txt','728 primary review; no independent agent review.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+'\n'+
      'Three checks passed. Rational interval evidence and analytic full-graph bound distinguished from star calibration. Original nonflat source not replaced. Goal remains open.\n')
v=(HERE/'verify_round727.py').read_text('utf8')
v=remap(v,{
 'joint_matter_ground_source':'joint_ground_gauss_lift',
 'joint_recorded_classical_wall':'joint_matter_ground_source',
 'round728':'round729','round727':'round728','round726':'round727',
 '_727':'_728','_726':'_727','range(584,727)':'range(584,728)',
 '3172':'3186','3186':'3200','==16':'==14',
 'round=727':'round=728','round=726':'round=727',
 '3386':'3389','1493':'1496','fermion_reference_entry':'macroscopic_source_entry'})
write('verify_round728.py',v)
post=(HERE/'postcheck_round727.py').read_text('utf8')
post=remap(post,{'range(584,728)':'range(584,729)','round727':'round728','_727':'_728',
                '3186':'3200','第727':'第728','round=727':'round=728'})
write('postcheck_round728.py',post)
pub=(HERE/'publish_round727.py').read_text('utf8')
pub=remap(pub,{'joint_matter_ground_source':'joint_ground_gauss_lift',
    'round728':'round729','round727':'round728','_727':'_728','第727':'第728','（727后':'（728后',
    '727／3386':'728／3389','726／3383':'727／3386','231—727':'231—728','231—726':'231—727',
    '3386':'3389','1493':'1496','3186':'3200','latest_round=727':'latest_round=728',
    'next_round=728':'next_round=729','第728轮研究索引':'第729轮研究索引',
    '## 373.':'## 374.','## 278.':'## 279.',
    '|727|':'|728|','原物质基态、同一经典来源与保留的量子涨落':'完整物质基态的Gauss提升与原来源分支'})
# Human-facing summaries are explicitly rewritten for the final result.
lines=pub.splitlines()
for i,line in enumerate(lines):
    if line.startswith('summary='):
        lines[i]="summary='**第728轮完成：** [完整物质基态的Gauss提升与原来源分支]({p}research_note_728.md)原非平坦源的基态字符障碍关闭；平坦内部链路分支另证全图正隙及完整一代电荷抵消，局部Gauss参考成立。三组、十四式通过，最新728／3389，1496份编号科学文件、3200份保护证据。[核验]({p}research_round_728_checks.json)、[全条件账]({p}unified_physics_condition_ledger_728.md)。两分支不混作同一Einstein来源。'"
    if line.startswith('order='):
        lines[i]="order='**当前执行顺序（728后，优先于下方历史安排）：** 接[729原非平坦来源与物质参考]({p}round729_drafts/STATUS.md)，恢复569—574实际规范链路、原字段及来源核费米谱；不以平坦证书替代原源。旧方差与响应结果复用，旧空间、649／699及统一目标保持。'"
pub='\n'.join(lines)+'\n'
write('publish_round728.py',pub)
