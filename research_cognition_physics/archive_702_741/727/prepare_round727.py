"""Publish scientific artifacts without replacing preceding evidence."""
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


ledger=(HERE/'unified_physics_condition_ledger_726.md').read_text('utf8')
_,rest=ledger.split('\n',1)
ledger='# 联合条件总账：727原物质基态、共同来源与保留涨落\n'+rest
ledger=ledger.replace('接[725全账](unified_physics_condition_ledger_725.md)，回填[726报告](research_note_726.md)。[结果](joint_recorded_classical_wall_results.json)、[核验](research_round_726_checks.json)。',
                      '接[726全账](unified_physics_condition_ledger_726.md)，回填[727报告](research_note_727.md)。[结果](joint_matter_ground_source_results.json)、[核验](research_round_727_checks.json)。')
marker='## 当前共同对象及仍存在的分支';assert marker in ledger
ledger=ledger.replace(marker,'**727当前增量：** 原完整质量与三方向Weyl传播共同决定的条件基态可消去能源本身的方差，却改变经典能源与来源，并保留带间来源噪声。局部有隙及Gauss提升条件明示；同次实际记录在该条件下集中于新有效能源。下一项核宏观任务真正需要的时空涂抹、共同响应与变化参考。\n\n'+marker)
updates={
 'C03':'727同一字段记录与基态投影逐点对易；完整能源可集中而真实来源残差趋非零',
 'C15':'727原质量和三方向Weyl星图共同基态保全部CAR与协变；不是完整手征连续',
 'C19':'727孤立局部带和Gauss等变截面是额外条件；完整热态存在不保证这个共同基态带',
 'C20':'727紧管带导数及隙常数仅固定图；微观瞬时零噪声不是宏观任务的自动要求',
 'C22':'727参考变动必须保e、Berry/Born-Huang项及带间来源；压缩源零方差不等于真实源零方差'}
rows=ledger.splitlines();seen=[]
for i,row in enumerate(rows):
    for key,value in updates.items():
        if row.startswith('|'+key+' '):
            parts=row.split('|');parts[2]+='；'+value;rows[i]='|'.join(parts);seen.append(key)
assert set(seen)==set(updates)
ledger='\n'.join(rows)+'\n'
ledger=ledger.rsplit('## 本轮合并与下一项',1)[0]+'''## 727原物质参考、经典能源及真实几何来源

- 591谱微扰／量子度量、602—603热参考、630—633连续真空和720参考运输直接继承，不重报一般基态公式。
- 紧管局部孤立费米带需另核Gauss等变截面；谱投影协变不等于固定构形费米态已成为独立Gauss态。
- 原完整动能作用于随场参考，精确形式产生联络、量子度量和原基态能量e。此压缩不自动给真实绝热时间误差。
- 同一实际字段instrument保带；在声明的固定图半经典窗口，完整能源集中于Hb+e，旧受约束Einstein来源须重核。
- 原几何来源保带间矩阵元。其瞬时方差与参考投影的实际变化共同受谱公式约束，不可直接换成P G P。
- 原全部32模式质量与三方向128模式星图均核；独立中性Fock及实际y后态积分核系数，不冒充全图Gauss求解。
- 微观非零瞬时噪声不构成宏观几何或统一目标的否定；下一步须用实际任务的分辨率和响应量词。

## 本轮合并与下一项

参考选择、完整能源、原传播及几何来源在同一条件带中连接；新增局部带／Gauss条件与旧来源变更明确列账，未减少原物理参数或证明引力生成。

接[728](round728_drafts/STATUS.md)，回查591及630—633亚隙谱／涂抹，核同一原参考的宏观来源、真实有限时间任务与背景变化。停止占据、单带函数或波包精度优化；旧空间、649／699范围和统一目标保持。
'''
write('unified_physics_condition_ledger_727.md',ledger)
write('round727_drafts/research_note_727_draft.md',(HERE/'research_note_727.md').read_text('utf8'))
write('round728_drafts/STATUS.md','''# 第728轮入口：同一参考的宏观来源与真实时间尺度

接[727](../research_note_727.md)、[全账](../unified_physics_condition_ledger_727.md)。原共同费米基态能消去能源方差，但改变来源并保留带间瞬时噪声；不能把微观瞬时零噪声偷偷升级为统一模型必须满足的公理。

1. 复用591亚隙响应、602状态／响应区分、630—633共同谱与原记录、704来源导数；不重做一般低通滤波或Berry公式。
2. 明确实际宏观任务读的是哪一份时空涂抹来源、响应和后态；同时保原均值与接触项，不能用P G P改掉观测。
3. 检验固定参考的时间涂抹是否已由旧谱解决；若是，只作入口记录，不计新轮。
4. 真正缺口是参考／背景会变化且发生真实记录时，同一过程能否有可用的尺度层级、有限噪声和来源反馈。不得先假设每次重新回基态。
5. 626—627、633因果、649引力及699正性边界照旧；不把星图或自由分支当原完整相互作用连续模型，不展开装置工程。
''')
write('round727_drafts/literature_scope_audit.json',json.dumps(dict(
 external_sources_newly_imported=[dict(title='Panati, Spohn, Teufel: The time-dependent Born-Oppenheimer approximation',
 url='https://arxiv.org/html/0712.4369',verified_sections='1 and2; equations13-16; projection is not an invariant evolution',
 use='Berry and Born-Huang projection method, independently derived with original curved kinetic form',
 not_imported='Molecular mass hierarchy, automatic second-order time accuracy or uniform continuum gap')],
 inherited='558,574,591,602,603,604,630,632,633,720,726; physical-right CAR state dictionary retained.',
 own_mapping='Same actual record semiclassical ground-band energy center, retained source-noise limit and original three-direction nonuniform full-matter calibration.',
 excluded='Counting standard spectral identities as new; Gauss invariance inferred from covariance; projected source substituted for actual source; original classical source claimed unchanged.'
),ensure_ascii=False,indent=2)+'\n')
write('round727_drafts/scope_and_dedup_review.md','''# 727范围与去重审查

主代理审查，无新增独立代理。上一完成单元726属于实质进展；本轮接续其参考缺口。

- 591已有量子度量与激发来源矩阵元身份，632已有真空对背景的改变。此次新连接是726同源读后误差及原全部质量和方向边，不重新编号一般公式。
- 条件费米带不同于603完整Gauss热态或591全玻色基态。稳定子平凡／等变截面另列；星图只核协变，不声称完成全图Gauss带。
- χ的导数均在固定紧管控制；全图一致、固定Planck连续、真实绝热长时估计没有由此得到。
- 精确形式压缩不是时间不变子空间。原a包含电动能；数值只给径向收缩，不把它当全Born-Huang值。
- 能源中心含e，实际瞬时源含带间部分；既不删除负基态能量，也不通过压缩观测消去真实噪声。
- 三方向星图含全部128物理模式，保原非均匀字段；256是Nambu冗余维数。它不修复604倍增。
- 16维Fock与Nambu独立核对、实际y核积分及差分交叉检查已通过；没有新增模拟相互作用全图或现实物理预测。
- 下一步以宏观观测／过程为准，禁止把每个微观源必须零方差当未声明验收条件。
''')
files=['research_note_727.md','joint_matter_ground_source.py','joint_matter_ground_source_results.json',
       'unified_physics_condition_ledger_727.md']
write('round727_drafts/final_review.txt','727 primary review; no independent agent review.\n'+
      '\n'.join(n+' '+hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in files)+'\n'+
      'Three scientific check groups passed. Conditional Gauss lift and isolated band explicitly retained. Original source and actual off-band noise preserved. No full goal completion.\n')
verification=(HERE/'verify_round726.py').read_text('utf8')
verification=remap(verification,{
 'joint_recorded_classical_wall':'joint_matter_ground_source',
 'joint_relational_normal_form':'joint_recorded_classical_wall',
 'round727':'round728','round726':'round727','round725':'round726',
 '_726':'_727','_725':'_726',
 'range(584,726)':'range(584,727)','3158':'3172','3172':'3186',
 '==18':'==16','round=726':'round=727','round=725':'round=726',
 '3383':'3386','1490':'1493',
 'common_wall_limit_entry':'fermion_reference_entry'})
write('verify_round727.py',verification)
post=(HERE/'postcheck_round726.py').read_text('utf8')
post=remap(post,{'range(584,727)':'range(584,728)','round726':'round727','_726':'_727',
                '3172':'3186','第726':'第727','round=726':'round=727'})
write('postcheck_round727.py',post)
