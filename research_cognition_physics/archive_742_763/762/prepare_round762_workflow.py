"""Prepare guarded 762 navigation publication from the previous frozen workflow."""
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
mapping={'762':'763','761':'762','760':'761','3470':'3473','3467':'3470',
         '1595':'1598','3637':'3653','3627':'3637','407':'408','312':'313',
         'joint_loop_scale_transport':'joint_gauss_first_order_process'}
def remap(text):
    return re.sub('|'.join(re.escape(k) for k in sorted(mapping,key=len,reverse=True)),
                  lambda m:mapping[m[0]],text)
pub=remap((HERE/'publish_round761.py').read_text('utf8'))
summary='**第762轮完成（共同首阶展开）：** [Gauss背景涨落与条件来源]({p}research_note_762.md)在761尺度分支内，同一Gauss准备的曲测度、宽度及量子荷与全部矩阵来源、原有限记录共同给首阶弱展开。固定图与时间，不是连续或量子GR。三组、十五式通过，最新762／3473，1598份编号科学文件、3653份保护证据。[核验]({p}research_round_762_checks.json)、[条件账]({p}unified_physics_condition_ledger_762.md)。'
order='**当前执行顺序（762后，优先于下方历史安排）：** 接[763同阶资料与受约束响应]({p}round763_drafts/STATUS.md)，先核内禀背景与费米诱导涨落的不同阶次、共同初值及完整约束；不继续局部积分精度。[范围审计]({p}round762_drafts/scope_and_dedup_review.md)。旧空间、604、649／699及总目标保持。'
pub=re.sub(r'^summary=.*$',lambda m:'summary='+repr(summary),pub,flags=re.M)
pub=re.sub(r'^order=.*$',lambda m:'order='+repr(order),pub,flags=re.M)
pub=pub.replace('第762轮完成（新增缩放的条件结果）','第762轮完成（共同首阶展开）')
pub=pub.replace('有限时间输运与同源反作用','Gauss背景涨落与条件来源')
pub=pub.replace('旧空间合同保持，新增阶次及同源响应','旧空间合同保持，共同涨落及来源')
pub=pub.replace('Publish762','Publish762')
post=remap((HERE/'postcheck_round761.py').read_text('utf8'))
post=post.replace('第762轮完成（新增缩放的条件结果）','第762轮完成（共同首阶展开）')
for name,text in (('publish_round762.py',pub),('postcheck_round762.py',post)):
    with (HERE/name).open('x',encoding='utf8',newline='\n') as stream:stream.write(text)
print('Prepared guarded round 762 publication.')

