"""Preserve731 counterflow entry; formal730 count stays fixed."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round730_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'730':'731','729':'730','3392':'3395','3214':'3228',
         'varying_reference_entry':'gauge_counterflow_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**731入口已执行，未计完成轮次：** [原弱场补足第三方向反流]( {p}round731_drafts/gauge_counterflow_entry.md)原规范不变曲率泛函产生精确保Gauss的电动量平移，补齐572径向菜单缺少的z总动量；有正电能代价。下一步联立同态来源与新引力约束，正式保持730／3395。'".replace(']( {p}',']({p}'))
with (HERE/'round731_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)

