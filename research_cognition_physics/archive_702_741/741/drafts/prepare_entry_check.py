"""Prepare the741 entry publisher; preserve all prior entry evidence."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
s=(HERE.parent/'round740_drafts/check_and_publish_entry.py').read_text('utf8')
m={'740':'741','739':'740','3417':'3419','3345':'3359','total_pole_entry':'absolute_source_order_entry'}
s=re.sub('|'.join(re.escape(k) for k in sorted(m,key=len,reverse=True)),lambda x:m[x.group()],s)
s=s.replace("checks['display_formulas']==5","checks['display_formulas']==7")
start=s.index('summary=heading+');end=s.index('\nsummary=summary.replace',start)
s=s[:start]+"summary=heading+' [实际绝对来源与共同圈阶]( {p}round741_drafts/research_note_741_working.md)原同一准备与记录的完整历史响应按统一阶次在二阶进入反作用；首阶只需领先背景上的绝对源。原对象校准通过，下一项接完整初始约束、背景族与真实参考的共同残差。'"+s[end:]
s=s.replace('pole evidence','source-order evidence')
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(s)
