"""Preserve730 entry; formal729 count stays fixed."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round729_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'729':'730','728':'729','3389':'3392','3200':'3214',
         'nonflat_reference_entry':'varying_reference_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**730入口已执行，未计完成轮次：** [原变化质量与连续参考]( {p}round730_drafts/varying_reference_entry.md)原h／s调制使质量本征框架随位置变化；630／633常量参考不能原样搬用。原复Y、BdG质量和框架联络已校准，下一项保变化连接核共同连续参考。正式保持729／3392。'".replace(']( {p}',']({p}'))
with (HERE/'round730_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)

