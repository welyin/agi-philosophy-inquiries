"""Preserve727 reference audit while formal726 count stays fixed."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round726_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'726':'727','725':'726','3380':'3383','3158':'3172',
         'common_wall_limit_entry':'fermion_reference_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**727入口已执行，未计完成轮次：** [替代费米参考与传播]( {p}round727_drafts/fermion_reference_entry.md)全内部单spin填充可保Gauss并消去原质量／spin标量跳跃，但不适用于全部原Weyl方向。642／709及720限制直接复用；下一步核随物质与传播共同变化的参考及来源，正式保持726／3383。'".replace(']( {p}',']({p}'))
with (HERE/'round727_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared727 entry checks.')
