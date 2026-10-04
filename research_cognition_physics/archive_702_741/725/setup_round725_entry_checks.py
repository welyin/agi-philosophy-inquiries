"""Publish725 entry without changing formal724 count."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round724_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'724':'725','723':'724','3374':'3377','3130':'3144',
         'relational_wall_entry':'wall_domain_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**725入口已执行，未计完成轮次：** [原壁位置与法向导数域]({p}round725_drafts/wall_domain_entry.md)556／563已有径向域警告；原651壁的目标流可有限时到轴，不能直接复用652光滑完备流证明。下一步核闭形式与原能源、真实区域的共同接法，不断言所有自伴实现失败；正式保持724／3377。'")
with (HERE/'round725_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared725 entry publication.')
