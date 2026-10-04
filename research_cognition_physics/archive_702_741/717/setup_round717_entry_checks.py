"""Prepare717 non-numbered entry publication without changing frozen status."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round716_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'716':'717','715':'716','3349':'3352','3016':'3031',
         'sterile_mode_entry':'postrecord_reference_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
p=p.replace(",'first_action_diagnostic.json'","")
line=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(line,"summary='**717入口已执行，未计完成轮次：** [实际读取后的原参考与反馈]({p}round717_drafts/postrecord_reference_entry.md)非选择CAR读取保玻色边缘与原时间交换子平方，却一般不保热平稳性；96模式校准通过。正式保持716／3352。下一步核原质量与玻色动量如何传递记录后的关联，不继续夹断或波包优化。'")
with (HERE/'round717_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:
    f.write(p)
print('Prepared717 entry verification and navigation snapshot publication.')
