"""Publish726 audit while formal725 remains the latest scientific round."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round725_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'725':'726','724':'725','3377':'3380','3144':'3158',
         'wall_domain_entry':'common_wall_limit_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**726入口已执行，未计完成轮次：** [原来源与共同尺度]( {p}round726_drafts/common_wall_limit_entry.md)574已给同源严格Gauss半经典态，直接复用；原节点体积进入位置读取的法向反作用，固定图极限不能直接扩为空间细化。下一步核同一实际后态的位置与法向集中；正式保持725／3380。'".replace(']( {p}',']({p}'))
with (HERE/'round726_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared726 entry checks.')
