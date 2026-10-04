"""Prepare721 entry publication without modifying frozen STATUS."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round720_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'720':'721','719':'720','3362':'3365','3074':'3088',
         'spin_propagation_entry':'material_velocity_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**721入口已执行，未计完成轮次：** [原物质参考与实际读取]({p}round721_drafts/material_velocity_entry.md)652已有参考自伴性与共同表示，直接复用；合法速度读口在原64维诊断通过，但完整H的能源形式域和来源仍须共同核验。正式保持720／3365；下一项检查原目标流、全部质量与记录资源，停止保护spin分类。'")
with (HERE/'round721_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared721 entry checks.')

