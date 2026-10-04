"""Prepare720 entry publication; frozen STATUS remains untouched."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round719_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'719':'720','718':'719','3359':'3362','3060':'3074',
         '2026-10-03':'2026-10-04','20261003':'20261004',
         'causal_accuracy_entry':'spin_propagation_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**720入口已执行，未计完成轮次：** [联合自旋参考与原空间传播]({p}round720_drafts/spin_propagation_entry.md)原自旋标量跳跃分支有全spin保护，但固定原spin／动量字典的内部对称矩阵与633 Weyl主部相距至少|k|。实际协变须同时旋转空间方向。正式保持719／3362；下一步核物质—方向联合参考，停止保护编码和倍增扫描。'")
with (HERE/'round720_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared720 entry checks.')
