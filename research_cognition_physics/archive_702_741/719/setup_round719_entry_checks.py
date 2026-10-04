"""Prepare719 entry publication without altering its frozen status."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round718_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'718':'719','717':'718','3356':'3359','3046':'3060',
         'reflected_history_entry':'causal_accuracy_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**719入口已执行，未计完成轮次：** [旧因果限制与有限精度]({p}round719_drafts/causal_accuracy_entry.md)复用633见证：在尚不可跨区传信的时间窗，原分布读取至少留下0.2165的统一概率误差；共同局部表示输送不改变此条件。正式保持718／3359。下一步核真实局部记录任务，停止重复瞬时反例，工程设计后置。'")
with (HERE/'round719_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared719 entry checks.')
