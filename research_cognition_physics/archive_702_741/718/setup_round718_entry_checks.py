"""Prepare718 entry publication, leaving its frozen status untouched."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round717_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'717':'718','716':'717','3352':'3356','3031':'3046',
         'postrecord_reference_entry':'reflected_history_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**718入口已执行，未计完成轮次：** [原比较历史与真实首次记录]({p}round718_drafts/reflected_history_entry.md)同一H与RHR可表示原非选择记录后的完整有限历史；保留首次结果还必须保交叉历史，原CAR校准通过。正式保持717／3356。下一步接原共同域、有限时间与几何来源，停止反射恒等式及高矩扫描。'")
with (HERE/'round718_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:
    f.write(p)
print('Prepared718 entry checks.')
