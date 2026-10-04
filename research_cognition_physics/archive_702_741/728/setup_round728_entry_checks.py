"""Preserve728 entry while formal727 count stays fixed."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round727_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'727':'728','726':'727','3383':'3386','3172':'3186',
         'fermion_reference_entry':'macroscopic_source_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**728入口已执行，未计完成轮次：** [宏观来源与真实记录后的参考]({p}round728_drafts/macroscopic_source_entry.md)旧亚隙谱可直接压低原基态的时间平均噪声；同一局部占据数记录却保留零频来源及能源涨落。591／602直接复用，下一步核记录后的共同动态参考及反作用，正式保持727／3386。'")
with (HERE/'round728_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared728 entry checks.')
