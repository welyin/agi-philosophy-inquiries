"""Publish723 entry while keeping formal count at722."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round722_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'722':'723','721':'722','3368':'3371','3102':'3116',
         'squared_reference_entry':'joint_reference_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**723入口已执行，未计完成轮次：** [四参考边缘与实际共同记录]({p}round723_drafts/joint_reference_entry.md)原径向诊断中，两态的T／s／Q_T／Q_s全部单参考谱分布及平均H相同，实际顺序记录概率却相差0.01287。完整Gauss提升和新仪器的共同来源映射继续审计，不把普通非交换机制另计新轮；正式保持722／3371。'")
with (HERE/'round723_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared723 entry checks.')
