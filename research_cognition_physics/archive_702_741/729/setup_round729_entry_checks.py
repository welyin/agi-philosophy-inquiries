"""Preserve729 nonflat entry while formal728 count stays fixed."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round728_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'728':'729','727':'728','3386':'3389','3186':'3200',
         'macroscopic_source_entry':'nonflat_reference_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**729入口已执行，未计完成轮次：** [原非平坦链路与参考谱]( {p}round729_drafts/nonflat_reference_entry.md)恢复574实际弱／圆链路后，星图仍有正隙但平坦扰动证书不足；能源与来源同步变化。下一步核原周期图及真实几何动能映射，正式保持728／3389。'".replace(']( {p}',']({p}'))
with (HERE/'round729_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared729 entry checks.')
