"""Publish724 entry while keeping formal count at723."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round723_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'723':'724','722':'723','3371':'3374','3116':'3130',
         'joint_reference_entry':'relational_wall_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**724入口已执行，未计完成轮次：** [同一选面与完整边界资料]({p}round724_drafts/relational_wall_entry.md)原h壁换成线性T壁可在一点保法向，却改变外曲率；因果平方还须保同一混合项。优先保原壁精确函数并核实际区域筛选，651／652直接推论不新增轮次；正式保持723／3374。'")
with (HERE/'round724_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared724 entry checks.')
