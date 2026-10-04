"""Prepare722 entry publication without modifying frozen STATUS."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
p=(HERE/'round721_drafts/check_and_publish_entry.py').read_text('utf8')
mapping={'721':'722','720':'721','3365':'3368','3088':'3102',
         'material_velocity_entry':'squared_reference_entry'}
p=re.sub('|'.join(re.escape(k) for k in mapping),lambda m:mapping[m.group()],p)
old=next(x for x in p.splitlines() if x.startswith('summary='))
p=p.replace(old,"summary='**722入口已执行，未计完成轮次：** [平方参考与原能源尾部]({p}round722_drafts/squared_reference_entry.md)原T流可转为平移坐标，平方流的tent准备产生严格多项式尾；是否构成原完整Gauss能源反例，仍待半密度提升、径向形式界及真实Kraus审计。正式保持721／3368，不以有限矩阵截断替代域检验。'")
with (HERE/'round722_drafts/check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(p)
print('Prepared722 entry checks.')

