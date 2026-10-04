"""Create only the non-numbered715 publication checker."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
text=(HERE/'round713_drafts/check_and_publish_entry.py').read_text('utf8')
text=text.replace('713','__ENTRY__').replace('712','714').replace('__ENTRY__','715')
text=text.replace('local_interaction_entry','observable_gain_entry')
text=text.replace('2971','3002').replace("display_formulas']==2","display_formulas']==4")
summary="**715入口已执行，未计完成轮次：** [复合观测重标定与实际读取]({p}round715_drafts/observable_gain_entry.md)原F_M仿射读口要求增益受效果正性限制；有界非线性可合法，但其精确最坏注能随增益平方增长。此为713直接推论，正式保持714／3346。下一步共同检验实际观测映射与原参考分布，不把最坏预算当真实热平均，不靠形式放大签收连续观测。"
text=re.sub(r"summary='[^\n]*'",lambda _:'summary='+repr(summary),text,count=1)
target=HERE/'round715_drafts/check_and_publish_entry.py'
with target.open('x',encoding='utf8') as f:f.write(text)
print('Created715 entry checker.')
