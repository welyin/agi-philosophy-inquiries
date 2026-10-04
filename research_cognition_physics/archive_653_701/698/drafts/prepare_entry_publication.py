"""Prepare the checked698 entry publication from the existing bounded protocol."""
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent
text=(HERE.parent/'round697_drafts/check_and_publish_entry.py').read_text('utf8')
text=re.sub(r'(?<!\d)(696|697)(?!\d)',lambda m:str(int(m.group())+1),text)
text=text.replace('center_holonomy_entry','polynomial_error_entry')
text=text.replace('original centre-holonomy sewing','original-source polynomial error contract')
a=text.index('    entry=');b=text.index('    head,rest=',a)
text=text[:a]+'''    entry=(marker+' [原全积分的有限多项式与严格截断误差]('+prefix+'round698_drafts/polynomial_error_entry.md)'
        '非对角中心项有精确有限Laurent界；静态项157阶逼近的全平均截断误差已给有理上界。'
        '实际多项式平均尚未求值，严格负号仍待认证；正式科学轮次仍为697，旧空间接口与目标不变。')
'''+text[b:]
text=text.replace("'primary_source_audit.json'","'polynomial_error_entry_first.py','prepare_entry_publication.py'")
text=text.replace("assert not TARGET.exists()","assert not TARGET.exists()\nassert core.text_checks(HERE/'polynomial_error_entry.md')['display_formulas']==8")
with (HERE/'check_and_publish_entry.py').open('x',encoding='utf8',newline='\n') as f:f.write(text)
print('698 entry publication prepared')
