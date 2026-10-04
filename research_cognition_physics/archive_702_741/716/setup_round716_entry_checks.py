"""Create the non-numbered716 entry checker from the frozen715 workflow."""
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent
text=(HERE/'round715_drafts/check_and_publish_entry.py').read_text('utf8')
text=text.replace('715','__ENTRY__').replace('714','715').replace('__ENTRY__','716')
text=text.replace('observable_gain_entry','sterile_mode_entry')
text=text.replace('3002','3016').replace("display_formulas']==4","display_formulas']==3")
text=text.replace("'sterile_mode_entry.md','check_and_publish_entry.py')",
                  "'sterile_mode_entry.md','check_and_publish_entry.py','first_action_diagnostic.json')")
summary="**716入口已执行，未计完成轮次：** [原sterile模与完整CAR读取]({p}round716_drafts/sterile_mode_entry.md)633占据数仪器在原有限图保持Gauss与有限能源域；实际差量由原质量／跳跃的至多三模交叉块决定，64模式校准通过。正式保持715／3349。下一步核平滑体积归一、原参考、传播和来源的共同尺度条件，保留633因果实现限制。"
text=re.sub(r"summary='[^\n]*'",lambda _:'summary='+repr(summary),text,count=1)
target=HERE/'round716_drafts/check_and_publish_entry.py'
with target.open('x',encoding='utf8') as f:f.write(text)
print('Created716 entry checker.')
