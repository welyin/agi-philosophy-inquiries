"""Freeze651: common material coordinates and moving boundary conditions."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        pat=re.escape(key)
        if key.isdecimal():pat=r'(?<!\d)'+pat+r'(?!\d)'
        pats.append(pat)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_gravity_boundary_flux':'joint_material_boundary_transport',
    '649':'650','650':'651','651':'652','3185':'3188','3188':'3191',
    '1259':'1262','1262':'1265','2155':'2165','2165':'2175',
    'boundary_flux_entry.md':'material_wall_joint_entry.md',
    'original_boundary_flux_probe':'material_wall_signature_probe'}
verify=replace((HERE/'verify_round650.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original gravity constraint completion and common boundary flux.',
    'Verify original causal material wall and full boundary-geometry transport.')
write('verify_round651.py',verify)
publish=replace((HERE/'publish_round650.py').read_text('utf8'),mapping|{'## 296.':'## 297.','## 201.':'## 202.'})
start=publish.index('summary=');end=publish.index('planned={}',start)
header="""summary=('**第651轮完成：** [原物质参考、因果区域与完整边界变分的共同连接]({p}research_note_651.md)'
         '同一原解的固定s面类空，原时钟组合给局部类时侧壁；'
         '原参考约束关系度规与F，同一嵌入输送完整边界及透射匹配。'
         '三组、十六式通过，最新651／3191，1265份编号科学文件、2175份保护证据。'
         '[核验]({p}research_round_651_checks.json)、[条件账]({p}unified_physics_condition_ledger_651.md)。'
         '限给定作用的局部经典接口；诊断族离壳，量子参考及全尺度统一仍开放。')
order=('**当前执行顺序（651后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[652原关系参考、完整量子过程与共同区域表示]({p}round652_drafts/STATUS.md)，'
       '回查原量子域、来源和参考的真实映射，不继续侧壁或坐标精度优化；目标不改。')
"""
publish=publish[:start]+header+publish[end:]
publish=publish.replace('面积匹配与实际物质—引力通量','原物质图册与共同边界输送')
publish=publish.replace('区域演化还须保持完整边界配对','物质标签须与因果边界共同相容')
publish=publish.replace('原引力约束、面积匹配与物质边界通量的联合条件','原物质参考、因果区域与完整边界变分的共同连接')
write('publish_round651.py',publish)
write('postcheck_round651.py',replace((HERE/'postcheck_round650.py').read_text('utf8'),mapping))
write('round651_drafts/research_note_651_draft.md',(HERE/'research_note_651.md').read_text('utf8'))
review="""651 primary-agent review; no independent agent review.
Previous goal turn is progress:650 completed and published with2165 protected
evidence files and8332 valid navigation links;651 original-wall probe saved.
Read root/research navigation,650 report/results,650 ledger and651 entry.
No Python process found. No task, automation or goal changes.
Read573,618,619,647,648. Search old reports for moving boundaries/relational
metrics.648 already has the same Gram data: not claimed as a new discovery.
Primary BFR2.5-2.6 and Harlow-Wu3.5 read. Tools are mature and input conditions
are mapped explicitly, not used to assert quantum or physical equivalence.
On the actual573 constrained solution, dh spatial=0 at the stated point.
Inherited continuous bounds imply B<0, whereas s-lambda_star h has N^2>0.
lambda_star is chosen ONCE from the original point and then held fixed under
field variation. Its value does not depend on numerical psi. No global wall.
Y=(h,s-lambda h,A,B) inherits rank. gbar^00=Y2 and the stated gbar^11 combination
equalsY3 by the original derivative definitions. Fbar is a known function of
Y0,Y1. These are chart identities, not field equations or degree-of-freedom
counts. They do not eliminate physical matter variations.
Same e_Phi pulls all fields, bulk action, GHY and corners to the fixed domain.
Normal and K variations use the total kappa=delta g+Lie_u g. The scalar normal
derivative includes delta normal and delta gauge connection, not only sigma.
Original nonminimal momentum formulas and their frame map are retained.
The action identity is exact locally; the difference between moving-domain
Theta+i_u L and dressed Theta has an on-shell dQ_u term. Delta u and corners
must remain. No claim that all regional symplectic potentials are identical.
Yang-Mills is retained in the analytic boundary form and full transmission.
Old fermion background is zero. No full chiral quantum gluing claimed.
Transmission is pulled back on restrictions of the SAME smooth solution with
the SAME embedding and gauge identification. It does not solve arbitrary
regional initial/boundary data or make each region autonomous.
Numerics separate original on-shell wall from a NEW OFF-SHELL chart family.
That family satisfies X_Phi(Y)=Y identically, with original F and target, and
tests geometric derivatives and old momentum conversion. It does not stand
in for573's unknown time derivatives or actual Einstein evolution.
Independent centered differences test n,K,v against tensor variations; steps
2e-6,1e-6,5e-7 give quadratic decrease, final normalized error<6.53e-6.
Initial implementation wrongly requested an absent coefficient key d; replaced
it by H*S/12. Larger difference steps failed truncation thresholds; steps
reduced, threshold unchanged. Only corrected successful outputs are frozen.
16 numbered formulas; a stray terminal slash removed before freeze. No image
inspection. Next returns to quantum/continuum common objects, not wall scans.
"""
for name in ('research_note_651.md','joint_material_boundary_transport.py',
             'joint_material_boundary_transport_results.json','unified_physics_condition_ledger_651.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round651_drafts/final_review.txt',review)
print('651 frozen; verification and publication prepared')
