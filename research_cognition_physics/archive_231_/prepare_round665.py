"""Prepare665 common Gauss/contact process and preserve all earlier evidence."""
import hashlib
import re
from pathlib import Path
HERE=Path(__file__).resolve().parent


def replace(text,mapping):
    pats=[]
    for key in sorted(mapping,key=len,reverse=True):
        p=re.escape(key)
        if key.isdecimal():p=r'(?<!\d)'+p+r'(?!\d)'
        pats.append(p)
    return re.sub('|'.join(pats),lambda m:mapping[m.group()],text)


def write(name,text):
    p=HERE/name;p.parent.mkdir(exist_ok=True)
    with p.open('x',encoding='utf8',newline='\n') as f:f.write(text)


mapping={'joint_connection_matter_matching':'joint_contact_gauss_history',
    '663':'664','664':'665','665':'666','3231':'3234','3234':'3237',
    '1301':'1304','1304':'1307','2299':'2306','2306':'2316'}
verify=replace((HERE/'verify_round664.py').read_text('utf8'),mapping)
verify=verify.replace('Verify original scalar/CAR connection branches, matching and common sources.',
    'Verify whole-graph contact scope, original auxiliary lift and actual spatial sources.')
a=verify.index('    import joint_spinorial_geometry_sources');b=verify.index('    result=core.read(model.TARGET)',a)
verify=verify[:a]+"""    import joint_connection_matter_matching as prior_model
    assert prior_model.run()==core.read(prior_model.TARGET)
    import importlib.util
    probe_path=HERE/'round665_drafts/spin_charge_auxiliary_probe.py'
    spec=importlib.util.spec_from_file_location('round665_probe',probe_path)
    probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
    assert probe.run()==core.read(probe_path.with_name('spin_charge_auxiliary_probe_results.json'))
"""+verify[b:]
verify=verify.replace("text['display_formulas']==16","text['display_formulas']==18")
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_665.md','round665_drafts/research_note_665_draft.md',
           'round665_drafts/final_review.txt','round666_drafts/STATUS.md',
           'round665_drafts/auxiliary_current_entry.md','round665_drafts/spin_charge_auxiliary_probe.py',
           'round665_drafts/spin_charge_auxiliary_probe_results.json')
"""+verify[b:]
write('verify_round665.py',verify)
publish=replace((HERE/'publish_round664.py').read_text('utf8'),mapping|{'## 310.':'## 311.','## 215.':'## 216.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第665轮完成：** [原Gauss量子过程、接触作用与辅助历史的共同条件]({p}research_note_665.md)'
         '原固定图的声明接触扩展保共同域、Gauss热态和原记录；辅助提升、归一及空间能源流须同步。'
         '三组、十八式通过，最新665／3237，1307份编号科学文件、2316份保护证据。'
         '[核验]({p}research_round_665_checks.json)、[全条件账]({p}unified_physics_condition_ledger_665.md)。'
         '限保原K的正几何有限图；连续测度及量子GR仍开放。')
order=('**当前执行顺序（665后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[666原CAR到全左手表示的排序与共同测度]({p}round666_drafts/STATUS.md)，'
       '核真实粒子—空穴字典、真空与原辅助对象；目标不改。')
"""+publish[b:]
publish=publish.replace('原物质、联络选择与共同几何来源的匹配条件','原Gauss量子过程、接触作用与辅助历史的共同条件')
publish=publish.replace('原物质与联络分支的共同匹配','原完整Gauss过程与接触辅助表示')
publish=publish.replace('联络选择同时约束标量几何、物质与来源','同一辅助表示须保真空提升、归一和空间来源')
write('publish_round665.py',publish)
write('postcheck_round665.py',replace((HERE/'postcheck_round664.py').read_text('utf8'),mapping))
write('round665_drafts/research_note_665_draft.md',(HERE/'research_note_665.md').read_text('utf8'))
review="""665 primary-agent proof/code audit, no independent agent review.
Previous goal turn completed/published664 and saved actual665 entry: progress.
Read navigation,664 report/results/postcheck,665 entry and saved probe; review
603,623/624,643,655,662 and original32 CAR, mass, hopping and record definitions.
No live Python process; no assumed waiting handle or unrequested restart.
No goal edits, app tasks, automation, images or subagents.
Full-graph statements apply to CE or a K-matched contact extension, NOT minimal
Jordan Cartan target. Original no-contact branch remains separately preserved.
V norm <=204 sum1/w follows full32 current norms and vacuum normal order.
Fixed positive geometry finite graph; bounded perturbation preserves core/domain,
Gauss/fermion parity and trace-class heat. Bounds from min-max, not exponential
operator monotonicity. No graph-uniform or w->0 claim.
Original scalar records commute V, preserve instantaneous energy injection;
probabilities, equilibrium states and subsequent dynamics generally differ.
All source statements include original H derivatives, not only a conditional trace.
Normalized auxiliary charge +five symmetric spin factors are explicitly given.
Finite spin split not exactly rotation invariant. Internal gauge exact per field.
Positive mean transfer does not imply positive integrands or classical ontology.
Noncommuting mass/hopping preserved in time order. No frozen spin/charge shortcut.
E A^dag A <=exp(1152a); iterative HS bound supplies uniform auxiliary mean
control. Whole bosonic bridge uses original confinement and Gauss-twisted endpoints.
Auxiliary mean limit defined before continuum path claims; raw differentiated
scores need not have uniformly bounded variation. Duhamel controls full sources.
Nambu GL lift uses inverse transpose, not conjugate. e^(2aN) gives det r=e^(64a)
and trace prefactor e^(32a). Actual32 quark determinant/lepton Fock validates it.
Lift sign/phase follows original Fock word, never absolute value or principal root.
Six physical auxiliary densities each have +3 normalization derivative. Missing
these yields -18 only for this specified step, not cosmological vacuum energy.
Two-node neutral modes are actual invariant conditional original sector. Spatial
hopping .23 independent of geometry is an allowed diagnostic input; not unique.
Contact current and its source checked with same mass/geometry/thermal state.
Three new groups; saved entry tests rerun but not recounted.18 formulas.
Next666 checks original particle-hole/left-handed ordering and actual auxiliary
measure before claiming same four-branch object. No cognitive system design.
"""
for name in ('research_note_665.md','joint_contact_gauss_history.py',
             'joint_contact_gauss_history_results.json','unified_physics_condition_ledger_665.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round665_drafts/final_review.txt',review)
print('665 report/review frozen; verification and publication prepared')
