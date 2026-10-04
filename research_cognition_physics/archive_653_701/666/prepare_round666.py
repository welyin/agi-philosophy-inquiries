"""Prepare666 reference/mass/measure integration without rewriting evidence."""
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


mapping={'joint_contact_gauss_history':'joint_reference_mass_stabilizer',
    '664':'665','665':'666','666':'667','3234':'3237','3237':'3239',
    '1304':'1307','1307':'1310','2306':'2316','2316':'2329'}
verify=replace((HERE/'verify_round665.py').read_text('utf8'),mapping)
verify=verify.replace('Verify whole-graph contact scope, original auxiliary lift and actual spatial sources.',
    'Verify reference/mass common stabilizer and transported actual auxiliary weights.')
a=verify.index('    import joint_connection_matter_matching');b=verify.index('    result=core.read(model.TARGET)',a)
verify=verify[:a]+"""    import joint_contact_gauss_history as prior_model
    assert prior_model.run()==core.read(prior_model.TARGET)
    import importlib.util
    probe_path=HERE/'round666_drafts/particle_hole_contact_probe.py'
    spec=importlib.util.spec_from_file_location('round666_probe',probe_path)
    probe=importlib.util.module_from_spec(spec);spec.loader.exec_module(probe)
    assert probe.run()==core.read(probe_path.with_name('particle_hole_contact_probe_results.json'))
"""+verify[b:]
verify=verify.replace("text['display_formulas']==18","text['display_formulas']==16")
verify=verify.replace("==(3,0,0)","==(2,0,0)").replace('fresh_tests=dict(run=3,','fresh_tests=dict(run=2,')
a=verify.index('    names=');b=verify.index('    preserved=',a)
verify=verify[:a]+"""    names=('unified_physics_condition_ledger_666.md','round666_drafts/research_note_666_draft.md',
           'round666_drafts/final_review.txt','round667_drafts/STATUS.md',
           'round666_drafts/particle_hole_ordering_entry.md','round666_drafts/particle_hole_contact_probe.py',
           'round666_drafts/particle_hole_contact_probe_results.json',
           'round666_drafts/transport_diagnostics.py','round666_drafts/transport_first_failure.json',
           'round666_drafts/joint_reference_mass_stabilizer_first_attempt.py')
"""+verify[b:]
write('verify_round666.py',verify)
publish=replace((HERE/'publish_round665.py').read_text('utf8'),mapping|{'## 311.':'## 312.','## 216.':'## 217.'})
a=publish.index('summary=');b=publish.index('planned={}',a)
publish=publish[:a]+"""summary=('**第666轮完成：** [原参考排序、质量与内部对称性的共同条件]({p}research_note_666.md)'
         '给定载体中的原参考分级与Majorana方向共同保留原SM连通商群；有限真空、辅助权重与来源须同步。'
         '两组、十六式通过，最新666／3239，1310份编号科学文件、2329份保护证据。'
         '[核验]({p}research_round_666_checks.json)、[全条件账]({p}unified_physics_condition_ledger_666.md)。'
         '载体及内部划分仍为输入；原手征测度、连续与量子GR仍开放。')
order=('**当前执行顺序（666后，优先于下方历史安排）：** 先整合共同条件，认知设计后置。'
       '接[667有限参考字典与跨尺度状态表示]({p}round667_drafts/STATUS.md)，'
       '核局部代数、参考态与尺度嵌入的共同连接；目标不改。')
"""+publish[b:]
publish=publish.replace('原Gauss量子过程、接触作用与辅助历史的共同条件','原参考排序、质量与内部对称性的共同条件')
publish=publish.replace('原完整Gauss过程与接触辅助表示','原参考、质量与同一内部对称性')
publish=publish.replace('同一辅助表示须保真空提升、归一和空间来源','参考字典须与状态、权重和来源共同输送')
write('publish_round666.py',publish)
write('postcheck_round666.py',replace((HERE/'postcheck_round665.py').read_text('utf8'),mapping))
write('round666_drafts/research_note_666_draft.md',(HERE/'research_note_666.md').read_text('utf8'))
review="""666 primary-agent proof/code audit; no independent agent review.
Previous goal turn completed665 and actual666 entry: progress. Read navigation,
latest note/results/entry and614,653,660-661,664-665. Search history before667.
No goal edits, app tasks, automation, images or subagents. Process query via
CIM denied access; Get-Process fallback returned no Python. No review escalation.
Existing finite32 CAR map is actual614, not a newly chosen left-handed vacuum.
N_a=16+Lambda number, Q_chi=16-N_b, spin unchanged. Contact reset differs by
34N_b-2Lambda number-288I. Entry lepton spectrum and logZ rerun, not recounted.
In original given carrier Lambda=-colour parity, mixed24 Gram16I is exact.
Reference commutant connected group is Spin6xSpin4/Z2. Majorana su5 result
inherited from614, not claimed new. Joint rank33, explicit SM12 basis verified.
Group intersection is mature Baez-Huerta theorem9, mapped to actual tensors.
Claim conditional on carrier,3+2,neutrino line,nonzero Majorana and normal order.
Pair preserves extra central -1; only connected stabilizer asserted as SM.
Do not equate Spin10 carrier or Pati-Salam stabilizer with new physical gauge.
Original gauge/Higgs/mass covariance still required; no fixed arbitrary full
mass-background symmetry claim. No continuum anomaly inferred from sorting.
Full32 trace uses four8 Fock blocks with colour diagonal closure, generic weak
and hypercharge and all original complex masses. No giant uncomputed Fock claim.
Right block determinants retained before multiplication. Total phase one.
Nonunitary Nambu transform uses inverse transpose. Reference factor exp32a+16s
retained; theta-source is -192a-48s. Conditional theta varies contact w only;
all other complete sources transported abstractly, not replaced by this number.
Normalized Gaussian variable translation gives exp288a and -32a N drift.
Initial32-point shifted quadrature failed. First attempt, diagnostic and failure
preserved.64/80 convergence passes unchanged3e-11 threshold; exact proof is
Gaussian completion of square, not numerical convergence alone.
Finite graph unitary does not automatically define infinite common Fock vacuum.
653 charge support obstruction survives. Gaussian and S9 auxiliaries distinct.
Four primary branches remain separate; no unique SM, continuous theory or GR.
Two fresh groups,16 equations. Next667 tests reference/scale compatibility,
not apparatus or cognitive design; full physics goal unchanged.
"""
for name in ('research_note_666.md','joint_reference_mass_stabilizer.py',
             'joint_reference_mass_stabilizer_results.json','unified_physics_condition_ledger_666.md'):
    review+=name+' '+hashlib.sha256((HERE/name).read_bytes()).hexdigest()+'\n'
write('round666_drafts/final_review.txt',review)
print('666 report/review frozen; verification and publication prepared')
