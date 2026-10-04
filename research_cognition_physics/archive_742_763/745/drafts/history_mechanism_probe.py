"""Diagnostic comparisons only; changed potentials are not the original model."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'round744_drafts'))
import history_density_probe as density
old=density.old
native_potential=old.old.potential
rows=[]
for strength in (0.,1.,2.):
    old.old.potential=lambda p:{k:strength*v for k,v in native_potential(p).items()}
    a=old.more_jets(0);b=old.more_jets((1<<30)|(1<<31))
    pa,pb=density.density(a,6),density.density(b,6)
    delta={k:pb.get(k,0)-pa.get(k,0) for k in set(pa)|set(pb)}
    delta={k:z for k,z in delta.items() if abs(z)>1e-8}
    av,bv=old.measure(a,64),old.measure(b,64)
    rows.append(dict(potential_strength=strength,difference_sixth=bv['sixth']-av['sixth'],
        coefficients=[dict(r=k[0],s=k[1],real=z.real,imag=z.imag) for k,z in sorted(delta.items())]))
old.old.potential=native_potential
result=dict(purpose='Identify coefficient dependence before analytic certification. Only strength1 is native.',
    rows=rows,original_H_used_for_physical_claim=True,other_rows_are_explicit_counterfactual_diagnostics=True)
with Path(__file__).with_name('history_mechanism_probe_results.json').open('x',encoding='utf8') as f:
    json.dump(result,f,ensure_ascii=False,indent=2)
print(json.dumps([{k:r[k] for k in ('potential_strength','difference_sixth')} for r in rows]))
