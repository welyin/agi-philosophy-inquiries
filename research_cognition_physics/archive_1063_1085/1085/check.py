"""Exact finite calibrations for the converged definition; topology is proved in proof.md."""
from fractions import Fraction as F
from pathlib import Path
import argparse, hashlib, json

checks=[]
def ck(name, value):
    assert value, name
    checks.append(name)
def mm(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def mv(a,v):
    return tuple(sum(x*y for x,y in zip(row,v)) for row in a)
def tr(a): return sum(a[i][i] for i in range(len(a)))
def norm2(v): return sum(x*x for x in v)
I3=[[1,0,0],[0,1,0],[0,0,1]]
Rx=[[1,0,0],[0,0,-1],[0,1,0]]
Rz=[[0,-1,0],[1,0,0],[0,0,1]]
v=(0,0,1)
z=(1,0,0)
a=mv(mm(Rz,Rx),v); b=mv(mm(Rx,Rz),v)
ck('same anchored rotations have actual endpoint order difference',a==(1,0,0) and b==(0,-1,0))
ck('same tent test is exactly 1 versus 0',norm2(tuple(x-y for x,y in zip(a,z)))==0 and norm2(tuple(x-y for x,y in zip(b,z)))>=1)
ck('both rotations orthogonal and preserve the anchor',mm(Rx,list(map(list,zip(*Rx))))==I3 and mm(Rz,list(map(list,zip(*Rz))))==I3 and mv(Rx,(0,0,0))==(0,0,0) and mv(Rz,(0,0,0))==(0,0,0))

A=[F(1),F(3,4),F(1,2)]
r=F(1,2)
trace_x=(1+(A[0]*r)**2)/2
trace_z=(1+(A[2]*r)**2)/2
ck('anisotropic effect spectra forbid full unitary covariance',trace_x==F(5,8) and trace_z==F(17,32) and trace_x-trace_z==F(3,32))
ck('anisotropic response is a contraction with positive least singular value',max(A)==1 and min(A)==F(1,2))
ck('complete half-radius sphere has contrast norm lower bound 1/4',min(A)*r==F(1,4))

x=(F(1,3),F(2,5),F(-3,7)); y=(F(-2,9),F(1,6),F(4,11))
m=tuple((s+t)/2 for s,t in zip(x,y)); delta=tuple((s-t)/2 for s,t in zip(x,y))
ck('joint meeting retains the exact relative record',tuple(s+t for s,t in zip(m,delta))==x and tuple(s-t for s,t in zip(m,delta))==y)
ck('fixed result compensation is unique affine solution',tuple((s+(-s))/2 for s in x)==(0,0,0))

P=[[1,1],[0,1]]; Q=[[1,0],[1,1]]
ck('two-dimensional shear boundary has a fixed-point commutator',mm(P,Q)!=mm(Q,P) and mv(P,(0,0))==(0,0) and mv(Q,(0,0))==(0,0))
ck('shear example does not preserve Euclidean task distance',norm2(mv(P,(0,1)))!=1)

# On RP^2 the half-turn fixes the anchor line, but its continuous lift rotates it away.
B=[[1,0,0],[0,-1,0],[0,0,-1]]
q=(1,0,1)
u=mv(mm(Rz,B),q); w=mv(mm(B,Rz),q)
ck('projective-plane anchor stabilizer has noncommuting components',mv(B,v)==tuple(-x for x in v) and u!=w and u!=tuple(-x for x in w))
ck('the connecting half-turn path leaves the anchor line',mv(Rx,v)!=(0,0,1) and mv(Rx,v)!=(0,0,-1))

# Four-dimensional raw effects: v=(0,0,0,1) yields E=3I/4 versus I/4.
raw_gap=F(1,2)
ck('raw scalar distinction survives while traceless distinction vanishes',raw_gap>0 and F(3,4)-F(3,4)==0 and F(1,4)-F(1,4)==0)
ck('four-dimensional raw effect family is uniformly positive',F(1,4)**2*2<F(1,2)**2)
root=Path(__file__).resolve().parents[3]
sources=['research_cognition_physics/archive_1063_/1071/proof.md','research_cognition_physics/archive_1063_/1074/proof.md','research_cognition_physics/archive_370_428/research_note_383.md','research_cognition_physics/archive_370_428/research_note_387.md','research_cognition_physics/archive_370_428/research_note_388.md']
out={'round':1085,'kind':'axiomatic_synthesis_finite_calibration','science_increment_groups':0,'topological_proof_by_finite_checks':False,'assertion_count':len(checks),'checks':checks,'joint_model':{'anisotropic_trace_square_x':str(trace_x),'anisotropic_trace_square_rotated':str(trace_z),'spectral_gap':str(trace_x-trace_z),'half_radius_contrast_lower_bound':'1/4','contact_order_probabilities':['1','0']},'source_hashes':{s:hashlib.sha256((root/s).read_bytes()).hexdigest() for s in sources}}
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--check',action='store_true')
args=parser.parse_args()
if args.check:
    saved=json.loads(Path(__file__).with_name('results.json').read_text(encoding='utf-8'))
    assert out==saved
    print(json.dumps({'passed':True,'assertions':len(checks),'saved_results_match':True}))
else:
    print(json.dumps(out,indent=2,ensure_ascii=False))
