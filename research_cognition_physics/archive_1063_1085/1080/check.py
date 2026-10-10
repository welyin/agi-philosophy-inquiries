"""1080: rational source/reader certificate for the common-drift angle witness.
No source scan. Two fixed analytic witnesses; stdout JSON only.
"""
from fractions import Fraction as F
from pathlib import Path
import importlib.util,json,math,hashlib
import numpy as np
BASE=Path(__file__).resolve().parent

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def run():
    src=load("source1077",BASE.parent/"1077/check.py")
    iv=load("interval1078",BASE.parent/"1078/check.py")
    reader=src.load_reader()
    _,_,rows,h=reader.system()
    degree=40;den=5**degree*math.factorial(degree)
    swap=np.array([6*(d^24 if ((d>>4)&1)!=((d>>3)&1) else d)+g for d in range(64) for g in range(6)])
    def ham(v):return np.sum(v[0][rows],axis=1),np.sum(v[1][rows],axis=1)
    def prop(v):
        r,i=v;orr=np.zeros(r.shape,dtype=object);oi=orr.copy()
        for k in range(degree+1):
            coeff=5**(degree-k)*(math.factorial(degree)//math.factorial(k))
            a,b=((r,i),(i,-r),(-r,-i),(-i,r))[k%4]
            orr+=coeff*a;oi+=coeff*b
            if k<degree:r,i=ham((r,i))
        return orr,oi
    *_,en,raw=reader.exact_reader()
    nr,ni,ns=raw
    br=[nr[:6,6:]+nr[6:,:6],ni[6:,:6]-ni[:6,6:],nr[:6,:6]-nr[6:,6:]]
    bi=[ni[:6,6:]+ni[6:,:6],nr[:6,6:]-nr[6:,:6],ni[:6,:6]-ni[6:,6:]]
    tail=F(9,5)**41/math.factorial(41)/(1-F(9,5)/42)
    tau=F(1,10**10)
    start=prop(src.initial_columns())
    pulse=(start[0][swap],start[1][swap]) # ideal theta=pi/2; global -i cancels
    y=prop(prop(pulse))
    points=[];bboxes=[];vboxes=[];centers=[]
    for usteps in (2,3):
        if usteps==3:y=prop(y)
        hr,hi=ham(y);dy=(hi,-hr)
        yr,yi=[v.reshape(64,6,28) for v in y]
        dr,di=[v.reshape(64,6,28) for v in dy]
        gr=np.zeros((6,6),dtype=object);gi=gr.copy();vr=gr.copy();vi=gr.copy()
        for a,b,c,d in zip(yr,yi,dr,di):
            gr+=a@a.T+b@b.T;gi+=b@a.T-a@b.T
            vr+=c@a.T+d@b.T+a@c.T+b@d.T
            vi+=d@a.T-c@b.T+b@c.T-a@d.T
        assert np.array_equal(gr,gr.T) and np.array_equal(gi,-gi.T)
        assert np.array_equal(vr,vr.T) and np.array_equal(vi,-vi.T)
        gs=384*den**(2*(1+usteps))
        def con(a,b,c,d):return F(int(np.sum(a*c.T-b*d.T)),gs*2*ns)
        bc=[con(a,b,gr,gi) for a,b in zip(br,bi)]
        vc=[con(a,b,vr,vi) for a,b in zip(br,bi)]
        e=(1+tail)**(1+usteps)-1
        ev=en+(1+en)*e*(2+e)+18*tau
        ej=18*en+(1+en)*18*e*(2+e)+324*tau
        assert ev<F(1,500000000) and ej<F(33,10**9)
        bb=[iv.interval(x,ev) for x in bc];vv=[iv.interval(x,ej) for x in vc]
        assert bb[2][0]>F(3,20)
        bboxes.append(bb);vboxes.append(vv);centers.append((bc,vc))
        points.append({"point":["1/5",str(F(usteps,5)),"pi/2"],
          "Taylor_propagator_count":1+usteps,"W_Taylor_error":str(e),
          "b_entry_error_upper":"1/500000000","v_entry_error_upper":"33/1000000000",
          "b_boxes":[iv.outward(a,10**14) for a in bb],
          "v_boxes":[iv.outward(a,10**14) for a in vv],
          "ideal_b_diagnostic":[float(x) for x in bc],
          "ideal_v_diagnostic":[float(x) for x in vc]})
    def dot(x,y):return iv.total(iv.mul(a,b) for a,b in zip(x,y))
    b1,b2=bboxes;v1,v2=vboxes
    s1,s2=dot(b1,b1),dot(b2,b2);c=dot(b1,b2)
    a=iv.add(dot(v1,b2),dot(b1,v2));d1,d2=dot(b1,v1),dot(b2,v2)
    witness=iv.sub(iv.mul(iv.mul(a,s1),s2),iv.mul(c,iv.add(iv.mul(d1,s2),iv.mul(s1,d2))))
    enclosure=iv.outward(witness,10**16)
    assert F(-1,10**8)<F(enclosure[0])<=F(enclosure[1])<F(-9,10**9)
    (b1,v1),(b2,v2)=centers
    dotc=lambda a,b:sum(x*y for x,y in zip(a,b))
    s1,s2=dotc(b1,b1),dotc(b2,b2);c=dotc(b1,b2)
    n=(dotc(v1,b2)+dotc(b1,v2))*s1*s2-c*(dotc(b1,v1)*s2+s1*dotc(b2,v2))
    assert F(enclosure[0])<n<F(enclosure[1])
    return {"round":1080,"status":"strict_certificate_passed",
      "source":"unchanged 480/1077 coherent W_tau(p), fixed 1076 reader",
      "finite_H_on_pulse_tau":str(tau),"reader_time":"7/10","Taylor_degree":40,
      "Taylor_step":"1/5","reader_error_bound":str(en),"points":points,
      "witness_definition":"N=(v1.b2+b1.v2)s1s2-(b1.b2)[(b1.v1)s2+s1(b2.v2)]",
      "finite_pulse_witness_enclosure":enclosure,"strict_upper":"-9/1000000000",
      "ideal_witness_diagnostic":float(n),
      "ideal_normalized_Gram_derivative_diagnostic":float(n)/float(s1*s2)**1.5,
      "source_hashes":{str(p.relative_to(BASE.parent)):hashlib.sha256(p.read_bytes()).hexdigest()
         for p in (BASE.parent/"1076/check.py",BASE.parent/"1077/check.py",BASE.parent/"1078/check.py")},
      "scope":{"whole_open_source_cube_exact_common_H_rotation_excluded":True,
        "positive_source_dependent_gain_allowed":True,"raw_bias_invariance_required":False,
        "output_must_remain_in_original_cube":False,"analytic_witness_is_new_control_permission":False,
        "all_subfamilies_excluded":False,"history_dependent_control_excluded":False,
        "finite_tolerance_macro_model_excluded":False,"actual_spatial_dimension_derived":False,
        "empirical_results":0,"new_accepted_cognitive_axioms":0}}
if __name__=="__main__":print(json.dumps(run(),ensure_ascii=False,indent=2))
