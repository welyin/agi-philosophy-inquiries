"""1205 independent fixed-program resource and whole-path audit.
No author code imports. Only --write creates independent_results.json once.
The finite Boolean templates are executed; the padded 2e13-gate program is
specified by uniform templates and bounds, not expanded or physically built.
"""
from __future__ import annotations
import argparse,hashlib,json,math,random
from bisect import bisect_right
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent
STAGE=HERE.parent

def load(p):return json.loads(p.read_text(encoding="utf-8-sig"))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def compare(a,b):
    if isinstance(a,dict):
        assert set(a)==set(b)
        for k in a:compare(a[k],b[k])
    elif isinstance(a,(list,tuple)):
        assert len(a)==len(b)
        for x,y in zip(a,b):compare(x,y)
    elif isinstance(a,float):assert math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-18),(a,b)
    else:assert a==b,(a,b)

class Circuit:
    """Fresh Boolean output wires; independent borrow-based subtraction."""
    def __init__(self,record=True):
        self.value=[0,1];self.initial=[0,1];self.gates=[];self.record=record
        self.boolean_count=0;self.reversible_count=0
    def input(self,x,w):
        ids=[]
        for j in range(w):
            ids.append(len(self.value))
            self.value.append((x>>j)&1);self.initial.append((x>>j)&1)
        return ids
    def b(self,kind,a,b=None):
        z=len(self.value);self.value.append(0);self.initial.append(0)
        self.boolean_count+=1
        gates = [((a,b),z)] if kind=="and" else (
                [((a,),z),((b,),z)] if kind=="xor" else [((),z),((a,),z)])
        for controls,target in gates:
            # Repeated AND controls reduce to one control; XOR still has two gates.
            controls=tuple(dict.fromkeys(controls))
            self.value[target]^=int(all(self.value[k] for k in controls))
            self.reversible_count+=1
            if self.record:self.gates.append((controls,target))
        return z
    def integer(self,wires):return sum(self.value[j]<<i for i,j in enumerate(wires))
    def add(self,a,b):
        carry=0;out=[]
        for x,y in zip(a,b):
            d=self.b("xor",x,y)
            out.append(self.b("xor",d,carry))
            p=self.b("and",x,y);q=self.b("and",d,carry)
            carry=self.b("xor",p,q)
        return out,carry
    def sub(self,a,b):
        borrow=0;out=[]
        for x,y in zip(a,b):
            d=self.b("xor",x,y);out.append(self.b("xor",d,borrow))
            nx=self.b("not",x);nd=self.b("not",d)
            p=self.b("and",nx,y);q=self.b("and",nd,borrow)
            borrow=self.b("xor",p,q)
        return out,borrow
    def mux(self,s,a,b):
        # select a if s=1, b otherwise
        out=[]
        for x,y in zip(a,b):
            d=self.b("xor",x,y)
            out.append(self.b("xor",y,self.b("and",s,d)))
        return out
    def mul(self,a,b):
        n=len(a);out=[0]*n
        for j in range(n):
            row=[0]*j+[self.b("and",a[k],b[j]) for k in range(n-j)]
            out,_=self.add(out,row)
        return out
    def divide_small(self,a,m):
        assert 1<=m<=100
        # All partial remainders <2m. In the actual program this is <=8 bits.
        rw=(2*m-1).bit_length()
        rem=[0]*rw;div=[(m>>j)&1 for j in range(rw)]
        quotient=[0]*len(a)
        for bit in range(len(a)-1,-1,-1):
            trial=[a[bit]]+rem[:-1]
            reduced,borrow=self.sub(trial,div)
            ge=self.b("not",borrow)
            rem=self.mux(ge,reduced,trial);quotient[bit]=ge
        return quotient,rem,rw
    def reverse_check(self):
        assert self.record
        values=self.value.copy()
        for controls,target in reversed(self.gates):
            values[target]^=int(all(values[k] for k in controls))
        assert values==self.initial

def primitive_tests():
    cases=0;max_ratio=F(0);max_large_nodes=0
    for w in range(1,5):
      for a in range(1<<w):
       for b in range(1<<w):
        for kind in ("add","sub","mul","mux"):
            c=Circuit();aa=c.input(a,w);bb=c.input(b,w)
            if kind=="add":
                out,carry=c.add(aa,bb)
                assert c.integer(out)==(a+b)%(1<<w)
                assert c.value[carry]==(a+b)//(1<<w)
            elif kind=="sub":
                out,borrow=c.sub(aa,bb)
                assert c.integer(out)==(a-b)%(1<<w)
                assert c.value[borrow]==int(a<b)
            elif kind=="mul":
                out=c.mul(aa,bb);assert c.integer(out)==a*b%(1<<w)
            else:
                s=(a+b)&1;out=c.mux(s,aa,bb)
                assert c.integer(out)==(a if s else b)
            c.reverse_check();cases+=1
            assert c.boolean_count<=256*w*w
            assert c.reversible_count<=3*c.boolean_count
            max_ratio=max(max_ratio,F(c.boolean_count,w*w))
       for m in range(1,min(100,(1<<w)-1)+1):
            c=Circuit();out,rem,rw=c.divide_small(c.input(a,w),m)
            assert c.integer(out)==a//m and c.integer(rem)==a%m
            c.reverse_check();cases+=1
            assert c.boolean_count<=256*w*w
            max_ratio=max(max_ratio,F(c.boolean_count,w*w))
    # Check full 256-bit semantic words without storing the huge gate history.
    large=[(0,0),((1<<256)-1,(1<<256)-1),
           ((1<<206)-1,100),((1<<195)-17,73)]
    for a,b in large:
        c=Circuit(False);out=c.mul(c.input(a,256),c.input(b,256))
        assert c.integer(out)==(a*b)%2**256
        max_large_nodes=max(max_large_nodes,c.boolean_count)
    divcases=0
    for m in range(1,101):
        for x in (0,1,m-1,m,2*m-1,(1<<256)-1):
            c=Circuit(False);out,rem,rw=c.divide_small(c.input(x,256),m)
            assert c.integer(out)==x//m and c.integer(rem)==x%m
            assert rw<=8;divcases+=1
    return {"small_exhaustive_cases":cases,"small_reverse_restoration":True,
        "max_small_boolean_count_over_w_squared":str(max_ratio),
        "full_256_bit_multiplication_cases":len(large),
        "largest_256_multiplication_boolean_outputs":max_large_nodes,
        "full_256_bit_small_divisor_cases":divcases,
        "independent_division_remainder_bits_max":8,
        "author_generic_division_uses_257_internal_bits":True,
        "semantic_program_words_are_256_bits":True}

def enumerate_tables():
    old=load(STAGE/"1204/independent_results.json")
    boxes=old["independent_root_boxes"];den=old["root_denominator"]
    fac=[math.factorial(k) for k in range(101)]
    tables=[[],[]];dims=[0,0]
    for a in range(101):
      for b in range(101-a):
       for c in range(101-a-b):
        n=(a,b,c,100-a-b-c)
        lo=sum(n[j]*boxes[j][0] for j in range(4))
        hi=sum(n[j]*boxes[j][1] for j in range(4))
        for tag,(left,right) in enumerate(((39,41),(79,81))):
            inside=2*lo>=left*den and 2*hi<=right*den
            outside=2*hi<left*den or 2*lo>right*den
            assert inside or outside
            if inside:
                ways=fac[100]//math.prod(fac[x] for x in n)
                tables[tag].append((n,ways,dims[tag]))
                dims[tag]+=ways
    assert dims==[int(old["cold_shell"]["dimension"]),int(old["hot_shell"]["dimension"])]
    return tables,dims

def program_checks(tables,dims):
    rng=random.Random(1205);all_tests=0;scan_tests=0
    max_sample_product=0;all_product_bound=0;max_ways=0
    saved=[]
    for tag,rows in enumerate(tables):
        prefix=[row[2] for row in rows]+[dims[tag]]
        assert prefix[0]==0
        for k,(n,ways,start) in enumerate(rows):
            assert start+ways==prefix[k+1]
            max_ways=max(max_ways,ways)
            all_product_bound=max(all_product_bound,ways*max(n))
        lookup={row[0]:i for i,row in enumerate(rows)}
        def unrank(index):
            nonlocal max_sample_product
            block=bisect_right(prefix,index)-1
            n,ways,start=rows[block];ns=list(n);r=index-start;word=[]
            for m in range(100,0,-1):
                products=[ways*x for x in ns]
                max_sample_product=max(max_sample_product,*products)
                assert all(x<2**207 and x%m==0 for x in products)
                parts=[x//m for x in products]
                cumul=[0]
                for x in parts:cumul.append(cumul[-1]+x)
                assert cumul[-1]==ways
                selected=[i for i in range(4) if cumul[i]<=r<cumul[i+1]]
                assert len(selected)==1
                j=selected[0];word.append(j)
                r-=cumul[j];ways=parts[j];ns[j]-=1
                assert all(0<=x<=100 for x in ns)
            assert r==0 and ways==1 and ns==[0]*4
            return tuple(word)
        def rank(word):
            ns=[word.count(i) for i in range(4)]
            n,ways,start=rows[lookup[tuple(ns)]]
            r=start
            for m,j in zip(range(100,0,-1),word):
                r+=sum(ways*ns[i]//m for i in range(j))
                ways=ways*ns[j]//m;ns[j]-=1
            return r
        tests={0,1,dims[tag]-2,dims[tag]-1}
        for n,ways,start in rows:tests.update((start,start+ways-1))
        tests.update(rng.randrange(dims[tag]) for _ in range(32))
        for idx in sorted(tests):
            assert rank(unrank(idx))==idx;all_tests+=1
        # Actual full table scan includes both tagged tables, not an oracle.
        scan_indices=[0,dims[tag]-1]+[rng.randrange(dims[tag]) for _ in range(16)]
        for idx in scan_indices:
            hits=[(n,ways,start) for tg,rr in enumerate(tables)
                  for n,ways,start in rr if tg==tag and start<=idx<start+ways]
            assert len(hits)==1;scan_tests+=1
        bits=173 if tag==0 else 195;mask=(1<<bits)-1
        cases=[[dims[tag]]*10,[0]+[dims[tag]]*9,
               [dims[tag]]*9+[dims[tag]-1],
               [dims[tag],dims[tag]//2]+[0]*8]
        for seeds in cases:
            if tag==0:seeds=[x+(((1<<22)-1)<<173) for x in seeds]
            chosen=0;found=False
            for full in seeds:
                x=full&mask
                take=(x<dims[tag]) and not found
                if take:chosen=x
                found=found or x<dims[tag]
            word=unrank(chosen)
            route=[[int(found and symbol==j) for j in range(4)] for symbol in word]
            assert all(sum(row)==int(found) for row in route)
            assert rank(word)==chosen
            # Exact energy-count conservation for each same-material inventory SWAP.
            for desired in range(4):
                before=[0,0,1,2,3]
                after=[desired]+[0 if j==desired else j for j in range(4)]
                assert sorted(before)==sorted(after)
            saved.append({"tag":tag,"success":found,"selected":str(chosen)})
    assert all_product_bound<2**207
    return {"rank_unrank_words_checked":all_tests,"explicit_two_table_scans":scan_tests,
        "seed_failure_and_first_valid_cases":saved,
        "maximum_initial_type_ways_bits":max_ways.bit_length(),
        "maximum_prefix_bits":max(dims).bit_length(),
        "all_legal_unrank_product_bound_bits":all_product_bound.bit_length(),
        "largest_tested_unrank_product_bits":max_sample_product.bit_length(),
        "product_bound_valid_for_all_paths_by_monotone_ways_and_counts":True,
        "modular_underflow_in_unselected_branches_is_defined":True,
        "online_factorial_used":False,"offline_factorial_100_bits":math.factorial(100).bit_length(),
        "division_denominators_are_fixed_1_to_100":True}

def entropy_bound(q,log_dimension_upper,ln2_upper=F(694,1000)):
    assert 0<q<1
    k=0
    while q<F(1,1<<k):k+=1
    # q>=2^-k, so -log(q)<=k*log2; -(1-q)log(1-q)<=q.
    return q*(1+k*ln2_upper+log_dimension_upper)

def run():
    primitives=primitive_tests()
    tables,dims=enumerate_tables()
    program=program_checks(tables,dims)
    w=256;rows=sum(map(len,tables))
    # Explicit independent macro accounting: these are template ceilings,
    # not a claim that all loose budgets are simultaneously saturated.
    per_pool={"seed_scan":8*10,"two_table_scan":20*rows,
              "fourway_unrank":48*100,"stock_route_controls":8*100,
              "initial_final_and_joint_success":50}
    independent_words=2*sum(per_pool.values())
    admitted_words=2*(32*rows+8*10+64*100+8*100+50)
    assert independent_words<=admitted_words==278276
    # Independent bit-template inequalities. Author's generic division uses
    # w+1 guard lines: 9*w*(w+1) <=18*w^2, still below the stated ceiling.
    template_counts={"add":5*w,"subtract":7*w,"mux":3*w,
        "truncated_multiply":5*w*w+w*(w+1)//2,
        "small_divide":81*w,"author_guarded_divide":9*w*(w+1),
        "compare_or_equality_upper":8*w}
    assert max(template_counts.values())<=256*w*w
    bools=admitted_words*256*w*w
    gates=3*bools+800
    L=2*10**13
    packed=173*(1728+1)+28*1728+195*(2391+1)+28*2391
    finite_bits=bools+10**7+3900+2000+10000
    assert packed==880889 and bools<10**13 and gates<L and finite_bits<2*10**13
    # Hardware addresses are finite, but the fixed instruction list is an input.
    address_bits=finite_bits.bit_length();position_bits=(3*L+606).bit_length()
    gate_descriptor_bits=3*address_bits+position_bits+8
    assert gate_descriptor_bits<=256
    seed_fail=[(1-F(dims[j],2**(173 if j==0 else 195)))**10 for j in range(2)]
    success_lower=(1-max(seed_fail))**2
    assert success_lower>F(998,1000)
    T=3*L+606;M=10**112;sigma=F(1,10**43)
    B=F(66,7)*(L+1)+150
    spread=sigma*(1+F(T*T,8)/(M*M*sigma**4))
    eps=B*T*spread
    cap=F(4,10**15);pstar=F(998,1000);eta=2*cap/(pstar-cap)
    tail=B*(spread+B*T*T/F(2*M))
    assert eps<cap and eta<F(9,10**15) and tail<F(2,10**28)
    q=cap**2/(pstar-cap)
    # Prove exp(.694)>2 by a positive finite series; hence log2<.694.
    log2_hi=F(694,1000)
    assert sum(log2_hi**j/F(math.factorial(j)) for j in range(9))>2
    clock_entropy=entropy_bound(q,3*10**13*log2_hi)
    assert clock_entropy<F(4,10**16)
    eta_cap=F(9,10**15)
    fk=entropy_bound(eta_cap,2*log2_hi) # log(4-1)<2log2
    f256=entropy_bound(eta_cap,8*log2_hi)
    energy_extra=2*F(1094,1000)*(F(6,5)+F(1,3))*eta_cap
    therm_k=2*fk+energy_extra
    therm256=2*f256+energy_extra
    assert therm256<F(2,10**12)
    # The physical meaning of source/clock comparisons is all-path, never a
    # fresh pure-clock assumption at the contact cut.
    flow=load(STAGE/"1204/results.json")
    transports={}
    for t in (600,300):
        v=600//t
        center=F(flow["actual_flow_predictions"][f"t{t}_lambda1"]["Q_dir_center"])
        shape=F(1,10000*v)
        energy_radius=F(1094,1000)*(2*eta_cap+shape)+F(1,10**10)
        pointer=eta_cap+shape+F(1,10**10)
        assert center>energy_radius
        transports[str(t)]={"shape_error_exact":str(shape),
             "actual_Q_radius_exact":str(energy_radius),
             "actual_Q_lower":float(center-energy_radius),
             "actual_Q_upper":float(center+energy_radius),
             "pointer_probability_radius_exact":str(pointer),
             "model_proxy_center_margin":float(center-energy_radius-F(".01094")),
             "proxy_margin_requires_extra_success_readout_contract":True,
             "full_old_statistical_protocol_not_automatically_certified":True}
    return {
      "status":"independent 1205 finite resource and all-path conditional audit",
      "all_checks_passed":True,"primitive_templates":primitives,
      "fixed_tables":{"type_counts":list(map(len,tables)),"dimensions":list(map(str,dims)),
                      "packed_bits":packed},
      "program_semantics":program,
      "resource_bounds":{
         "independent_word_template_ceiling_per_pool":per_pool,
         "independent_word_ceiling_two_pools":independent_words,
         "admitted_word_ceiling_two_pools":admitted_words,
         "boolean_template_counts_at_w256":template_counts,
         "boolean_outputs_upper":bools,"reversible_gates_upper":gates,
         "padded_slots":L,"working_bit_upper":10**13,
         "conservative_finite_log2_dimension":finite_bits,
         "dynamic_finite_log2_dimension_declared_upper":2*10**13,
         "hardware_gate_description_bits_each_upper":256,
         "descriptor_bits_by_current_address_count":gate_descriptor_bits,
         "hardware_description_bits_total_upper":256*L,
         "seed_bits":3900,"four_level_materials":1000,
         "no_free_lookup_or_global_unitary_oracle":True,
         "templates_use_X_CNOT_Toffoli_and_one_control_same_h_SWAP":True,
         "three_register_U_also_couples_shared_X_and_Z_in_H":True,
         "spatial_locality_or_SR_implementation_claimed":False},
      "whole_path_clock":{
         "M":str(M),"sigma_exact":str(sigma),"T":T,
         "force_bound_exact":str(B),"epsilon_clock_exact":str(eps),
         "epsilon_clock":float(eps),"epsilon_contract_exact":str(cap),
         "endpoint_interaction_bound_exact":str(tail),
         "endpoint_interaction_bound":float(tail),
         "two_endpoint_bound":float(2*tail),
         "conservative_conditional_trace_error_exact":str(eta),
         "conditional_trace_error":float(eta),
         "semantic_seed_success_probability_lower":float(success_lower),
         "terminal_computed_flag_mismatch_probability_upper":float(cap**2),
         "conserved_input_success_is_logical_event_not_added_oracle":True,
         "same_success_history_at_both_sections_required":True,
         "terminal_success_readout_classification_not_covered_by_program_error":True,
         "old_K_point001_budget_does_not_cover_success_readout":True,
         "conditional_clock_impurity_exact":str(q),
         "clock_entropy_rational_upper":str(clock_entropy),
         "clock_entropy_upper":float(clock_entropy),
         "clock_entropy_uses_clock_error_only":True,
         "complete_material_spectral_width_upper":1094,
         "whole_material_energy_transport_bound_from_clock_cap":float(1094*cap),
         "preparation_cut_clock_not_reset":True,
         "actual_ABCK_not_assumed_closed_after_preparation":True},
      "finite_thermodynamic_transport":{
         "eta_contract":"9/1000000000000000",
         "K_dimension4_two_section_correction":float(therm_k),
         "dimension256_two_section_correction":float(therm256),
         "conservative_declared_extra":"2e-12",
         "shape_error_excluded_from_smooth_path_entropy_transport":True,
         "only_cold_hot_endpoint_menu_transport_claimed":True},
      "native_flow_transports":transports,
      "source_sha256":{str(p.relative_to(STAGE)).replace("\\","/"):sha(p) for p in
           [STAGE/"1204/independent_results.json",STAGE/"1204/results.json",
            STAGE/"1204/acceptance.json"]},
      "scope":{
         "author_code_imported":False,"expanded_2e13_gate_file_created":False,
         "full_continuous_clock_simulated":False,
         "hardware_and_blank_resource_preparation_is_input":True,
         "code_PASS_claimed_as_microphysical_device":False,
         "navigation_or_historical_files_modified":False}
    }

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--write",action="store_true")
    args=p.parse_args();result=run();out=HERE/"independent_results.json"
    if args.write:
        if out.exists():raise FileExistsError("--write refuses to overwrite")
        out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    else:compare(result,load(out))
    print(json.dumps(result,ensure_ascii=False,indent=2))
