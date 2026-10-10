"""Independent exact audit of the original fixed-slot port model.
Standard library only. Does not import the author's construction or results.
The 24-state sector is not asserted to replace the full 17496-state source.
Default: stdout only. --write creates independent_results.json exclusively.
--check reruns and compares the existing independent result without writing.
"""
from fractions import Fraction as F
from itertools import permutations, combinations
from math import factorial
from pathlib import Path
import argparse
import json

def multiply(A, B):
    columns=list(zip(*B))
    return [[sum(x*y for x,y in zip(row,column)) for column in columns] for row in A]

def trace(A):
    return sum(A[i][i] for i in range(len(A)))

def transpose(A):
    return list(map(list, zip(*A)))

def matvec(A, v):
    return [sum(a*b for a,b in zip(row,v)) for row in A]

def parity(values):
    return -1 if sum(values[i]>values[j] for i in range(len(values)) for j in range(i+1,len(values)))%2 else 1

def run():
    assertions=0
    def require(condition):
        nonlocal assertions
        assert condition
        assertions+=1

    leaves=(0,3,4,5)
    assignments=list(permutations(leaves))
    address_index={s:i for i,s in enumerate(assignments)}
    b_subsets=list(combinations(leaves,2))
    graph_index={s:i for i,s in enumerate(b_subsets)}
    # Slot order: b1,b2,c1,c2. The central b-c port is held fixed.
    flip=[[0]*24 for _ in assignments]
    for col,s in enumerate(assignments):
        for i,j in ((0,2),(0,3),(1,2),(1,3)):
            out=list(s)
            out[i],out[j]=out[j],out[i]
            flip[address_index[tuple(out)]][col]+=1
    require(flip==transpose(flip))
    require(all(sum(row)==4 for row in flip))

    graphs=[]
    for ss in b_subsets:
        other=tuple(x for x in leaves if x not in ss)
        edges=tuple(sorted({(1,2)} | {tuple(sorted((1,x))) for x in ss}
                           | {tuple(sorted((2,x))) for x in other}))
        require(len(edges)==5)
        graphs.append(edges)
    initial_edges=((0,2),(1,2),(1,4),(1,5),(2,3))
    initial_graph=graphs.index(initial_edges)
    require(b_subsets[initial_graph]==(4,5))

    compressed={}
    reduced_h={}
    source_columns={}
    spectral={}
    for label in ("uniform","center2_sign"):
        # V/2 is the isometry. In the sign source, owner labels at the
        # central-2 ports are (1,c1,c2), in physical port order.
        V=[[0]*6 for _ in assignments]
        for row,s in enumerate(assignments):
            g=graph_index[tuple(sorted(s[:2]))]
            V[row][g]=1 if label=="uniform" else parity((1,s[2],s[3]))
        gram=multiply(transpose(V),V)
        require(gram==[[4*int(i==j) for j in range(6)] for i in range(6)])
        numer=multiply(transpose(V),multiply(flip,V))
        require(all(x%4==0 for row in numer for x in row))
        A=[[x//4 for x in row] for row in numer]
        require(multiply(flip,V)==multiply(V,A))
        require(A==transpose(A))
        require(all(abs(A[i][j])==int(len(set(b_subsets[i])&set(b_subsets[j]))==1)
                    for i in range(6) for j in range(6)))
        compressed[label]=A
        source_columns[label]=[row[initial_graph] for row in V]
        A2=multiply(A,A)
        A3=multiply(A2,A)
        A4=multiply(A2,A2)
        if label=="uniform":
            require(all(A3[i][j]-2*A2[i][j]-8*A[i][j]==0
                        for i in range(6) for j in range(6)))
            require([trace(M) for M in (A,A2,A3,A4)]==[0,24,48,288])
        else:
            require(A2==[[4*int(i==j) for j in range(6)] for i in range(6)])
            require([trace(M) for M in (A,A2,A3,A4)]==[0,24,0,96])
        return_fourth=F(A2[initial_graph][initial_graph]**2,4)+F(A4[initial_graph][initial_graph],12)
        spectral[label]={"traces":[trace(M) for M in (A,A2,A3,A4)],
                         "pure_flip_return_fourth_coefficient":str(return_fourth)}

        # Verify the two local slot-relabel characters on every column.
        for slots,sgn in (((0,1),1),((2,3),1 if label=="uniform" else -1)):
            for s,row in address_index.items():
                out=list(s); i,j=slots
                out[i],out[j]=out[j],out[i]
                require(V[address_index[tuple(out)]]==[sgn*x for x in V[row]])

        h=[[0]*36 for _ in range(36)]
        for g,edges in enumerate(graphs):
            for d in range(6):
                for gg in range(6):
                    h[6*g+d][6*gg+d]+=A[g][gg]
            for i,j in edges:
                h[6*g+i][6*g+i]-=1
                h[6*g+j][6*g+j]-=1
                h[6*g+i][6*g+j]+=1
                h[6*g+j][6*g+i]+=1
            # Independently check the original six-qubit SWAP sum on
            # each one-excitation column: sum S_e = 5 I - L_G.
            for d in range(6):
                raw=[0]*64
                bit=1<<(5-d)
                for i,j in edges:
                    bi=(bit>>(5-i))&1; bj=(bit>>(5-j))&1
                    target=bit if bi==bj else bit^(1<<(5-i))^(1<<(5-j))
                    raw[target]+=1
                expected=[0]*64
                expected[bit]=5
                for i,j in edges:
                    if d==i:
                        expected[1<<(5-i)]-=1; expected[1<<(5-j)]+=1
                    elif d==j:
                        expected[1<<(5-j)]-=1; expected[1<<(5-i)]+=1
                require(raw==expected)
        require(h==transpose(h))
        require(max(sum(abs(x) for x in row) for row in h)==10)
        reduced_h[label]=h

    require([x*x for x in source_columns["uniform"]]==[x*x for x in source_columns["center2_sign"]])
    require(sum(x*x for x in source_columns["uniform"])==4)
    require(spectral["uniform"]["pure_flip_return_fourth_coefficient"]=="8")
    require(spectral["center2_sign"]["pure_flip_return_fourth_coefficient"]=="16/3")

    t=F(7,10)
    order=40
    # ||h||_2 <= sqrt(||h||_1 ||h||_infty) = 10.
    # The remaining exponential terms have ratios <= 7/42.
    norm_bound=10
    x=t*norm_bound
    vector_tail=x**(order+1)/factorial(order+1)/(1-x/F(order+2))
    probability_tail=vector_tail*(2+vector_tail)
    require(vector_tail>0 and probability_tail<F(1,10**10))
    reports={}
    powers_by_label={}
    receiver_rows=[6*g+5 for g in range(6)]
    initial=6*initial_graph+0
    for label,h in reduced_h.items():
        v=[0]*36; v[initial]=1
        powers=[v]
        for k in range(1,order+1):
            powers.append(matvec(h,powers[-1]))
        powers_by_label[label]=powers
        real=[F(0)]*36; imag=[F(0)]*36
        for k,power in enumerate(powers):
            coefficient=t**k/factorial(k)
            target=real if k%2==0 else imag
            sign=(1,-1,-1,1)[k%4]  # (-i)^k
            for i,val in enumerate(power):
                target[i]+=sign*coefficient*val
        p=sum(real[i]**2+imag[i]**2 for i in receiver_rows)
        norm2=sum(real[i]**2+imag[i]**2 for i in range(36))
        require(abs(norm2-1)<=probability_tail)
        reports[label]={"taylor_probability":str(p),
                        "probability_decimal":float(p),
                        "strict_lower":str(p-probability_tail),
                        "strict_upper":str(p+probability_tail)}
    gap=F(reports["center2_sign"]["taylor_probability"])-F(reports["uniform"]["taylor_probability"])
    gap_lower=gap-2*probability_tail
    require(gap_lower>F(1,400))
    dt=F(1,100000)
    # Each probability derivative is <= 2||h||; their difference <= 40.
    window_lower=gap_lower-4*norm_bound*dt
    require(window_lower>F(1,500))

    # Independent exact short-time coefficients of the original receiver.
    coefficients={}
    for label,powers in powers_by_label.items():
        vals=[]
        for n in range(9):
            total=F(0)
            if n%2==0:
                for k in range(n+1):
                    j=n-k
                    # (-i)^k i^j = (-1)^k i^n, a real sign for even n.
                    sign=(-1)**(k+n//2)
                    inner=sum(powers[k][r]*powers[j][r] for r in receiver_rows)
                    total+=F(sign*inner,factorial(k)*factorial(j))
            vals.append(total)
        coefficients[label]=[str(v) for v in vals]
        require(vals[:6]==[F(0)]*6)
        require(vals[6]==F(7,36))
    require(F(coefficients["uniform"][8])==F(-607,1440))
    require(F(coefficients["center2_sign"][8])==F(-443,1440))
    require(F(coefficients["center2_sign"][8])-F(coefficients["uniform"][8])==F(41,360))

    return {
        "passed":True,
        "assertions":assertions,
        "method":"independent 24-slot integer construction; exact Fraction Taylor40",
        "imports_author_code":False,
        "fixed_slot_address_dimension":24,
        "character_graph_dimension":6,
        "one_excitation_dimension":36,
        "full_17496_source_verified_by_this_algorithm":False,
        "initial_graph_edges":[list(e) for e in initial_edges],
        "source":0,"receiver":5,"t":"7/10",
        "center2_character":"parity of neighbor owner labels (1,c1,c2) in port order",
        "same_initial_graph_and_fixed_block_computational_populations":True,
        "row_norm_bound":norm_bound,
        "vector_tail_bound":str(vector_tail),
        "probability_tail_bound":str(probability_tail),
        "probability_tail_decimal":float(probability_tail),
        "reports":reports,
        "gap_decimal":float(gap),
        "gap_strict_lower":str(gap_lower),
        "gap_greater_than":"1/400",
        "time_radius":str(dt),
        "window_gap_strict_lower":str(window_lower),
        "window_gap_greater_than":"1/500",
        "receiver_coefficients_0_through_8":coefficients,
        "pure_flip_spectral_checks":spectral,
        "port_phase_is_spatial_axis":False,
        "same_report_effect_implies_full_instrument_equality":False
    }

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--write",action="store_true")
    group.add_argument("--check",action="store_true")
    args=parser.parse_args()
    result=run()
    path=Path(__file__).with_name("independent_results.json")
    if args.write:
        with path.open("x",encoding="utf-8",newline="\n") as f:
            json.dump(result,f,ensure_ascii=False,indent=2)
            f.write("\n")
    elif args.check:
        assert result==json.loads(path.read_text(encoding="utf-8"))
    print(json.dumps({"passed":result["passed"],"assertions":result["assertions"],
                      "mode":"write" if args.write else ("check" if args.check else "read_only"),
                      "gap_decimal":result["gap_decimal"],
                      "probability_tail_decimal":result["probability_tail_decimal"],
                      "strict_center_lower":"1/400","strict_window_lower":"1/500"},
                     ensure_ascii=False,indent=2))
