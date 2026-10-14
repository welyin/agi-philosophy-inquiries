# 微正则参考准备预筛：固定能窗、记录与原接触

2026-10-14。**未编号准入记录，不是新轮验收。** 本页仅筛查一个固定准备，不改目标，不调整N或能窗以追求热谱匹配。计数代码已在内存实际运行；整数能窗归属已核，热拟合与接触数值为复用1201结果的诊断，未另签新实验误差合同。

## 1. 历史去重与准备对象

相关旧入口：

|旧结果|可复用范围|尚非本准备|
|---|---|---|
|[973 §3—4](../../archive_956_989/research_note_973.md)|已给Gibbs模式的总能量分布及逐能源壳诊断|未从独立能窗准备原材料参考|
|[977 §3](../../archive_956_989/research_note_977.md)|能量分箱与量子块的共同有效表示|明确不是实际谱测量或后选制备|
|[637 §7](../../archive_629_652/research_note_637.md)、[875 §3—5](../../archive_873_893/research_note_875.md)|已给Gibbs的谱截断及过程/熵控制|不产生独立参考|
|[980](../../archive_956_989/research_note_980.md)|原材料、完整QQ接触、有限环境和开关账|初始热副本直接输入Gibbs|
|[另一路准备预筛](reference_preparation_screen.md)|三材料产品准备和有限交换诊断|没有本次微正则准备链|

按全项目微正则、典型性、能窗/能壳等关键词及相关旧轮回查，未找到“980原材料的均匀准备→独立总能窗→单/双边缘→原QQ联合预测”旧证书；这是检索结论，不是绝对不存在证明。微正则组合计数、近似热谱或拟合一个β本身均不作为新原理。

本次明确选择**均匀非相干能级准备**，每份边缘为$I/4$，各份初始独立。固定原980的$h$，接触关闭时$H_N=\sum_a h_a$。先给定

$$
N=100,\qquad u=0.2,\qquad W=[Nu-1/2,Nu+1/2]=[19.5,20.5].
$$

宽度1是总能窗全宽，不是每份宽度；这些数值未从目标β反选。记$P_W=\mathbf1_W(H_N)$，理想成功分支为

$$
\rho_W=P_W/D_W,\qquad p_{\rm acc}=D_W/4^N.
$$

成功与失败的标志及全部准备/测量记录都留在整体中，不能把归一化成功态当作无条件免费通道。

**等幅纯态不能混用：** 对$|+\rangle^{\otimes N}$只测粗能窗，成功态一般仍是跨不同总能量的相干纯态，未必平稳；即使测精确总能量并保留记录，等能量简并内的排列相干仍存在。例如两份的$|01\rangle+|10\rangle$与对应经典混合有相同单边缘但不同双边缘。本页以下双边缘只适用于所声明的均匀非相干准备。

## 2. 有限组合计数及实际代价

对占据数$n_0+\cdots+n_3=N$，设

$$
m(n)=\frac{N!}{\prod_i n_i!},\quad
D_W=\sum_{n:\,n\cdot e\in W}m(n),
$$

$$
r_i=\frac{\sum_Wm(n)n_i}{ND_W},\qquad
R_{ij}=\frac{\sum_Wm(n)n_i(n_j-\delta_{ij})}{N(N-1)D_W}.
$$

只枚举$\binom{103}{3}=176851$种占据型，不建立$4^{100}$维矩阵。原三个非中心能级满足

$$
32000x^3-31200x^2-2095x-17=0.
$$

用整数二分隔离每根至分母$10^{18}$的有理区间，再逐占据型做整数能窗比较。全部归属明确；最近边界的保守距离大于$1.83859\times10^{-5}$。

结果：

- 选中1728种占据型，$D_W=9399690799043690121938015947926455391852712761214920$。
- $p_{\rm acc}=5.84944194496\times10^{-9}$。
- 每成功池平均约$1.70956\times10^8$次尝试、$1.70956\times10^{10}$份材料。平均数不等于有限最坏预算；若硬截为$5.2\times10^8$次尝试，成功概率约$.952247$，材料上限$5.2\times10^{10}$份，仍保留失败标志。
- 输入总裸能量均值为49.375，成功池均值为19.969725。理想能量测量与$H_N$对易，选择成功分支并不意味着无条件免费降能；失败材料、记录及装置仍在总账。
- 成功边缘熵$\log D_W\simeq119.672517$自然单位；成功的惊奇量$-\log p_{\rm acc}\simeq18.956920$不是已证明的提功或擦除功。

本次没有构造巨大的总能量投影装置。若均匀准备已经保留每份可信能源标签，能窗成功标志是标签串能量和的函数，在该明确态族上可由记录计算给出同一分支统计；能源基准备、标签真实性、比较计算、时序和控制仍须独立说明。它不是等幅纯态的一般替代，更不是已完成内部仪器的声明。

## 3. 单材料边缘与温度检验的界

计数给

$$
r=(.392125469729,\ .387527434887,\ .110847426422,\ .109499668962),
\qquad e\cdot r=.1996972500.
$$

仅用此均能唯一拟合得到$\beta_{\rm mean}=1.19991273521$。相应Gibbs人口约

$$
\tau_{\rm mean}=(.400045238502,\ .379193898469,\ .112931715213,\ .107829147816).
$$

最大人口偏差$.00833353642$。三个名义能隙人口比给逆温度约$(1.19861572,1.23990202,.31747934)$。这些是有限壳偏差诊断；**不能仅由三隙差或相对于均能拟合点的偏差，宣布1201排除共同温度。** 当前1201人口预算为$.015701$，大于上述单点人口偏差；均能拟合点也不等于到全部Gibbs族的最小距离。是否能在完整观测误差下拒绝全β族，须另作实际距离检验。

复用[1201原接触结果](../1201/independent_results_1201.json)，$g=.1$、A仍准备中心态：

|输出|真实有限壳|均能匹配Gibbs|
|---|---|---|
|600时长A四项概率|(.375998163,.195379771,.294749306,.133872760)|(.383575928,.191196593,.292270900,.132956579)|
|300时长A四项概率|(.157872075,.057838726,.735485409,.048803791)|(.161045517,.056610363,.733967737,.048376383)|
|300时长$\Delta E_A$|−.222913897690|−.225035419219|
|300时长$\Delta E_B$|.223674841214|.225808314756|
|300时长关断均值功|.000760943524|.000772895537|

600的最大概率偏差约$.007577765$仍只是对该拟合点的差；不能代替针对全部β与实际误差的证书。这里温度是可检查的模型参数，并未从“微正则”这个名称获得实际测温身份。

## 4. 双材料关联与留出联合任务

同一成功池的两份不同材料具有对角双边缘

$$
R\simeq
\begin{pmatrix}
.153404943640&.151408616191&.043956053448&.043355856450\\
.151408616191&.149904706338&.043327660459&.042886451899\\
.043956053448&.043327660459&.011869930263&.011693782253\\
.043355856450&.042886451899&.011693782253&.011563578360
\end{pmatrix}.
$$

它与$r\otimes r$的总变差约$.00346353863$。对两份分别用独立中心探针做原600、300接触，保留两次记录，精确的操作连接是

$$
p_{jk}^{\rm joint}=(M_{600}R M_{300}^{T})_{jk}.
$$

这与误换成独立同边缘材料的预测$(M_{600}r)_j(M_{300}r)_k$总变差约$.000695029271$。该非零值不自动表示在当前统计合同下可分辨；下一步若采用，必须事前给联合读数样本量、误差及其高昂准备成本。

条件于第一读数$j$，第二源人口为

$$
r_k^{(j)}=\frac{\sum_i(M_{600})_{ji}R_{ik}}{(M_{600}r)_j}.
$$

其均能依$j=0,1,2,3$分别约$.2021157763,.2019483739,.1988926645,.1913905975$，从而给可复算的第二接触预测。**这是保留记录的条件更新，不能冒称未读取第一结果时的远端源自动冷却。** 对第一材料的非选择操作不改变第二材料边缘；真正重复使用已碰撞材料还需保其后态及新关联。

池内材料不是独立新鲜副本。若统计合同要求独立联合试验，可每次使用独立成功池中的一对，或另证适用于相关池内抽样的界；不能直接套独立样本公式。若R保留最初所有能级标签，成功后B池与R仍相关，给定完整标签时源状态甚至已确定。局部边缘温度性质不等于与R解关联，也不涵盖访问这些标签后的全部许可未来。

## 5. 方向判断与停止线

这条准备规则不输入目标β，且单/双边缘、后选代价和原接触能源/联合记录能在同一有限对象上计算，具备进一步审查的接口。当前最有价值的是**独立准备的资源账与保关联的实际未来预测**，不是再画一条接近Gibbs的谱。

成功率极低；总能窗实现、标签制备、失败材料处理、全部记录、时钟控制及实际区域身份尚未共同实现。本页不据此否定一般微正则准备，也不自动扩大N、改窗或把构造该装置列成新阶段门槛。

## 6. 内嵌复算代码

在项目根目录，用已装Python与NumPy执行以下代码；只读取1201独立结果，不写文件。能窗计数使用整数根区间，β和接触联接使用浮点；后者不替代新联合任务的误差认证。

```python
import json, math
from pathlib import Path
import numpy as np
N, D = 100, 10**18
def poly(x):
    return 32000*x**3 - 31200*x*x*D - 2095*x*D*D - 17*D**3
def bracket(a, b):
    a, b = a*(D//10**6), b*(D//10**6)
    sa = poly(a)
    assert sa*poly(b) < 0
    while b-a > 1:
        c = (a+b)//2
        if poly(c)*sa > 0:
            a, sa = c, poly(c)
        else:
            b = c
    return a, b
boxes = [bracket(-54072,-54071), bracket(-9461,-9460),
         (D,D), bracket(1038532,1038533)]
e = np.array([(a+b)/(2*D) for a,b in boxes])
fac = [math.factorial(k) for k in range(N+1)]
Z = types = 0
one, two = [0]*4, [[0]*4 for _ in range(4)]
margin = 10**100
for a in range(N+1):
    for b in range(N-a+1):
        for c in range(N-a-b+1):
            n = (a,b,c,N-a-b-c)
            lo = sum(n[j]*boxes[j][0] for j in range(4))
            hi = sum(n[j]*boxes[j][1] for j in range(4))
            inside = 2*lo >= 39*D and 2*hi <= 41*D
            outside = 2*hi < 39*D or 2*lo > 41*D
            assert inside or outside
            margin = min(margin, *[abs(2*x-y*D)
                         for x in (lo,hi) for y in (39,41)])
            if not inside:
                continue
            w = fac[N] // math.prod(fac[x] for x in n)
            Z += w
            types += 1
            for i in range(4):
                one[i] += w*n[i]
                for j in range(4):
                    two[i][j] += w*n[i]*(n[j]-(i==j))
r = np.array([x/(N*Z) for x in one])
pair = np.array([[x/(N*(N-1)*Z) for x in row] for row in two])
assert max(abs(pair.sum(axis=0)-r)) < 1e-14
p, mean = Z/4**N, float(e@r)
def thermal(beta):
    v = np.exp(-beta*(e-e[0]))
    return v/v.sum()
lo, hi = -50., 50.
for _ in range(100):
    mid = (lo+hi)/2
    if e@thermal(mid) > mean:
        lo = mid
    else:
        hi = mid
beta = (lo+hi)/2
tau = thermal(beta)
path = Path("research_cognition_physics/archive_1201_/1201/independent_results_1201.json")
old = json.loads(path.read_text(encoding="utf-8"))
M = np.array(old["population_readout"]["matrix_rows_A_output_columns_B_population"])
L = np.array(old["predeclared_holdout"]["matrix"])
joint, product = M@pair@L.T, np.outer(M@r,L@r)
coeff = old["predeclared_holdout"]["energy_coefficients_by_initial_B_population"]
conditional = []
for j in range(4):
    post = (M[j]@pair)/(M[j]@r)
    conditional.append(dict(first_outcome=j, second_populations=post.tolist(),
                           second_mean=float(e@post),
                           second_300_probability=(L@post).tolist()))
out = dict(
    input="uniform incoherent levels, not equal-amplitude coherent preparation",
    N=N, total_energy_window=[19.5,20.5],
    root_boxes=boxes, root_denominator=D,
    boundary_margin_lower=margin/(2*D),
    total_count_types=math.comb(N+3,3), selected_count_types=types,
    shell_dimension=str(Z), probability=p, mean_attempts=1/p,
    mean_materials=N/p, initial_energy_mean=N*sum(e)/4,
    accepted_energy_mean=N*mean, shell_entropy=math.log(Z),
    surprisal_success=-math.log(p),
    hard_attempt_budget=520000000,
    hard_budget_success=-math.expm1(520000000*math.log1p(-p)),
    one=r.tolist(), two=pair.tolist(),
    pair_TV_from_product=float(abs(pair-np.outer(r,r)).sum()/2),
    mean_matched_beta=beta, mean_matched_gibbs=tau.tolist(),
    max_population_deviation=float(max(abs(r-tau))),
    population_beta_triplet=[float(math.log(r[i]/r[2])/(e[2]-e[i]))
                            for i in (0,1,3)],
    actual600=(M@r).tolist(), fitted600=(M@tau).tolist(),
    actual300=(L@r).tolist(), fitted300=(L@tau).tolist(),
    joint_600_300=joint.tolist(),
    joint_TV_from_independent=float(abs(joint-product).sum()/2),
    conditional_second=conditional,
    holdout_energy={k:dict(shell=float(np.array(v)@r),
                          mean_matched_gibbs=float(np.array(v)@tau))
                    for k,v in coeff.items()},
    temperature_rejection_certified=False,
    global_projection_device_implemented=False)
print(json.dumps(out, ensure_ascii=False, indent=2))

```
