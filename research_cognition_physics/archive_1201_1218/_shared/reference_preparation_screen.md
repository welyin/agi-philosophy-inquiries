# 独立有限参考准备：原生交换网络的未编号预筛

> 2026-10-14，未编号、未验收的方向预筛。本文不新增1202或其它研究轮次，不改变目标。解析部分只对明确的共振比较网络成立；数值部分只是事前固定的原64维过程诊断，未经区间误差或实际装置认证。

## 1. 去重：此次不再解决什么

- [305](../../archive_301_341/research_note_305.md)已证明给定满秩参考、物理H及足够方向后，模能量等于$\beta H$当且仅当参考为该H的Gibbs态；平稳性不足。其参考准备没有产生。
- [1020](../../archive_1009_1043/research_note_1020.md)已核完全被动性在给定H及操作族下的参考选择范围；它不选择H、温度，也不构造有限准备。
- [980](../../archive_956_989/research_note_980.md)输入独立Gibbs副本，再用原生接触获得热复位。把这些副本重命名为“内部参考”不能补上独立准备。
- 已有522/529的谐参考、972的热模式、873—893的联合热核等，分别给定参考或Gibbs族；全库针对能壳、微正则、典型性、本征态热化、独立准备的检索，未找到原$h_\eta,Q\otimes Q$从能级产品态产生可复用非零温度参考的完整证书。此检索结论只作为去重记录，不断言文献或所有历史附件中不存在相关思想。

共同资源及读口的范围沿用[1201在研协议审计](../1201/protocol_resource_audit.md)。本预筛不重写全局熵恒等式，不先把准备定义为目标Gibbs态。

## 2. 候选与预先固定的最小检验

候选是同一原生材料的封闭有限链：

$$
H_N=\sum_{j=1}^N h_{\eta,j}
 +g\sum_{j=1}^{N-1}Q_jQ_{j+1}.
$$

初态为已知能级产品态，环境不是另给的热浴。若某一局部材料在事前指定有限窗内成为稳定且通过1201温度留出检验的参考，就可能用一次完整准备接上实际温度记录。其余材料和任何等待时间标签仍保留，不能把局部约化当作环境被复位。

主代理在本组计算前固定：

$$
N=3,\qquad g=1/10,\qquad
|\Psi_0\rangle=|0,2,3\rangle,
\qquad t\in[300,600].
$$

使用原$h_\eta$与$Q$，完整Hilbert空间为$4^3=64$维。不改变能级、耦合、初态或窗口挑选近似Gibbs结果；报告全部三处局部材料。平均窗比较固定为$[300,450]$和$[450,600]$。额外301个整数时点仅是振荡诊断，不是连续时间极值的严格界。

## 3. 共振交换人口守恒：决定性的短证

令$c=2$，$\Pi_i^{(j)}=|i\rangle\langle i|$作用在第$j$份材料，定义各能级总人数

$$
\mathcal N_i=\sum_{j=1}^N\Pi_i^{(j)}.
$$

旧980的共振平均作用在任一相邻对上为

$$
\bar V_{jk}=\sum_{i\ne c}w_i
 \left(|i,c\rangle\langle c,i|+|c,i\rangle\langle i,c|\right)_{jk}.
$$

每个交换项的输入与输出都各有一份$i$和一份$c$，故对全部$r$，

$$
[\mathcal N_r,|i,c\rangle\langle c,i|_{jk}]=0.
$$

裸H对这些人数同样对易，所以在任意有限网络、任何交换时序下，

$$
[\mathcal N_r,\bar H_N]=0,
\qquad \bar H_N=\sum_jh_{\eta,j}+g\sum_{\langle jk\rangle}\bar V_{jk}.
$$

这不依赖图的连通程度、是否长时间混合或是否把更多相同副本加入网络。

**准确结论：** 若要求全部材料最终具有同一局部人口$p_r$，则必有

$$
p_r=\frac{\langle\mathcal N_r\rangle_{\rm in}}{N}.
$$

该人口由初始准备决定，一般不满足单一Gibbs温度。对本组$(0,2,3)$，它是$(1,0,1,1)/3$：有一个精确零人口，其余三个为正，既不是有限温度全秩Gibbs，也不是此非简并H的单一极端能级态。因此原共振交换不能把这份准备转成所有材料共享的Gibbs人口。

本结论不禁止个别局部材料暂时恰为Gibbs而由其它材料补偿人数，也不禁止初始人数已经具有热比例的准备。它反对的是从任意初始人口仅靠这套同种交换自动产生共同温度的推断。

真实$H_N$仍含非共振过程，通常$[H_N,\mathcal N_r]\ne0$；例如$|i,j\rangle\leftrightarrow|c,c\rangle$可改人数。因此上述短证不能直接升级成完整母过程不热化定理。下一项若要有净增，必须检验这些真实过程是否在声明有限域内产生足够大且稳定的参考变化，而非延长理想交换过程后宣称出现Gibbs。

## 4. 固定窗口的完整64维诊断

令$H_N=R\operatorname{diag}(E_\alpha)R^\dagger$，$c=R^\dagger\Psi_0$。无需数值时间积分即可计算有限窗口平均：

$$
\bar\rho_{[a,b]}=R M R^\dagger,
\qquad M_{\alpha\beta}=c_\alpha c_\beta^*
 e^{-i(E_\alpha-E_\beta)(a+b)/2}
 \operatorname{sinc}\!\left((E_\alpha-E_\beta)(b-a)/2\right),
$$

其中$\operatorname{sinc}(x)=\sin x/x$，零点取连续值。这里的平均首先是一个诊断统计量；把它作为实际准备，还需等待时间抽样/时钟及其记录合同。不能用平均态平稳代替每个实际时点稳定，也不能免费删除时间标签。

以原能级基读各局部平均人口。列出的三条$\beta_i$是用于暴露单温度不一致的对数比值：

$$
\beta_i=\frac{\log(p_i/p_c)}{e_c-e_i},\qquad i=0,1,3.
$$

本表只诊断数学人口，不是已实现的有限样本温度估计；小概率的误差不能在实际仪器中省略。

| 站点 | 平均人口$(p_0,p_1,p_2,p_3)$ | 三条边的$\beta_i$ | 前后半窗局部半迹距 |
|---|---|---|---:|
| 0 | $(.9988036853,.0000292914,.0010248512,.0001421721)$ | $(6.52898,-3.52169,51.26263)$ | .000231510 |
| 1 | $(.0005046525,.0004975559,.4957609454,.5032368462)$ | $(-6.53654,-6.83944,-.388430)$ | .0569951 |
| 2 | $(.0002818773,.0000931200,.5035042889,.4961207137)$ | $(-7.10376,-8.51490,.383392)$ | .0569528 |

负的候选$\beta$对有限有界谱不是代数矛盾；本例的问题是同一站点三条边互不一致。不能把负温度本身当作本轮的拒绝理由。

对每站再取与其平均能源相同的Gibbs态作**诊断比较**，人口总变差分别约$.02769697,.05058522,.05924858$。该拟合并非全Gibbs族的总变差最小化，不能把这些数值误称为“离任意Gibbs态的严格下界”。

### 4.1 平均谱掩盖实际振荡

中间站点的$p_2(300)\simeq.62657148$，$p_2(600)\simeq.06528807$；301个预设采样时点的$p_2$最大值减最小值约$.99875932$。末站对应波幅约$.99800419$。因此窗口平均接近二能级各半并不表示可在该窗随时提取同一稳定参考。

平均态本身也不严格与孤立$h_\eta$对易：三个站点的$\|[h_\eta,\bar\rho_j]\|$约为$.00137986,.000889448,.000800743$。它们是原耦合系统的局部约化；关断后不能不加审计就称其为原孤立材料的精确平稳态。

### 4.2 完整能源与人数变化保留

窗口平均总人数为

$$
(\langle\mathcal N_0\rangle,
\langle\mathcal N_1\rangle,
\langle\mathcal N_2\rangle,
\langle\mathcal N_3\rangle)
\simeq(.9995902152,.0006199673,1.0002900854,.9994997320).
$$

相对于初始$(1,0,1,1)$，真实非共振项产生了小而非零的变化。它不能被说成精确人口守恒；它也尚未把这组准备变成共同稳定的单温度参考。

窗口平均裸能源变化为$-.0002131662503451$，平均关断功为$-.0002131662503454$，完整$H_N$能源残差约$2.22\times10^{-16}$。各站裸能源变化约为

$$
(.00123691139,.01835664435,-.01980672199).
$$

因此本诊断没有丢弃接触交互能源。以上浮点结果没有作区间认证；约$-2.7\times10^{-16}$的平均态最小计算特征值属于舍入尺度，不被解释成物理负概率。

## 5. 准入判断

**本组不建议立新编号来继续修这个小链。** 已有精确比较机制说明弱交换保留额外人口自由；固定完整母过程诊断则显示真实人数改变很小，中间两站仍有大振荡，三站都未通过单温度/孤立平稳资格。

本结论只拒绝把当前准备、当前链和当前窗口直接认作已经获得的稳定热参考。它不证明大系统、不同耦合强度、不同初态或其它实际稳定机制均不热化，也不要求微观必须连续或存在无限热浴。

以后再进入“独立参考生成”研究，至少应出现一个不同于重复交换混合的实际机制或有限对象，并在事前说明：它如何改变当前保留的人口信息、如何保留完整能源及准备资源，以及为什么该机制有机会同时满足温度留出和稳定窗口。更长时间或更多副本本身不是新增机制；若引入新相互作用，也不得直接按目标Gibbs比例设置跃迁率再称温度已生成。

## 附录：本组可复算代码

下面代码只依赖现有Python/NumPy和原968材料，运行于项目根目录，只输出诊断，不写文件。完整态计算为64维，没有作RWA或丢弃非共振项。

```python
from pathlib import Path
import importlib.util
import json
import numpy as np

path = Path('research_cognition_physics/archive_956_989/968/internal_relay.py')
spec = importlib.util.spec_from_file_location('old968', path)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
_, h, Q, F, _, _ = old.material()
e, U = np.linalg.eigh(h + .05 * F)
Q = U.T @ Q @ U
I = np.eye(4)
N, g, levels = 3, .1, (0, 2, 3)

def tensor(parts):
    out = np.ones((1, 1))
    for a in parts:
        out = np.kron(out, a)
    return out

def embed(a, j):
    parts = [I] * N
    parts[j] = a
    return tensor(parts)

H0 = sum(embed(np.diag(e), j) for j in range(N))
Qs = [embed(Q, j) for j in range(N)]
V = sum(Qs[j] @ Qs[j + 1] for j in range(N - 1))
H = H0 + g * V
E, R = np.linalg.eigh(H)
psi = np.eye(4**N)[:, np.ravel_multi_index(levels, (4,) * N)]
c = R.conj().T @ psi

def average(a, b):
    d = E[:, None] - E[None, :]
    M = (np.outer(c, c.conj()) * np.exp(-1j * d * (a + b) / 2)
         * np.sinc(d * (b - a) / (2 * np.pi)))
    return R @ M @ R.conj().T

def state(t):
    v = R @ (np.exp(-1j * E * t) * c)
    return np.outer(v, v.conj())

def local(rho, j):
    letters = 'abc'
    expr = letters + letters[:j] + 'd' + letters[j + 1:]
    return np.einsum(expr + '->' + letters[j] + 'd',
                     rho.reshape([4] * (2 * N)))

def distance(a, b):
    return float(np.abs(np.linalg.eigvalsh(a - b)).sum() / 2)

def energy_matched_gibbs(p):
    mean = p @ e
    def gibbs(beta):
        x = -beta * e
        x -= x.max()
        x = np.exp(x)
        return x / x.sum()
    lo, hi = -1000., 1000.
    for _ in range(120):
        mid = (lo + hi) / 2
        if gibbs(mid) @ e > mean:
            lo = mid
        else:
            hi = mid
    beta = (lo + hi) / 2
    return beta, gibbs(beta)

rho = average(300, 600)
halves = [average(300, 450), average(450, 600)]
rows = []
for j in range(N):
    r = local(rho, j)
    p = np.diag(r).real
    beta, thermal = energy_matched_gibbs(p)
    samples = np.array([np.diag(local(state(float(t)), j)).real
                        for t in np.linspace(300, 600, 301)])
    rows.append(dict(
        site=j, populations=p.tolist(),
        energy_change=float(p @ e - e[levels[j]]),
        edge_beta=[float(np.log(p[i] / p[2]) / (e[2] - e[i]))
                   for i in (0, 1, 3)],
        energy_matched_beta=beta,
        matched_gibbs_population_distance=float(np.abs(p-thermal).sum()/2),
        half_window_distance=distance(local(halves[0],j),local(halves[1],j)),
        commutator_norm=float(np.linalg.norm(np.diag(e)@r-r@np.diag(e),2)),
        endpoint_300=np.diag(local(state(300),j)).real.tolist(),
        endpoint_600=np.diag(local(state(600),j)).real.tolist(),
        sampled_ranges=(samples.max(axis=0)-samples.min(axis=0)).tolist()))
counts = []
for k in range(4):
    Nk = sum(embed(np.diag(np.eye(4)[k]), j) for j in range(N))
    counts.append(float(np.trace(Nk @ rho).real))
print(json.dumps(dict(
    model=dict(N=N,g=g,levels=levels,window=[300,600]),
    energies=e.tolist(), sites=rows, mean_counts=counts,
    mean_bare_energy_change=float(np.trace(H0@rho).real-psi@H0@psi),
    mean_switch_off_work=float(-g*np.trace(V@rho).real),
    full_energy_residual=float(np.trace(H@rho).real-psi@H@psi)), indent=2))
```
