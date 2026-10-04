# 706入口：保边界局部通道的双向能量域

2026-10-02。接[705](../../research_note_705.md)，先完成[CAR宇称补充](parity_completion_705.md)，回查[623](../../../archive_585_628/research_note_623.md)、[625](../../../archive_585_628/research_note_625.md)、[637](../../../archive_629_652/research_note_637.md)。[代码](local_domain_entry.py)、[结果](local_domain_entry_results.json)、[核验](entry_checks.json)。本入口不计完整科学轮次。

## 1. 分级重数与固定参照

固定一份637的区域比较能源K_A，整个来源邻域都使用它，不让谱投影或重置向量随来源变化。623已给原H的共同算符域；是否可接真实二阶历史须分别核正向和伴随域。

保留原边界λ，并将重数按本地CAR宇称σ分块。两侧宇称不被强制匹配；σ仅用于选择偶Kraus。取同块最低向量：

$$
M_{A,\lambda}=\bigoplus_{\sigma=\pm}M_{A,\lambda,\sigma},\qquad
K_Ae_{\lambda\sigma j}=a_{\lambda\sigma j}e_{\lambda\sigma j},\quad
K_Ag_{\lambda\sigma}=\mu_{\lambda\sigma}g_{\lambda\sigma}.
\tag{1}
$$

## 2. 同一辅助空间及精确交织

取来源和N无关的数学辅助空间，基为|ok〉及全部|λ,σ,j〉；令bλσj=aλσj−μλσ≥0。使用705同一低能主投影，并仅对a>N的尾加入K_{λσj}=|gλσ〉〈eλσj|⊗I_V：

$$
S_N\psi=P_{A,N}\psi\otimes|ok\rangle
+\sum_{a_{\lambda\sigma j}>N}
K_{\lambda\sigma j}\psi\otimes|\lambda,\sigma,j\rangle,\qquad
S_\infty\psi=\psi\otimes|ok\rangle,\qquad S_N^\dagger S_N=I .
\tag{2}
$$

环境标签对边界群及CAR宇称均取平凡作用，因为每个Kraus自身为规范不变偶算符。它们保存的是分析用的能量标签，不是新带荷补偿装置，也不表示免费物理复位。

以D_E|ok〉=0、D_E|λ,σ,j〉=bλσj|λ,σ,j〉，在包含B的共同空间取

$$
R=K_A+K_B+c\ge1,\qquad
R^+=R_{\rm system}\otimes I+I\otimes D_E .
\tag{3}
$$

P_A,N与R对易，尾输出的系统K_A能量为μ，加环境差值即还原a。因此在相应算符域精确有

$$
R^+S_N=S_NR,\qquad
RS_N^\dagger=S_N^\dagger R^+ .
\tag{4}
$$

这正是625双来源证明需要的反向域合同；只知道通道降低平均能量并不足够。S_N是整体Gauss交织，同时只在A和其数学辅助副本作用。

不同辅助分量正交给出

$$
\|R^+(S_N-S_\infty)\psi\|^2
=2\|R(I-P_{A,N})\psi\|^2\longrightarrow0,\qquad \psi\in D(R).
\tag{5}
$$

S_N†的强收敛另由有限个固定辅助基向量验证：每个固定尾标签最终不再使用，ok分量P_N趋I；再用一致范数≤1延拓。式(4)把它提升为伴随的图范数收敛，不能仅凭S_N强收敛就猜测伴随强收敛。

## 3. 与原H的域连接

637的K_A+K_B在匹配空间经J回到固定K₀。623已证明原T、W及全部相互作用的共同域和图范数等价；固定正系数改变不会改变该域。故在原匹配空间有

$$
D(R)=J D(H_F),\qquad
c_1\|R\psi\|\le\|(JH_FJ^\dagger+c_H)\psi\|
\le c_2\|R\psi\|,\qquad
\|G_a\psi\|+\|C_{ab}\psi\|\le c_{ab}\|R\psi\|.
\tag{6}
$$

充分大固定c_H及正c₁,c₂只用于图范数，非物理能源减除。R包含各节点正束缚势及全部群Casimir；不能只保一个链路的电矩来套式(6)。

这一连接复用623原算符结论，不从637单边形式下界直接推出。它为后续同一真实过程的来源比较准备了所需域；本入口尚未构造新有效Hamiltonian、其自身Gibbs态或完整二阶历史极限。

## 4. 复算

复用705原ν、L、u两切口载体，各重数块现有两偶、两奇向量，声明相对能阶(0,1,7,8)。区域维数56，共同辅助维数13；它们是代数诊断，不是全K_A谱。

N=10、24、36、45均核等距、式(4)正反向、偶Kraus与式(5)。在此对角夹具中相应身份为零残差；式(5)两端分别为2687.40342001、2572.11048162、1887.19987768、0。直接将奇输入回填到偶最低向量会给宇称交换子范数2；正确分块后为0。

一般N→∞与原H共同域由解析承担。下一步将这一来源无关局部等距接入原含时演化与真实仪器，核同一状态／记录／来源的联合系数；不把它叫整体有限维过程，不重算空间维数或705容量反例。目标及四分支范围不变。

