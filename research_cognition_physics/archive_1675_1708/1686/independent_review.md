# 第1686轮独审　两个轻代谱隙怎样加强有限强角的上界

若每部门的有限阈值 $Y_f'=(I+P_f)Y_f$ 只读取完整的 $A=Y_uY_u^\dagger,B=Y_dY_d^\dagger$，沿用主文的实系数与收敛合同，则两个轻代谱隙给出全部阶次的上界 $|\Delta\bar\theta|<1.1367\times10^{-22}$ rad。若阈值能按原算符分别读取同一带电矩阵元的两条路径，取 $\kappa=1$ 的对照给出 $1.95086\times10^{-8}$ rad。两条实际路径的干涉还有一项可通过标量顶角加权和表达的独立读数。[可复算程序](independent_check.py) · [精确有理数上界及矩阵检验](independent_results.json)

## 闭迹与开词的统一界

每个 $P_f$ 是由 $A,B$ 构成的实系数单项之和。每个单项可有一个开词和任意多个归一迹 $\operatorname{tr}_3=\operatorname{Tr}/3$，总词长至少为一，所有系数的绝对值之和至多 1。设 $a=\|A\|_2$、$b=\|B\|_2$，则 $\|P_f\|_2\le\max(a,b)<1$。因而连续从 $P_f=0$ 选出的行列式相位，可写作 $\operatorname{Im}\sum_{k\ge1}(-1)^{k+1}\operatorname{Tr}P_f^k/k$，每阶展开后的系数绝对和仍至多 1。

一个 $\operatorname{Tr}P_f^k$ 中的单项是 3 乘以若干归一迹之积。若全部字母中总数 $n_A<3$ 或 $n_B<3$，则每个迹都实。例如只有两个 $A$ 时，该迹循环等于 $\operatorname{Tr}(A B^p A B^q)$，其复共轭是同一循环词。故只需估计 $n_A,n_B\ge3$。证书既覆盖开放 Yukawa 指标上的乘积，也覆盖全部字母都在闭迹里的标量阈值。

记 $A_0$ 为将 $A$ 的前两个特征值并为 $a_u$ 后的矩阵。$A-A_0=(a_c-a_u)P_c$ 为秩一算子，写 $\Delta a=a_c-a_u$。同理 $B-B_0=\Delta b P_s$。任一矩阵前两代简并时，两矩阵可同时实化：先对角化另一矩阵，再把该简并矩阵剩余唯一谱线的矢量逐分量重定相为实。故在 $(A_0,B),(A,B_0),(A_0,B_0)$ 三种情形，所考察乘积的虚部都为零。在 $A_0+u\Delta a P_c$、$B_0+v\Delta b P_s$ 上作二重积分，得到混合导数的精确表达。共有 $n_A n_B$ 项，每项至少有一个归一迹含秩一导数。用 $|\operatorname{Tr}(P_c X)|\le\|X\|_2$，该迹贡献的 $1/3$ 恰好抵消总的前因子 3。无论两次求导落在同一迹还是不同迹，总有

$$
|\operatorname{Im}(3\prod_j\operatorname{tr}_3 W_j)|
\le n_A n_B\,\Delta a\Delta b\,a^{n_A-1}b^{n_B-1}.
$$

令 $z=4\max(a,b)/3$。对整数 $m\ge3$ 有 $m/3\le(4/3)^{m-3}$。第 $k$ 个对数展开项的总词长至少为 $k$，两个部门合计满足

$$
|\Delta\bar\theta|\le18\,\Delta a\Delta b\,a^2b^2
\left(H_6+\frac{z}{7(1-z)}\right),\qquad H_6=\sum_{j=1}^6\frac1j .
$$

第1397轮的 $10^{16}\,\mathrm{GeV}$ 两宽边界给 $a\le0.207025$、$b\le3.81924\times10^{-5}$、$\Delta a\le2.2801\times10^{-6}$、$\Delta b\le1.7689\times10^{-8}$。复算最后一步用有理数，不依赖某个 CKM 相位或阈值截断。对第1681轮具体非简并点，最低六次 CP 词的绝对虚部可由 $|\det[A,B]|/2$ 算得 $7.3575\times10^{-29}$，与实际 $|J|$ 乘两谱的 Vandermonde 因子一致。该单点数量提供实现核验，全部阶次的结论来自上面的不等式。[第1397轮输入](../../archive_1364_1409/1397/inputs.json) · [Jarlskog 原始不变量](https://doi.org/10.1103/PhysRevLett.55.1039)

## 分路径阈值与可读的标量顶角

第1681轮的 CP 基下，带电下夸克边 $Z$ 有实系数的两项 $z_R,z_D$，味真空幂为 1、5。用原夸克数据恢复有量纲标度前的完整 Yukawa 单位，取总实符号使树级 $\bar\theta=0$。在固定裸强角、相同规范动能并排除额外补偿阈值的比较中，令

$$
Y_d'=Y_d\exp[\kappa(D_R^\dagger D_D-D_D^\dagger D_R)] .
$$

该右酉乘子保持整个 $Y_dY_d^\dagger$，故保持六谱与全部 CKM 读数。行列式给出 $\Delta\bar\theta=2\kappa\operatorname{Im}(z_R^*z_D)\pmod{2\pi}$。这是对指定 Wilson 系数的物理修改。若将同一右酉只用作场变量更换，手征异常要求同时转移胶子拓扑角，物理强角随之保持。若另有重有色物质、独立拓扑耦合或动力学弛豫，则应按完整作用重新计入。

在固定径向参数、其余各边均为单独定模单项的前件下，味真空相位变化时，$\operatorname{Tr}Y_dY_d^\dagger$ 的所有变化都来自 $Z$ 的干涉。设两幂差 $n=4$，得到

$$
\partial_\phi\sum_j y_{d_j}^2=-2n\operatorname{Im}(z_R^*z_D),\qquad
\Delta\bar\theta=-\frac{\kappa}{n}\,\partial_\phi\sum_j y_{d_j}^2 .
$$

若 $G=U_L^\dagger(\partial_\phi Y_d)U_R$ 是与真实质量同基的标量顶角，则右边的导数亦等于 $2\sum_j y_{d_j}\operatorname{Re}G_{jj}$。故不必将两路分别当成直接可探测的外线，也能检验这个单一干涉假设。该读数需一种实际角向标量、已知的规范化与混合投影，且其它矩阵元的模长独立于这一角变量。额外的模长响应或别的路径干涉会进入同一总和。

在冻结点，导数是 $-7.80345\times10^{-8}$，$\kappa=1$ 时强角移动 $1.95086\times10^{-8}$ rad。沿用第1389轮单一强角来源的中子电偶极矩条件界，靠近零的分支要求 $|\kappa|<0.00896$。角向顶角加权和对右酉阈值前后保持相同，因而单凭它不能测得 $\kappa$。若进一步测量完整强角的角向响应，所写阈值在两种方案间还给出 $\partial_\phi\Delta\bar\theta=2\kappa n\operatorname{Re}(z_R^*z_D)=4.50532\times10^{-8}\kappa$；若存在独立角向胶子耦合，须把它共同计入。[第1681轮树级方案](../1681/proof.md) · [第1389轮实验约定](../../archive_1364_1409/1389/results.json)
