# 1275：完整热参考的两时共同校准及测量能源账

2026-10-17。依据[准入](admission.md)，作者核验1组。本轮是[1237原完整载体](../1237/proof.md)及[1274完整Gibbs参考](../1274/proof.md)的新有限动态推论，不是新的普遍KMS定理。本稿将完整参考和两次理想Lüders读取明确作为数学测试合同；未声称其实际准备、原中途读口和装置来源已获证。

## 1. 固定对象、参考与所比较的时间

保留原全部边界$B$、工作$W$、720维旗标$F$及三站程序$C$。令
$$
P=P_4\otimes I_{WF},\quad H_B=4eP_4+6eP_6,\quad
H_W=2e|1\rangle\langle1|_W,\quad H_0=H_B+H_W.
$$
完整生成元仍为
$$
H_{\rm tot}=I_C\otimes H_0+
g\bigl(|1\rangle\langle0|\otimes U_1+
|2\rangle\langle1|\otimes U_2+\mathrm{h.c.}\bigr),\qquad
e,g>0,\quad 4e>\sqrt2g.
\tag{1}
$$
未丢弃相互作用或旗标。原路径装饰酉
$$
\mathcal D=\operatorname{diag}(I,U_1,U_2U_1)
$$
满足
$$
\mathcal D^\dagger H_{\rm tot}\mathcal D
=gA_3\otimes I+I_C\otimes H_0.
\tag{2}
$$
所以正常、忠实的有限维参考精确为
$$
\mathcal D^\dagger r_\beta\mathcal D
=r_{C,\beta}\otimes r_{H_0F,\beta},\quad
r_{C,\beta}=\frac{e^{-\beta gA_3}}{Z_C},\quad
r_{H_0F,\beta}=\frac{e^{-\beta H_0}}{Z_0}\otimes\frac{I_F}{720}.
\tag{3}
$$
这是完整相互作用的Gibbs参考，绝非将裸参考直接替换全参考。1274已证其$B$边缘为$e^{-\beta H_B}/Z_B$。特别地
$$
\beta_*=\frac{\log(5/2)}{2e},\qquad
r_{B,*}=\frac{P_4}{4}+\frac{P_6}{10},\qquad
-\log r_{B,*}=L+\log2\,I.
\tag{4}
$$
本轮比较同一父参数$t$下的两种预测：完整式(1)演化，以及孤立边缘$H_B$演化。后者也代表$L$模流经固定参数换算后的预测，因为两者均保持$P_4$。不因此认证父参数为实际观察者固有时，也不认定所有区域的模流本就应等于完整物理时间。完整$r_\beta$的Gibbs/KMS身份没有被否定。

## 2. 同能源伙伴块与实际扇区投影

令$S$为$\mathcal H_4\otimes|1\rangle_W$上的投影，$T$为$\mathcal W\mathcal H_4\otimes|0\rangle_W$上的投影，均张量$I_F$。它们维数相同、正交，裸能源均为$6e$；$U_1$互换这两个子空间。原$U_2$保持$P$，因此
$$
Q:=U_1^\dagger P U_1=P-S+T,\quad \Delta=S-T,\quad
E_0=|0\rangle\langle0|_C,
$$
$$
\widetilde P:=\mathcal D^\dagger(I_C\otimes P)\mathcal D
=I_C\otimes Q+E_0\otimes\Delta.
\tag{5}
$$
在$S$上$\widetilde P=E_0$，在$T$上$\widetilde P=I_C-E_0$；其余数据子空间上$\widetilde P$恒为0或1。完整参考中的两块等权：
$$
w:=\operatorname{Tr}(r_{H_0F,\beta}S)
=\operatorname{Tr}(r_{H_0F,\beta}T)
=\frac{2e^{-6\beta e}}
{(2e^{-4\beta e}+5e^{-6\beta e})(1+e^{-2\beta e})}.
\tag{6}
$$
定义边缘扇区概率
$$
p_B=\frac{2e^{-4\beta e}}{2e^{-4\beta e}+5e^{-6\beta e}}.
$$
在$\beta_*$时，$p_B=1/2$且$w=1/7$。等权是同一完整热参考与原保能交换的结果；不是另外指定两个不相干准备。

## 3. 精确动态关联与完整详细平衡

记
$$
\omega=\sqrt2g,\quad z=\beta\omega,\quad Z_C=1+2\cosh z,\quad
a_C=\operatorname{Tr}(r_CE_0)=\frac{1+\cosh z}{2Z_C},
$$
$$
c(t)=\operatorname{Tr}(r_C E_0(t)E_0),\qquad
E_0(t)=e^{igA_3t}E_0e^{-igA_3t}.
$$
完整非对称关联为
$$
C_P(t):=\operatorname{Tr}(r_\beta P(t)P)
=p_B-2wa_C+2wc(t).
\tag{7}
$$
三站链的能源为$+\omega,0,-\omega$，$E_0$对这三条本征向量的重叠概率为$1/4,1/2,1/4$。按约定$c(t)=\sum_\nu c_\nu e^{i\nu t}$，直接得到
$$
c_0=\frac{\cosh z/8+1/4}{Z_C},\quad
c_{+\omega}=\frac{1+e^{-z}}{8Z_C},\quad
c_{-\omega}=\frac{1+e^{z}}{8Z_C},
$$
$$
c_{+2\omega}=\frac{e^{-z}}{16Z_C},\qquad
c_{-2\omega}=\frac{e^{z}}{16Z_C}.
\tag{8}
$$
故正负频率之比为$e^{-\beta\omega}$和$e^{-2\beta\omega}$；这是所选傅里叶符号下的完整Gibbs详细平衡，不是温度选择结果。非零频率在孤立$H_B$对$P_4$的预测中不存在。

取有限时刻
$$
t_m=\frac{\pi}{2\sqrt2g}.
$$
式(8)给
$$
\operatorname{Re}c(t_m)=\frac1{4Z_C},\qquad
a_C-\operatorname{Re}c(t_m)=\frac14.
\tag{9}
$$
因此$C_P(0)-\operatorname{Re}C_P(t_m)=w/2$。这里不以非对称关联直接冒充两次读数的概率；下一节显式包含测量反作用。

## 4. 两次理想Lüders读取的全部联合记录

令$P_1=I_C\otimes P$、$P_0=I-P_1$，$U(t)=e^{-itH_{\rm tot}}$。本测试合同的完整、保全部结果仪器为
$$
\mathcal J_t(\rho)=\sum_{\alpha,\gamma=0}^1
|\alpha,\gamma\rangle\langle\alpha,\gamma|\otimes
P_\gamma U(t)P_\alpha\rho P_\alpha U(t)^\dagger P_\gamma.
\tag{10}
$$
对任意旁挂伙伴$R$以$\mathcal J_t\otimes{\rm id}_R$定义：仪器不直接作用$R$，非选择输出的$R$边缘保持；系统及其与$R$的关联可受测量反作用，条件支的$R$边缘也可改变。保留全部结果及相应条件后态，不声称未知输入完整可逆恢复。经典化的Lüders记录并不保留跨结果相干；另行保留相干指针的全等距实现不由式(10)认证。本节概率针对系统边缘恰为$r_\beta$的输入及其任意扩张，不推广为所有输入同一数表。完整后态并未保持$r_\beta$。

令$u_{00}(t)=\langle0|e^{-igA_3t}|0\rangle$。在式(5)的$S$块上，两次结果1的贡献为$wa_C|u_{00}|^2$；在$T$块上为
$$
w\bigl[1-2a_C+2\operatorname{Re}c(t)-a_C|u_{00}|^2\bigr].
$$
其余恒为1的数据子空间权重是$p_B-w$。相加得到
$$
p_{11}(t)=p_B-2wa_C+2w\operatorname{Re}c(t).
\tag{11}
$$
首次非选择读取在$S,T$上均实施同一个$E_0$去相干通道；两块权重相等、末效果互补，所以它们对末边缘的反作用贡献相消。两次记录的各自边缘均为$p_B$，对所有$t$成立。因此
$$
p_{10}(t)=p_{01}(t)=p_B-p_{11}(t),\qquad
p_{00}(t)=1-2p_B+p_{11}(t).
\tag{12}
$$
式(11)恰等于$\operatorname{Re}C_P(t)$，是本模型等能源、等权互补块的特殊恒等式；不能当作一般Lüders关联规则，也不表示非选择测量没有物理影响。

在$\beta_*$及$t_m$，完整联合表为
$$
\boxed{\begin{array}{c|cc}
 &\gamma=0&\gamma=1\\ \hline
\alpha=0&3/7&1/14\\
\alpha=1&1/14&3/7
\end{array}}
\tag{13}
$$
孤立$H_B$或$L$流保持$P_4$，在同一边缘参考、相同理想读取规则下给$p_{00}=p_{11}=1/2$、$p_{01}=p_{10}=0$。全部四结果分布的总变差为
$$
\boxed{D_{\rm TV}=\frac17.}
\tag{14}
$$
“结果发生翻转”这一联合事件已有概率$1/7$对0的区别，无须删去任何分支。此为指定数学仪器的有限差，不是原1237已授的中途物理读口。

## 5. 首次读取的能源改变必须入账

首次非选择Lüders通道为$\mathcal M(\rho)=\sum_\alpha P_\alpha\rho P_\alpha$。因$[P_\alpha,H_0]=0$，它不改变裸$H_B+H_W$期望。但在每个$S,T$块中，它删除程序链0与1间的热关联。由
$$
(r_C)_{01}=-\frac{\sinh z}{\sqrt2 Z_C}
$$
直接得到
$$
\boxed{\Delta E_{\rm first}
=\operatorname{Tr}\!\left[H_{\rm tot}(\mathcal M(r_\beta)-r_\beta)\right]
=\frac{2w\sqrt2g\,\sinh z}{1+2\cosh z}>0.}
\tag{15}
$$
在$e=g=1,\beta=\beta_*$时约为$0.08166756$。式(15)已经足以否决“读数边缘不变，所以读取免费或保持全部能源”的解释；不是某个实现的最小功定理。若兑现式(10)，必须登记仪器、控制、两份记录和全部来源，并核后一次读取及等待的能源账。这里没有这些实现证书，也不声称原自治过程已经自行实施两读。平衡参考在不插入读取时保持平衡；插入读取后全态、关联和相互作用期望改变，与式(12)的两边缘相同并不矛盾。

## 6. 有限误差与结论边界

若两份完整联合分布分别与上述全父预测和孤立边缘预测总变差至多$\epsilon_{\rm actual},\epsilon_{\rm bare}$，则三角不等式给
$$
D_{\rm TV}(\widetilde p_{\rm actual},\widetilde p_{\rm bare})
\ge\frac17-\epsilon_{\rm actual}-\epsilon_{\rm bare}.
\tag{16}
$$
这是事先给定的分布误差合同。它不能仅由各次单时边缘误差领取；实际准备、读取、时序和记录所需误差须另行推导。在两项各为$1/1000$的条件例中，下界为$493/3500>0.14$。

可选的纯数学时序护栏是$D_{\rm TV}(p(t+\delta t),p(t))\le\min(1,\sqrt2g|\delta t|)$：$H_0$与相互作用$K=H_{\rm tot}-H_0$及末效果对易，故两读概率只依赖$e^{-itK}$，而$\|K\|=\sqrt2g$；酉Duhamel界和测量收缩给该式。它不是实际钟误差认证，也不把完整输出通道中的$H_0$相位删掉。

本轮精确地排除：仅凭式(4)的静态边缘模身份，就让孤立$L/H_B$流复现同一父时间下全部边界接续记录。没有排除全父Gibbs/KMS、区域模理论、物理观察者温度或用户共同治理猜想。原$a,e,g$及时间归一未被选择；独立Casimir仍未认证为Lorentz面积。本结果不生成几何、Clausius或Einstein，也不拼接1270/1272的时钟和来源权限。

1237及1274给非平衡轨迹与静态参考；本轮新量是同一完整平衡参考中的有限两读联合记录及首次读取能源变化，不重报一般对易子、KMS或仪器运输定理。完成后停止此载体扩建。

## 7. 复算与来源

[check.py](check.py)只读复用1237的144维Gauss载体、实际$\mathcal W$及全部720置换，再逐旗标重建42维完整Hamiltonian。直接矩阵Gibbs、完整Lüders顺序和能源变化与式(6)—(15)比较，覆盖多个$\beta$及有限时间；另用精确分数核式(9)、(13)、(14)、(16)。连续$t$及全部$\beta$量词由本稿解析证明承担，浮点样本不是其替代。

默认执行仅比对[results.json](results.json)；仅首次显式--write排他创建。只计本轮作者1组，不把继承载体或720块拆成新组。
