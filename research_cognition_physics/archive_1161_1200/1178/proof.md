# 1178证明：由共同钟响应选出的耦合及同一自主来源回读

2026-10-13。条件性有限模型定理；新增公理及经验校准均为0。响应类、角色合同、量子实施与实际仪器权限均以下文为准。

## 1. 选择定理复用旧工具，固定条件通道不等于非线性背景映射

固定已校准关系钟 $h_t=-\Omega X_t$、原物质公共钟率 $N_s=1+\epsilon z$。在事前声明的线性Hermitian响应类

$$
h_t(z)=h_t+zB_t,\qquad B_t=B_t^\dagger\in M_2(\mathbb C)
\tag{1}
$$

中，候选共同钟角色要求：条件z在非空开区间内，时间在原点邻域内，完整qubit条件通道满足

$$
\operatorname{Ad}_{e^{-it(h_t+zB_t)}}
=\operatorname{Ad}_{e^{-it(1+\epsilon z)h_t}}.
\tag{2}
$$

复用[1170的生成元与可达支持](../1170/proof.md)、[1171的中央核](../1171/proof.md)，在t=0微分给 $[B_t-\epsilon h_t,A]=0$ 对所有 $A\in M_2$。故式(2)当且仅当

$$
B_t=-\epsilon\Omega X_t+\epsilon E_0I,\qquad E_0\in\mathbb R.
\tag{3}
$$

非中央系数 $\lambda=-\epsilon\Omega$ 被相对钟率选定，Y/Z响应被声明的全矩阵合同排除。中央项只改变条件通道的相位，不能由其选择。若仅限定 $B_t=\lambda X_t+\xi I$，也可用已知非零Z钟振幅准备、开放z域与短时Z曲率得到同一系数；那种窄菜单不自动选择一般Y/Z项。

这是额外的**理想数学响应合同**，不是从1177两个有限概率无误差估计出唯一实数。可达域外、隐藏因子、非线性／多体／空间响应仍自由。条件z固定时的通道是线性的；在下面自主背景中改变已知准备q又改变z历史，这份自洽映射不得称为未知Gt上的同一线性通道。全部复制背景都是已知准备，不复制未知Gt、不领取其新增实际全通道权限；实际统一保护的是原未知探针及任意R。

## 2. 同一完整父Hamiltonian实现互易来源

保留[1176完整父H](../1176/proof.md)、CAR源每站4模、探针每站12模及原48准备码、双关系自旋、全部空间边和原探针反作用。新增

$$
\begin{aligned}
H'_{BG,K}&=H_{BG,K}+H_{\mathrm{rec},K},\\
H_{\mathrm{rec},K}&=\frac1K\sum_{r<s,i}
\left(Z_{si,r}B_{ti,s}+B_{ti,r}Z_{si,s}\right),\\
H'_{\mathrm{act}}&=C'_KI+H'_{BG,K}+H_{P,0}
+\sum_{i,c=s,t}F^P_{ci}\overline Z_{ci}.
\end{aligned}
\tag{4}
$$

只作用原同站，复制是内部因子；不含 $r=s$。这份量子组合规则是明示新增输入，不由单个形式通道赠送实际器件。源CAR真空仍精确不变，旧CAR背景W在其上严格消失，但新互易项不消失。所有关系自能、探针相互作用及细记录均保留。

若 $b_*=\|B_t\|$，则 $\|\Delta W\|\le2Vb_*$、$\|H_{\mathrm{rec},K}\|\le V(K-1)b_*$。取

$$
C'_K=C_K+V(K-1)b_*
\tag{5}
$$

即可沿用原完整Fock正谱保证。比较器始终保同一 $C'_K$ 相位。式(4)有限、自伴、自治，总能源精确守恒；不给它的控制、准备或末读免费资源。

对同站四维cluster $G_s\otimes G_t$，Hartree有效作用分离为两单体场。新增每站交互能 $z_s\langle B_t\rangle$ 的两个变分同时给

$$
h_s^{\mathrm{eff}}=-\Omega X_s+bZ_s,\quad
h_t^{\mathrm{eff}}=-\Omega X_t+z_sB_t,\quad
b=\langle B_t\rangle=-\epsilon\Omega x_t+\epsilon E_0.
\tag{6}
$$

因此在式(3)被选定后，非中央回写来源 $-\epsilon\Omega x_t$ 不能再独立填写。正谱平移没有消除 $E_0$：在条件Gt通道它是中央项，内化后却乘 $Z_s$，会改关系后态。复用[997 §3—4](../../archive_990_1008/research_note_997.md)及1171的区分，不把本有限基线称为已定Λ。

## 3. 两份有能源回写、又有钟振幅的背景准备

固定原参数 $(\epsilon,\nu,g,\Omega)=(.15,.35,.2,.5)$、$\ell=16/65$，不优化；$K=10^{14}$、$T=\pi/2$ 沿1177固定。取公开见证分支 $E_0=0$，于是 $\lambda=-.075$。每站、每份复制准备

$$
\mathbf g_s(0)=(1,0,0),\qquad
\mathbf g_t(0)=(q,0,\sqrt3/2),\qquad q=\pm1/2.
\tag{7}
$$

两份准备有不同关系能源，不称同能谱；制备能源与后态都在合同中。Gt的 $x_t=q$ 在Hartree中守恒，且非零Z振幅让它仍可作钟。令 $b=\lambda q$、$r_b=\sqrt{\Omega^2+b^2}$，常场进动给

$$
z_s(t)=-\frac{\Omega b}{r_b^2}[1-\cos(2r_bt)],\qquad
\tau_b(t)=t-\frac{\epsilon\Omega b}{r_b^2}
\left[t-\frac{\sin(2r_bt)}{2r_b}\right].
\tag{8}
$$

因为选定式(3)，关系钟的条件生成元正是 $-\Omega[1+\epsilon z_s(t)]X_t$，故

$$
z_t(t)=\frac{\sqrt3}{2}\cos[2\Omega\tau_b(t)].
\tag{9}
$$

原探针具体见证仍为p=0、$\beta=+1$、质量 $.3,.5$ 的 $|+\rangle$，原 $P_{-Y}$ 概率为

$$
p_D(t)=\frac{1+\sin[.2\tau_b(t)]}{2}.
\tag{10}
$$

该钟感到的N来自同一Gt能源回写，不是给每份准备独立设置外部曲线。p=0空间分布仍均匀，原全部站点细记录保留；本轮不另宣称位置边缘产生传播差或获得始终共位的世界线。

两支的b反号、$r_b$相同，$\tau_+>\tau_-$。其差

$$
\Delta\tau(T)=-\frac{2\epsilon\Omega\lambda(1/2)}
{\Omega^2+\lambda^2/4}
\left[T-\frac{\sin(2\sqrt{\Omega^2+\lambda^2/4}\,T)}
{2\sqrt{\Omega^2+\lambda^2/4}}\right]>0
\tag{11}
$$

使同一原钟的概率差严格正。平均每站完整关系能 $-\Omega x_s-\Omega q+\lambda z_sq$ 守恒。裸Gt能源 $-\Omega q$ 本身守恒；变化的是进动及交互能源分配，未声称抽取或消耗该裸能源。整个实际父过程的总能源守恒由式(4)保证，Hartree能源恒等式不替代它。

## 4. 必要C对照：共同钟修复与非零回写分开检验

另取相同新H、$G_s=+Z,G_t=+Z,E_0=0$。此时 $q=0,b=0$，故Hartree中

$$
z_s(t)=\cos t,\quad \tau_C(t)=t+.15\sin t,\quad
p_{G,C}(t)=\frac{1+\cos\tau_C(t)}2.
\tag{12}
$$

原物质钟仍用式(10)。C与1177旧父违例采用相同已知准备及原末读；区别是新互易H。它验证关系钟按照同一物质钟参数运行，而不是仍按t。C因 $x_t=0$ 没有非中央平均回写，不能单凭C声称完整反馈；式(7)—(11)专门承担这一项。

同站关系读口沿1177：每站第1复制Gt当地CNOT到空白指针，原站点／钟与全部指针细记录保留，后续因果汇集按探针末端站点选相应结果。这是末端齐次同站合同，不是预设单条共位轨迹。准备、正窗、同步、触发、后态及资源的实际误差每份不超过 $\eta_{\mathrm{inst}}=2\times10^{-5}$，是未制造的有限实施合同。

## 5. 有限K：精确中心化、同站相关、完整未知P/R

真空扇区使**不同空间站的精确BG**独立演化，但同站s/t和复制间通常相关。每站四维cluster用[既有计数证明](../_admission/after1174_direction/mean_field_bound.md)，$w_{\mathrm{site}}\le2\|B_t\|=.15$。设Hartree单cluster为 $\phi(t)$、精确BG为 $\Phi_K(t)$，则

$$
\alpha(t)\le\frac{e^{12w_{\mathrm{site}}t}-1}{K},\qquad
\|(\overline Z_c-z_c)\Phi_{K,i}(t)\|
\le u(t)=\frac{4e^{6w_{\mathrm{site}}t}}{\sqrt K}.
\tag{13}
$$

不能说围绕Hartree均值的不同站涨落正交。令精确均值 $\mu_c=\langle\overline Z_c\rangle$，按

$$
\overline Z_{ci}-z_c=(\overline Z_{ci}-\mu_c)+(\mu_c-z_c)
\tag{14}
$$

拆分。零均值项不同站正交；同站s/t用两项Cauchy，保留其相关。偏置空间均匀，不对各站范数求和。复用1177一粒子局部支撑和全局来源界

$$
a_s=4.340625,\quad a_t=1.625,\quad
B^2=7a_s^2+5a_t^2,\quad
\left\|\sum_iF^P_{ci}\right\|\le a_c.
\tag{15}
$$

其中最后两式对原完整一粒子Hilbert空间成立，非仅所选p=0。残差的中心部分至多 $\sqrt2uB$，偏置部分至多 $u(a_s+a_t)$。用真实传播作Duhamel左传播，得到全部原48输入、跨质量相干及任意R的准备等距界

$$
\left\|U'_{\mathrm{act}}(T)E-
\bigl[U'_{BG,K}(T)\Phi_0\bigr]\otimes U_P(T)\right\|
\le\eta_{P,K}
\le\frac{4Te^{6w_{\mathrm{site}}T}}{\sqrt K}
(\sqrt2B+a_s+a_t)<6\times10^{-5}.
\tag{16}
$$

右侧按未知输入作用，带相同 $C'_K$ 相位。比较背景始终是**精确BG整体**，不是全产品近似；真实父过程包含探针反作用。后接同一完整终端记录CP保持该误差。任意新中途干预、非均匀准备或未知跨复制背景不得直接继承此界。

若读取单站单复制关系钟，纯态保真度和计数式(13)另外给其半迹距不超过 $\sqrt\alpha$。因此定义

$$
\delta_D=6\times10^{-5}+2\times10^{-5}=8\times10^{-5},\qquad
\delta_G=\delta_D+4\times10^{-7}=8.04\times10^{-5},
\tag{17}
$$

其中 $\sqrt{\alpha(T)}<4\times10^{-7}$。按实际探针站点汇总的关系效果同样成立：比较背景与探针分离后各站相同，单站误差以位置概率权重相加为1。这个局部保证不等于全部指针联合分布接近Hartree产品；全部指针和BG仍完整保留。

## 6. 有限证书及数值责任

两份回写准备原钟差的严格下界为

$$
p_{D,+}^{\mathrm{act}}-p_{D,-}^{\mathrm{act}}
\ge\frac{\sin(.2\tau_+)-\sin(.2\tau_-)}2-2\delta_D
>0.00106.
\tag{18}
$$

数值由有理Machin区间、整数平方根、三角Taylor余项和正指数级数余项认证，不依赖步差或浮点本征精度。固定粗误差足够，不再收紧或调参数。

对C，复用1177事前单调窗内 $R(p)=[1+\cos(5\arcsin(2p-1))]/2$、Lipschitz常数6。代码以有理端点验证三份任务的整个 $p_D\pm\delta_D$ 区间严格处于 $[1/2,(1+\sin(.38))/2]$ 内，故实际概率仍在同一个 $0\le\tau\le1.9$ 单调域，而非只检查理想中心。理想C残差为0，所以真正有限过程满足

$$
|p_{G,C}^{\mathrm{act}}-R(p_{D,C}^{\mathrm{act}})|
\le\delta_G+6\delta_D=0.0005604.
\tag{19}
$$

这与1177旧父单支残差下界 $.0745650662$ 分离；不能将修复同率的C误当非零能源反馈证明。两份q准备的关系钟振幅是 $\sqrt3/2$，其对应角色曲线也按同一振幅校准，不错误沿用振幅1。

任何实际采样误差另加；式(18)每个频率若还有 $10^{-4}$ 估计误差，下界再减 $.0002$ 仍正。本轮是有限理论证书，没有经验实验。小矩阵只核生成元、非零回写及记录符号；不声称模拟巨大完整父态。

## 7. 剩余自由与停止

额外输入包括共同钟角色的理想量词、线性Hermitian响应类、同站互易的量子实施及实际仪器合同。被选出的只有相对非中央系数及由同一H导出的互易项。$E_0=0$仅为见证分支；非零基线改变来源及有限误差，不能沿用零分支数据。

尚未建立Gs自身及全部交互作为共同几何参与者的适用身份、完整来源和响应；不能预设全部裸自由流都必须乘同一个N。空间推进、非线性约束和健康模式未被选择。不能把本项机制称为普适自耦或Einstein方程，不提前删除额外模式、暗部门，也不以继续修本候选作为完成整体研究的唯一道路。
