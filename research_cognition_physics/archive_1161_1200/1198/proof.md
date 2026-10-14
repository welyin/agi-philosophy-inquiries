# 1198证明：保源反馈的两关系低频响应界

2026-10-14。限定于[1196自洽填带参考](../1196/derivation.md)的Hartree受驱线性响应。原有限K父Hamiltonian不变；不把其Hartree参考升级为精确有限K真空。

## 1. 要判定的命题与输入

沿[1176半边来源](../1176/proof.md)和[1178互易项](../1178/proof.md)，保固定 $\ell=16/65$、$m=2$、$\epsilon=3/20$、$\nu=7/20$、$g=1/5$、$\Omega=1/2$、$E_0=0$。1196的固定点证书给 $z_s>49/50$，$z_t=0$，$x_t=1$，并且

$$
N=1+\epsilon z_s,\quad v=1+\nu z_s,\quad
h_p=mN\beta+v\boldsymbol\alpha\cdot k_p,\quad
E_p=\sqrt{m^2N^2+v^2|k_p|^2},\quad P_\pm(p)=\tfrac12(I\pm h_p/E_p).
$$

源占据为 $P_-(p)$，不丢任何原CAR模式。关系准备为各自有效场的低能纯态。数学微扰是同一Hartree方程在该参考的受驱线性化；它不授予一个新实际驱动器、时钟或测量权限。

**结论：** 对原盒内所有 $q$ 和实频率 $|\omega|\le1/4$，在声明的关系权重范数中，两关系与完整源虚跃迁的反馈算符范数严格小于 $.835$，故反馈逆有大于 $.165$ 的余量。这些原关系方向在该频带没有受驱集体极点。

界不含最小非零动量。下面进一步证明：保持相同格距和耦合、每轴格点数至少3的任意周期盒，其同类唯一各向同性填带固定点都满足 $z_s>.98$。无需逐盒搬用原65³固定点数值。

### 固定点下界对所有这些盒一致

对 $n_b\ge3$，离散几何级数给 $\langle\sin^2(2\pi j/n_b)\rangle=1/2$，故 $\langle|k|^2\rangle=3/(2\ell^2)$。在 $z_0=49/50$、$N_0=1+\epsilon z_0$、$v_0=1+\nu z_0$ 处，

$$
E_{\max}^2=(mN_0)^2+3(v_0/\ell)^2
=\frac{24208434691}{256000000}<100,
$$

$$
|b(z_0)|=2\left\langle\frac{\epsilon m^2N_0+\nu v_0|k|^2}{E_p}\right\rangle+\epsilon\Omega
\ge B_0:=\frac{2[\epsilon m^2N_0+\nu v_0\,3/(2\ell^2)]}{10}+\epsilon\Omega
=\frac{130044843}{51200000}.
$$

精确算术给

$$
B_0^2(1-z_0^2)-(\Omega z_0)^2
=\frac{100735097898570251}{6553600000000000000}>0.
$$

因此原固定点映射 $F(z_0)=|b(z_0)|/\sqrt{\Omega^2+b(z_0)^2}>z_0$。1196全区间 $[-1,1]$ 的压缩界直接继承：其推导只有逐动量上界及归一平均，不依赖盒格点数。$F$连续、导数小于1且 $F(1)<1$，故每个此类盒的唯一固定点都大于 $z_0$。这是有限盒族的一致定理；不宣称另有已完成的无限体积量子真空构造。

## 2. 来源顶角由原边定义确定

每站 $F_{si}=\epsilon m_i+\nu\sum_b t_{bi}$，$F_{ti}=g(t_{xi}-t_{yi})$。每条跃迁边在两端各分一半，因此位置算符形式的变分是

$$
D_s(q)=\epsilon Z_q M+\frac\nu2\{Z_q,T_x+T_y+T_z\},\qquad
D_t(q)=\frac g2\{Z_q,T_x-T_y\},\quad (Z_q)_{ii}=e^{iq\cdot x_i}.
$$

$M$是质量单粒子矩阵，$T_b$是原轴向跃迁。令 $p'=p+q$，得到

$$
D_s(p,q)=\epsilon m\beta+\frac\nu2\boldsymbol\alpha\cdot(k_p+k_{p'}),\qquad
D_t(p,q)=\frac g2\{\alpha_x(k_{px}+k_{p'x})-\alpha_y(k_{py}+k_{p'y})\}.
\tag{1}
$$

这是真实来源的变分，不是事后给定的应力权重。$D(q)^\dagger=D(-q)$；单个非零q算符一般非Hermitian，不能只保一个跃迁方向。

## 3. 两组跃迁Gram给正半定响应

由 $i\dot C=[h,C]$，在能量本征基中

$$
(\delta C_b)_{mn}=\frac{f_n-f_m}{\omega-(e_m-e_n)}(D_b)_{mn},\qquad
\Pi_{ab}=V^{-1}\operatorname{tr}(D_a^\dagger\delta C_b).
\tag{2}
$$

两条能带的占据固定为0和1；带内项的分子为零。对每p定义 $\Delta_p=E_p+E_{p'}$，以及

$$
M^{\sigma'\sigma}_{ab}(p,q)
=\operatorname{tr}[D_a^\dagger P_{\sigma'}(p')D_bP_\sigma(p)].
$$

每个M都是Hilbert–Schmidt向量 $P_{\sigma'}D_bP_\sigma$ 的Gram矩阵，正半定。因而实频率 $|\omega|<\Delta_0:=2mN$ 时

$$
Q(q,\omega):=-\Pi(q,\omega)
=\frac1V\sum_p\left[
\frac{M^{+,-}(p,q)}{\Delta_p-\omega}
+\frac{M^{-,+}(p,q)}{\Delta_p+\omega}\right]\succeq0,
$$

$$
0\preceq Q(q,\omega)\preceq
\frac{Q(q,0)}{1-|\omega|/\Delta_0}.
\tag{3}
$$

这是含全部原动量和两通道交叉关联的响应；没有令 $Q_{st}=0$。它不是对所有初态改变的演化算符：初始占据、混合度或其它守恒标签的改变可产生另外的非齐次来源，不由式(2)消除。

## 4. 全q的静态上界

秩2投影给 $M_{aa}^{\sigma'\sigma}\le2\|D_a\|^2$。Clifford反对易及式(1)给

$$
\|D_s\|^2\le\epsilon^2m^2+\frac{\nu^2}2(|k_p|^2+|k_{p'}|^2),\qquad
\|D_t\|^2\le\frac{g^2}2(|k_p|^2+|k_{p'}|^2).
$$

记 $k_{\max}=\sqrt3/\ell$。因为 $E_p\ge mN$，且 $|k_p|^2\le k_{\max}E_p/v$，逐项求和可得

$$
Q_{ss}(q,0)\le B_s:=\frac{2\epsilon^2m}{N}+\frac{2\nu^2\sqrt3}{v\ell},\qquad
Q_{tt}(q,0)\le B_t:=\frac{2g^2\sqrt3}{v\ell}.
\tag{4}
$$

这一步不依赖小q展开、格点数或源与关系解耦的假设。

## 5. 关系易感率与互易项

原有效关系场为 $h_s=-\Omega X_s+bZ_s$、$h_t=-aX_t$，其中 $r=\sqrt{\Omega^2+b^2}$、$a=\Omega N$。直接线性化Bloch方程，Z场到Z均值的响应为

$$
\chi_s(\omega)=-\frac{\Omega^2/r^3}{1-\omega^2/(4r^2)},\qquad
\chi_t(\omega)=-\frac{1/a}{1-\omega^2/(4a^2)}.
\tag{5}
$$

互易项 $-\epsilon\Omega z_s x_t$ 给 $\delta b_s=\delta\langle F_s\rangle-\epsilon\Omega\delta x_t$、$\delta h_{t,x}=-\epsilon\Omega\delta z_s$。在纯 $+X_t$ 的受驱切向中 $\delta x_t=0$，沿X的场变化也与原参考对易，所以没有额外一阶交叉项。互易项已经进入a和r，未被从父H中删去。任意初始混合度变化若使 $\delta x_t\ne0$，则是另外的初始来源，不能归入本受驱核的无极点结论。

## 6. 完整反馈逆与严格余量

令 $C(\omega)=-\operatorname{diag}(\chi_s,\chi_t)>0$，$\eta$是关系的外加数学测试源。原源变化 $\delta f=\Pi\delta z=-Q\delta z$，因此

$$
(I-CQ)\delta z=-C\eta.
\tag{6}
$$

$CQ$相似于 $C^{1/2}QC^{1/2}$。由正半定性，其算子范数不大于迹；式(3)—(5)给对 $|\omega|\le w=1/4$ 一致的上界

$$
\|C^{1/2}QC^{1/2}\|
\le\frac1{1-w/(2mN)}
\left[
\frac{(\Omega^2/r^3)B_s}{1-w^2/(4r^2)}
+\frac{B_t/(\Omega N)}{1-w^2/(4\Omega^2N^2)}
\right].
\tag{7}
$$

保守取 $N\ge1147/1000$、$v\ge1343/1000$、$\sqrt3<17321/10000$、$\sqrt{1-z_s^2}<199/1000$。代码以精确分数给

$$
\text{式(7)}\le
\frac{373901331379513018111433}{448151477922597240000000}
<\frac{835}{1000},
$$

$$
\lambda_{\min}(I-C^{1/2}QC^{1/2})
\ge\frac{74250146543084221888567}{448151477922597240000000}
>\frac{165}{1000}.
\tag{8}
$$

分母在本频带非零，各源谱项与裸关系响应有限。因此式(6)的受驱响应没有极点。这里的余量属于加权反馈矩阵，不是全H能隙，也不是任意原始坐标下的最小奇异值。

原操作比较器的 $\delta N=\epsilon\delta z_s$、$\delta v_x=\nu\delta z_s+g\delta z_t$、$\delta v_y=\nu\delta z_s-g\delta z_t$、$\delta v_z=\nu\delta z_s$ 是固定线性读出；该读出不能从有限响应创造新极点。尚未给任意有限K实际钟概率的全频认证。

## 7. 判别力和停止线

本参考的同一受驱关系响应不能直接充当趋近零频的无质量几何传播部门。全q界使这不是仅由有限盒排除零模式的观察，提供停止继续精调此参考boost谱的独立理由。

但本结论不排除其它准备、不同关系规则、有限频带的受控近似、初态来源、有限K新集体效应或其它量子引力路线。尤其没有用无限长程要求否定所有有限有效GR任务。若未来要声明某一有限GR误差任务失败，仍须对同一实际源/钟/端点及预算比较预测，不能用式(8)替代。

暗物质、暗能量身份与整体认知公理的蕴涵命题均未由此判定。旧1178的有限时间反馈和1197的完整量子二阶差仍成立，无极点与非零响应不矛盾。

## 8. 核验

[check.py](check.py)的第一组以精确有理数复算式(8)，读取冻结1196根区间并绑定哈希。第二组独立从3³实空间108阶单粒子矩阵构造半边来源，再与4阶动量块比较顶角和源线性响应；后者只是数值代数核查，不替代全q解析证明或原65³准备。

独立审阅及结果见同目录正式记录；当前是否验收以acceptance.json为准。
