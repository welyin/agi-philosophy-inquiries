# 1281直接接续：多方向比较与实际准备域

2026-10-17。非编号直接推论；独立推导`/root/curved_contact_scope`，主线核对。没有新增科学轮次或公理，不把1016已知的反对称标签冗余重新记为发现。

## 1. 完整正比较算符的交叉项

在 $(n+1)d$ 维直和上，设H左上块A、交叉块 $-iC_a$、下部块K_ab，且整个H≥0。给定各R_a Hermitian，令Y_a上槽为iR_a，第a下槽为I，其他为0。每方向独立满秩零校准使 $HY_a=0$，故

$$
C_a=AR_a,\qquad K_{ab}=R_aAR_b.
$$

有限误差直接继承1281。若 $Y_a^\dagger HY_a\le\tau_aI$、$\delta_a=\sqrt{\|A\|\tau_a}$，则

$$
\|K_{ab}-R_aAR_b\|
\le\sqrt{\tau_a\tau_b}+\|R_a\|\delta_b+\|R_b\|\delta_a.
$$

证明用H的Cauchy–Schwarz控制 $Y_a^\dagger HY_b$，与1281同一展开即可。没有一项要求R_a的线性张成已对交换子闭合；因此不能据此宣称选出Lie代数或规范群。

## 2. 只在实组合上非负不能升级为整体正性

取d=1、n=2、A=R_1=R_2=1：

$$
H_t=\begin{pmatrix}1&-i&-i\\i&1&1+it\\i&1-it&1\end{pmatrix}.
$$

对所有 $v,\phi\in\mathbb C$、$b_1,b_2\in\mathbb R$，在 $w_a=b_a\phi$ 上其二次型恒为 $|v-i(b_1+b_2)\phi|^2$。每个单方向校准为零。但归一相干向量

$$
\psi_\pm=\tfrac12(i(1\pm i),1,\pm i)^T
$$

给 $\psi_\pm^\dagger H_t\psi_\pm=\mp t/2$，所以t≠0时H_t并非整体正。不能仅在受限实组合上检查正性就套用1281。允许方向间相对相位可探测这个差别，但必须实际登记该菜单。

可直接复算：

```python
import numpy as np
t = 1/3
H = np.array([[1,-1j,-1j],[1j,1,1+1j*t],[1j,1-1j*t,1]])
for sign in (1,-1):
    p = np.array([1j*(1+sign*1j),1,sign*1j])/2
    assert np.isclose(np.vdot(p,H@p).real,-sign*t/2)
```

## 3. 与场论接口的区别

1016已经说明反对称规范标签与对称实字段乘积收缩恒零；不能将上述扩大后的相干槽自动认作原场论已有观测。完整PSD与受限字段正性是不同模型范围。

所以此处不另开1282去重复块推广。下一真正有用的问题是：能否在**旧的同一物质载体**中，凭独立校准选择原来单独输入的比较作用，并直接继承其来源和读口证书。557的实际规范链路是待查重候选，尚未准入。
