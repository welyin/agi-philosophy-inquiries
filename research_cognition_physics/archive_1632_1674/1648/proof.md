# 区分两条质量路径的 Higgs 径向响应

## 物理输入

第1489轮在正对角上夸克质量基内构造下夸克矩阵。两条路径到达同一端口。采用原轮的质量单位，

$$
Y=\begin{pmatrix}0&A&0\\B&C&Z\\0&E&F\end{pmatrix},\qquad Z=R+D e^{i\eta},\qquad H=YY^\dagger .
$$

七个具名系数均可取有号实数，径向小扰动保持相位 $\eta$ 及纹理零点。取正则归一、沿电弱真空 $v$ 增大的单一物理 Higgs 模，其符号可由 $W$ 质量响应定标。将 $v\,\partial_v M$ 用同一单位记作 $Y+\delta Y$。每条边的额外对数斜率为实数 $n_e$，即 $\delta Y_e=n_eY_e$，故 $\delta Y_{23}=X=n_RR+n_DD e^{i\eta}$。一条路径若比线性质量多含一份规范不变因子 $H_{\rm EW}^\dagger H_{\rm EW}$，则有 $n_e=2$。这些物理建模输入只用于树级匹配；应用至实验还要核算真实场混合、尺度运行和极点顶角。[原矩阵](../../archive_1453_1496/1489/proof.md) · [六维算符](https://arxiv.org/abs/1008.4884)

在 $U_L^\dagger YU_R={\rm diag}(y_i)>0$ 的质量基中定义

$$
p_i=\operatorname{Im}[U_L^\dagger(Y+\delta Y)U_R]_{ii}/y_i,\qquad\lambda_i=y_i^2 .
$$

$p_i$ 是有号对角赝标量耦合除以该味的线性质量响应。质量固定正实后，剩余场换相不影响它。单个无极化部分宽度一般只给标量和赝标量的模方和；有号读数需要干涉、自旋或其它 CP 灵敏观测。[赝标量 Yukawa 的物理意义](https://arxiv.org/abs/1310.1385)

## 两个可反演的矩

$\sum_i p_i=\operatorname{Im\,tr}(Y^{-1}\delta Y)=0$，因为 $\det Y=-ABF$ 的径向对数导数为实。再令 $I=\operatorname{Im}(H_{12}H_{23}H_{31})=A^2CEF\operatorname{Im}Z$，$T=H_{22}+F^2$，则

$$
\begin{aligned}
\sum_i\lambda_i p_i &= R\operatorname{Im}Z\;(n_D-n_R),\\
\sum_i\lambda_i(\lambda_i-T)p_i &=\frac{I}{H_{11}}\;(n_D+n_E-n_C-n_F).
\end{aligned}
$$

第一矩是 $\operatorname{Im\,tr}(\delta Y Y^\dagger)$，只有两条干涉路径的交叉项留下。计算 $\operatorname{Im\,tr}(H\delta Y Y^\dagger)$ 时，第23元贡献 $T\operatorname{Im}(XZ^*)+CEF\operatorname{Im}X$。实边 $C,E,F$ 分别贡献 $(-n_C,+n_E,-n_F)CEF\operatorname{Im}Z$，另外两边对虚部为零。扣除 $T$ 倍第一矩即得第二式。三个 $\lambda_i$ 两两不同时，三行矩组成可逆的 Vandermonde 系统，故两个斜率差能预报全部三个 $p_i$。若 $R\operatorname{Im}Z$ 与 $I$ 均非零，也可从有号 $p_i$ 反推两个差。$I$ 可由三代弱流不变量及质量平方差定标。把 $Z$ 分解为两条路径仍需路径相位的身份或独立校准。

这两个组合穷尽本读口能分辨的自由。三左三右端口由七条边连通。若每条边的斜率增加 $l_i+r_j$，则 $\delta Y$ 增加 $L Y+Y R_{\rm port}$，其中 $L,R_{\rm port}$ 实对角。换到质量基，新增的对角顶角均为实。此类加法有五个独立方向，两个已测组合又确实独立，因此该五维空间就是一般点的完整核。它们对对角赝标量不可见，仍可改变 CP 偶顶角或换味读数。质量简并、弱环消失或某条路径为零时，反向识别须重算。[直接导数及奇异值分解](check.py)

## 相同质量下的标量预报

在第1489轮已用夸克谱及弱相位定标的中央矩阵中取 $\eta=-120^\circ$。只改变各路径的径向依赖，保持其在真空处的值。最后一列用无夸克质量近似的树级部分宽度，两个末态合计给 $(|g_{bs}|^2+|g_{sb}|^2)/|g_{bb}|^2$。

| 额外斜率 | $p_d$ | $p_s$ | $p_b$ | $\Gamma(b\bar s+\bar b s)/\Gamma(b\bar b)$ |
|---|---:|---:|---:|---:|
| 仅相位路径 $n_D=2$ | +1.42137 | −1.42228 | +0.000904 | 0.00108948 |
| 仅实路径 $n_R=2$ | +0.48416 | −0.48393 | −0.000234 | 0.00023869 |
| 两路径同为2 | +1.90553 | −1.90620 | +0.000670 | 0.00183536 |
| 所有边线性 | 0 | 0 | 0 | 0 |

复共轭全部振幅时，有号 $p_i$ 全部反向，两体宽度比保持。最后两行的第一矩同为零，可由第二矩区分。数表描述指定的局部多项式质量模型；平方振幅保留了高次项，不能充作普适的线性维六截断。第1534轮已分析行列式相位导数和不同味的抵消；本轮新增的是三零纹理中的两个可反演矩和五维核。现实的路径斜率仍待检验。弱相位用于定标，新预报在标量部门。[原行列式结果](../../archive_1497_1544/1534/proof.md) · [当前结果](results.json) · [独审](independent_review.md)
