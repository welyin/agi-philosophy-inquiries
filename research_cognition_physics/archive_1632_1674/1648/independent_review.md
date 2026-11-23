# 第1648轮独审　三零质量纹理中可由对角径向赝标量顶角分辨的路径差

沿用第1489轮在一个共同匹配基中保存的三零矩阵及两条相干路径。该矩阵的 CKM 相位已经校准。本轮检验的是在相同质量与相位下，改变实径向斜率会怎样改变一个实际标量的夸克顶角。物理 Higgs 必须与所选径向方向按登记方式匹配；不同方向混合、动量依赖和可探测率的额外传递仍待确定。[原质量重建](../../archive_1453_1496/1489/results.json)

设 $H=YY^\dagger$，奇异质量平方 $\lambda_i=y_i^2>0$ 互异，左右奇异列为 $L_i,R_i$，$\delta Y_{23}=X$。定义 $V=\delta Y Y^\dagger$。以 $R_i=Y^\dagger L_i/y_i$ 消去右列，有 $p_i=\operatorname{Im}(L_i^\dagger V L_i)/\lambda_i$。三项对角矩于是成为弱基中的迹：

$$
\sum_i p_i=\operatorname{Im}\operatorname{tr}(Y^{-1}\delta Y)=0,\qquad
A_1=\sum_i\lambda_i p_i=\operatorname{Im}(XZ^*).
$$

第一式来自行列式 $-ABF$ 与 $Z$ 无关。只改 $23$ 时，$V_{22}=XZ^*$、$V_{23}=XF$，其余全零。$H_{32}=CE+FZ^*$，所以 $\operatorname{Im}\operatorname{tr}(HV)=(H_{22}+F^2)A_1+CEF\operatorname{Im}X$。写 $Z=R+D e^{i\eta}$、$X=q_RR+q_DD e^{i\eta}$，则

$$
A_1=(q_D-q_R)R\operatorname{Im}Z,\qquad
I=\operatorname{Im}(H_{12}H_{23}H_{31})/H_{11}=CEF\operatorname{Im}Z .
$$

从第三矩减去 $(H_{22}+F^2)A_1$ 得 $q_D I$。消去另两味的 Vandermonde 方程可逐味解得

$$
p_i=\frac{q_D I+(\lambda_i-A^2-E^2)A_1}{(\lambda_i-\lambda_j)(\lambda_i-\lambda_k)},\qquad \{i,j,k\}=\{1,2,3\}.
$$

上述 $A,C,E,F$ 指原纹理中的实质量元；$A_1$ 是一阶矩。$I$ 和 $R\operatorname{Im}Z$ 同时非零时，三项 $p_i$ 可以反解 $q_D$ 与 $q_R$。$q_D=q_R$ 只使 $A_1$ 消失，束缚的第二个回路仍可给各味相反的赝标量耦合。若径向插入同时给整张矩阵共同缩放，则所有 $p_i$ 均为零。

| 第1489轮中央矩阵的额外斜率 $(q_R,q_D)$ | $p_d$ | $p_s$ | $p_b$ |
|---:|---:|---:|---:|
| $(0,2)$ | 1.421374 | −1.422277 | 0.000903712 |
| $(2,0)$ | 0.484159 | −0.483926 | −0.000233617 |
| $(2,2)$ | 1.905533 | −1.906203 | 0.000670094 |
| $(0,0)$ | 0 | 0 | 0 |

数值是相对于各味自身标量质量耦合的无量纲系数。直接奇异值分解、矩迹和逐味有理表达式相符。36组非简并质量及各自共轭另作复算，原相位共轭会翻转全部 $p_i$。[独立程序](independent_check.py) · [完整结果](independent_results.json)

七条质量路径均可有独立实斜率。按 $A,B,C,E,F,R,D$ 排序，直接展开 $\operatorname{Im}\operatorname{tr}(HV)$ 的 $23$ 与 $32$ 项，第二个矩的系数变成 $n_D+n_E-n_C-n_F$；第一个矩的系数为 $n_D-n_R$。这两个整系数行向量秩为2。把两组端口的六个数相加到连接它们的每条边，得到的 $7\times6$ 矩阵秩5，其每一列都被这两行消掉。因此在上述非退化条件下，两矩零空间恰好由端口斜率生成。独立程序对七条单位边、25组一般实斜率及一组端口斜率验证了全过程。端口例的三个赝标量读数消失，同时仍改变 CP 偶标量顶角。这个核只描述这组三项对角观测。
