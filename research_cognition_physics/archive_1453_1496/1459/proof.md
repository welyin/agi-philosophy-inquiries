# 从同一凝聚的两道部分宽度到绝对寿命和端点质量

## 归一和消元

取一个正则复弱三重态，电荷取 $Q=T_3+Y$ 且 $Y=1$，中性分量 $\langle\Delta^0\rangle=v_\Delta/\sqrt{2}$。在带电轻子质量基中 $M_\nu=\sqrt 2 v_\Delta h$。规范顶角为 $i\sqrt2g^2v_\Delta g_{\mu\nu}$，可从 $\operatorname{tr}(D_\mu\Delta)^\dagger D^\mu\Delta$ 的两次弱升降作用取得。双电荷态质量 $M>2m_W$，两个 $W$ 取在壳。横向的两个偏振给2，纵向给 $(M^2/(2m_W^2)-1)^2$，再乘相同粒子相空间。实测 $G_F$ 的使用要求标准缪子衰变树级匹配 $g^2/(8m_W^2)=G_F/\sqrt2$，其他四轻子贡献或圈匹配须在实际比较时校准。得到

$$
x=m_W^2/M^2,\quad f(x)=\sqrt{1-4x}(1-4x+12x^2),\qquad
\Gamma_W=\frac{G_F^2 v_\Delta^2 M^3}{2\pi}\, f(x).
$$

令 $B=M_\nu^\dagger M_\nu$，在以下宽度公式中质量用 GeV，则 $\Gamma_\ell=M\operatorname{tr} B/(16\pi v_\Delta^2)$，$r_e=B_{ee}/\operatorname{tr} B$。$\Gamma_\ell$ 只求和六个同号带电轻子道。独立文献已分别推导两种宽度、双 $W$ 的偏振成分和衰变长度，采用其领先相互作用。[Fileviez Pérez 等，式 26–28](https://arxiv.org/abs/0805.3536)

定义 $\kappa=G_F^2 M^4 f(x)/(32\pi^2)$，为无量纲。三个直接的正、负或有界观测组合为

$$
\Gamma_W\Gamma_\ell=\kappa\sum_i m_i^2,\qquad
\Gamma_W\Gamma_\ell r_e=\kappa m_\beta^2,\qquad
\Gamma_W\Gamma_\ell(1-3r_e)=\kappa G,\quad G=\operatorname{tr}B-3B_{ee}.
$$

$G$ 完全由两条振荡质量平方差及电子行混合决定。KATRIN 五批数据给出 $m_\beta<0.45$ eV（90% 置信水平），因此若未来能分别测量双 $W$ 分宽以及电子腿加权分宽，可与独立的氚端点结果核对。单独电子对的宽度涉及另一个矩阵元素；这里必须计入 $e\mu,e\tau$。三态幺正轻混合，狭窄极点、相同低能质量源以及共同弱流都是输入。[KATRIN](https://doi.org/10.1126/science.adq9592) [第1456轮振荡量与符号](../1456/proof.md)

## 把寿命留作独立观测

写 $b_W=\Gamma_W/\Gamma_\mathrm{tot}, b_\ell=\Gamma_\ell/\Gamma_\mathrm{tot}$。对 $G\ne0$，$\Gamma_\ell,\Gamma_W>0$，

$$
\Gamma_\mathrm{tot}^2=\frac{\kappa G}{b_W b_\ell(1-3r_e)},\qquad
c\tau=\hbar c\sqrt{\frac{b_W b_\ell(1-3r_e)}{\kappa G}}.
$$

两个有符号量应同号。可仅用分支比预报寿命，但 $b_W,b_\ell$ 各以**全部**衰变为分母。在只知道两道相对比例时，未知的其余分支会随总宽缩放；这个测量条件独立于相对分味筛选。若 $r_e$ 接近 $1/3$，提取须传播味效率和相关误差，符号穿零时不能给精确寿命。数值例明确加入32%的第三通道，独立分支与振荡量恢复同一总宽。[可复算例](results.json)

无需测任何分支比也有 $\Gamma_\mathrm{tot}\ge\Gamma_W+\Gamma_\ell\ge 2\sqrt{\kappa\,\operatorname{tr}B}$。前一界在无其他通道时取等，后一界可调 $v_\Delta$ 使两道相等。正常序中 $\operatorname{tr}B=3m_{\min}^2+\Delta_{21}+\Delta_{31}$，倒序中 $\operatorname{tr}B=3m_{\min}^2+\Delta_{21}+2|\Delta_{31}|$。采用已登记的矩形取下端后，分别至少为 $0.002531281,0.004936281$ eV$^2$，给报告中各质量点的最长固有长度。上界对额外衰变通道仍成立。这些是树级在壳有限范围内的理论界，未赋予矩形联合置信度。$WW^*$ 或双离壳区须重算相空间。[JUNO](https://arxiv.org/html/2609.37895v1) [输入及单位](inputs.json)

## 能被什么打破

引入另一套独立的轻中微子质量，例如额外重单态的低能贡献。极限 $h,v_\Delta\to0$ 时保持 $M_{\rm other}$ 给出观测振荡，两个部分宽度均可趋零；若没有开放的快速级联，寿命即可远超上述界。即使只将实际 $WW$ 顶角乘实因子 $k$，三式的右边皆乘 $k^2$，长度上界放宽为原来的 $1/|k|$。自由的高能匹配、额外标量混合和强烈味依赖运行需要独立标定。该反例定位共同作用的可证伪前提，同时保留基准标准模型的其他读数。
