# 第5章（续）声学方程的有限体积离散与时间推进

声学计算要求同时控制幅值误差、相位误差和边界反射。与低速流动求解相比，声波传播对时间项、对流项和人工边界更加敏感：局部残差很小并不必然意味着远场声压准确，长时间推进中的微小数值色散也会逐渐积累为明显的相位差。

本章以课件采用的格点型有限体积法为主线，离散含平均流对流效应的声学波动方程。未知声压存储在网格节点，介质密度、声速和马赫数存储在原始单元内；围绕节点构造对偶控制体，通过形函数计算单元内声压梯度。内容依次包括空间离散、边界条件、显式与 Newmark 时间积分、稳定性、数值色散和程序验证。节号与公式编号承接本章前三篇。

本章用 $p'$ 表示声学扰动压力，用 $\mathbf M=\mathbf u_0/c_0$ 表示背景马赫数矢量；为避免混淆，离散质量、阻尼和刚度矩阵分别记为 $\mathsf M$、$\mathsf C$ 和 $\mathsf K$。

## 5.20 声学控制方程与离散变量

### 5.20.1 含平均流的声学波动方程

在小振幅扰动、平均流场已知且熵波和涡波影响可按模型假设处理的条件下，课件采用统一形式的声学波动方程：

$$
\frac1{\rho_0c_0^2}\frac{\partial^2p'}{\partial t^2}
=\nabla\cdot\left(\frac1{\rho_0}\nabla p'\right)
-\frac2{\rho_0c_0}\frac{\partial(\mathbf M\cdot\nabla p')}{\partial t}
-\frac1{\rho_0}\mathbf M\cdot\nabla(\mathbf M\cdot\nabla p').
\tag{5.227}
$$

其中 $\rho_0$ 和 $c_0$ 为背景密度与声速，$\mathbf u_0$ 为平均流速度，且

$$
\mathbf M=\frac{\mathbf u_0}{c_0}
=(M_x,M_y,M_z).
\tag{5.228}
$$

右端第一项表示静止介质中的扩散型声学算子，第二项是混合空间—时间导数，第三项对应平均流引起的二阶对流修正。静止介质中 $\mathbf M=0$，方程退化为经典波动方程。

下面从线性化 Euler 方程推导式（5.227）。推导的关键是先消去扰动速度，得到沿平均流物质导数表示的波动方程，再展开二阶物质导数。严格的消元和求导交换针对局部均匀、定常平均流；空间变系数背景的处理见本步末的适用性说明。

**第一步：变量分解与基本假设。** 将瞬时变量分解为平均量与声学扰动量：

$$
p=p_0+p',\qquad
\rho=\rho_0+\rho',\qquad
\mathbf u=\mathbf u_0+\mathbf u'.
$$

假设流体无黏，忽略黏热耗散和外力；声学扰动等熵，不单独求解熵波；扰动振幅足够小，可忽略二阶及以上扰动项。等熵关系给出

$$
p'=c_0^2\rho',\qquad
c_0^2=\left(\frac{\partial p}{\partial\rho}\right)_s.
$$

定义沿平均流的物质导数算子

$$
\mathcal D_0=\frac{\partial}{\partial t}
+\mathbf u_0\cdot\nabla,
\qquad
\mathcal D_0^2p'=\mathcal D_0(\mathcal D_0p').
$$

**第二步：线性化连续方程。** 无黏可压缩流的连续方程为

$$
\frac{\partial\rho}{\partial t}
+\nabla\cdot(\rho\mathbf u)=0.
$$

代入变量分解，扣除平均场连续方程，并忽略二阶项 $\nabla\cdot(\rho'\mathbf u')$，得到

$$
\frac{\partial\rho'}{\partial t}
+\nabla\cdot(\rho_0\mathbf u'+\rho'\mathbf u_0)=0.
$$

对于局部均匀、定常平均流，展开上式可得

$$
\frac{\partial\rho'}{\partial t}
+\mathbf u_0\cdot\nabla\rho'
+\rho_0\nabla\cdot\mathbf u'=0,
$$

即

$$
\mathcal D_0\rho'+\rho_0\nabla\cdot\mathbf u'=0.
$$

利用 $\rho'=p'/c_0^2$，再除以 $\rho_0$，得到

$$
\frac1{\rho_0c_0^2}\mathcal D_0p'
+\nabla\cdot\mathbf u'=0.
$$

**第三步：线性化动量方程。** 忽略黏性力后的动量方程为

$$
\rho\left(
\frac{\partial\mathbf u}{\partial t}
+(\mathbf u\cdot\nabla)\mathbf u
\right)=-\nabla p.
$$

代入变量分解，扣除平均场方程并忽略二阶扰动项。对于均匀平均流，$\nabla\mathbf u_0=0$，于是

$$
\rho_0\left(
\frac{\partial\mathbf u'}{\partial t}
+(\mathbf u_0\cdot\nabla)\mathbf u'
\right)=-\nabla p',
$$

即

$$
\rho_0\mathcal D_0\mathbf u'=-\nabla p',
\qquad
\mathcal D_0\mathbf u'=-\frac1{\rho_0}\nabla p'.
$$

**第四步：消去扰动速度。** 对线性化连续方程的压力形式作用一次 $\mathcal D_0$：

$$
\frac1{\rho_0c_0^2}\mathcal D_0^2p'
+\mathcal D_0(\nabla\cdot\mathbf u')=0.
$$

在局部常系数、均匀平均流条件下，物质导数与空间散度可以交换：

$$
\mathcal D_0(\nabla\cdot\mathbf u')
=\nabla\cdot(\mathcal D_0\mathbf u').
$$

代入线性化动量方程，得到

$$
\frac1{\rho_0c_0^2}\mathcal D_0^2p'
-\nabla\cdot\left(\frac1{\rho_0}\nabla p'\right)=0,
$$

即紧凑的对流波动方程

$$
\frac1{\rho_0c_0^2}\mathcal D_0^2p'
=\nabla\cdot\left(\frac1{\rho_0}\nabla p'\right).
$$

**第五步：展开二阶物质导数。** 由 $\mathcal D_0$ 的定义，

$$
\begin{aligned}
\mathcal D_0^2p'
&=\left(\frac{\partial}{\partial t}
+\mathbf u_0\cdot\nabla\right)
\left(\frac{\partial p'}{\partial t}
+\mathbf u_0\cdot\nabla p'\right)\\
&=\frac{\partial^2p'}{\partial t^2}
+\frac{\partial}{\partial t}(\mathbf u_0\cdot\nabla p')
+\mathbf u_0\cdot\nabla\left(\frac{\partial p'}{\partial t}\right)
+\mathbf u_0\cdot\nabla(\mathbf u_0\cdot\nabla p').
\end{aligned}
$$

因为 $\mathbf u_0$ 定常且局部均匀，第二项与第三项相等，所以

$$
\mathcal D_0^2p'
=\frac{\partial^2p'}{\partial t^2}
+2\mathbf u_0\cdot\nabla\left(\frac{\partial p'}{\partial t}\right)
+\mathbf u_0\cdot\nabla(\mathbf u_0\cdot\nabla p').
$$

因此，平均流同时产生一个混合空间—时间导数项和一个二阶空间导数修正项。

**第六步：用背景马赫数矢量表示。** 利用 $\mathbf u_0=c_0\mathbf M$，并在局部单元内把 $c_0$ 和 $\mathbf M$ 视为常数，有

$$
\mathbf u_0\cdot\nabla\left(\frac{\partial p'}{\partial t}\right)
=c_0\frac{\partial}{\partial t}(\mathbf M\cdot\nabla p')
$$

以及

$$
\begin{aligned}
\mathbf u_0\cdot\nabla(\mathbf u_0\cdot\nabla p')
&=c_0\mathbf M\cdot\nabla
\left(c_0\mathbf M\cdot\nabla p'\right)\\
&=c_0^2\mathbf M\cdot\nabla(\mathbf M\cdot\nabla p').
\end{aligned}
$$

代回第四步的对流波动方程，得到

$$
\begin{aligned}
&\frac1{\rho_0c_0^2}\frac{\partial^2p'}{\partial t^2}
+\frac2{\rho_0c_0}\frac{\partial(\mathbf M\cdot\nabla p')}{\partial t}
+\frac1{\rho_0}\mathbf M\cdot\nabla(\mathbf M\cdot\nabla p')\\
&\qquad=\nabla\cdot\left(\frac1{\rho_0}\nabla p'\right).
\end{aligned}
$$

将混合导数项和二阶对流项移到右端，即得式（5.227）。

**第七步：静止介质与平面波检验。** 当 $\mathbf M=0$ 时，式（5.227）退化为

$$
\frac1{\rho_0c_0^2}p'_{tt}
=\nabla\cdot\left(\frac1{\rho_0}\nabla p'\right),
$$

对均匀介质进一步化为 $p'_{tt}=c_0^2\nabla^2p'$。对于均匀平均流中的平面波

$$
p'=\widehat p\,e^{\mathrm i(\mathbf k\cdot\mathbf x-\omega t)},
$$

第四步的对流波动方程给出频散关系

$$
(\omega-\mathbf u_0\cdot\mathbf k)^2=c_0^2|\mathbf k|^2,
\qquad
\omega=\mathbf u_0\cdot\mathbf k\pm c_0|\mathbf k|.
$$

这表明声波传播由介质中的声速和平均流携带速度共同决定，也验证了混合导数项和二阶对流项的符号。

式（5.227）不是任意可压缩流动的替代品。若平均流非均匀，完整的线性化动量方程还包含 $\rho_0(\mathbf u'\cdot\nabla)\mathbf u_0$ 和 $\rho'(\mathbf u_0\cdot\nabla)\mathbf u_0$，连续方程也包含 $\mathbf u'\cdot\nabla\rho_0$ 和 $\rho'\nabla\cdot\mathbf u_0$；若 $\rho_0$、$c_0$ 和 $\mathbf M$ 随空间变化，交换求导以及展开 $\mathbf u_0=c_0\mathbf M$ 时还会出现额外的系数梯度项。大振幅、强非均匀平均流、显著黏热效应或激波问题因此需要更完整的控制方程。即使在本章模型范围内，也应明确平均场是否定常、是否满足相应的连续方程，以及各项空间变系数在离散中的处理方式。

### 5.20.2 对流项的张量形式

若 $\mathbf M$ 在单元内为常数，或更一般地满足 $\nabla\cdot\mathbf M=0$，则纯空间算子满足

$$
\mathbf M\cdot\nabla(\mathbf M\cdot\nabla p')
=\nabla\cdot[(\mathbf M\otimes\mathbf M)\nabla p'],
\qquad
\mathbf M\otimes\mathbf M=
\begin{bmatrix}
M_x^2&M_xM_y&M_xM_z\\
M_yM_x&M_y^2&M_yM_z\\
M_zM_x&M_zM_y&M_z^2
\end{bmatrix}.
\tag{5.229}
$$

在单元内 $\rho_0$ 也为常数时，式（5.227）可写成较紧凑的局部形式：

$$
\frac1{\rho_0c_0^2}p'_{tt}
=\nabla\cdot\left(\frac1{\rho_0}\mathsf A\nabla p'\right)
-\frac2{\rho_0c_0}\frac{\partial(\mathbf M\cdot\nabla p')}{\partial t},
\qquad
\mathsf A=\mathsf I-\mathbf M\otimes\mathbf M.
\tag{5.230}
$$

矩阵 $\mathsf A$ 的特征值为 $1$、$1$ 和 $1-|\mathbf M|^2$。亚声速背景流 $|\mathbf M|<1$ 时它为正定矩阵；接近声速时，沿平均流方向的有效扩散趋近于零，离散和边界处理更加敏感。若 $\nabla\cdot\mathbf M\ne0$，直接使用式（5.230）会遗漏平均流散度相关项；若 $\rho_0$ 变化，把 $\rho_0^{-1}$ 移入散度也会产生额外梯度项。这些情况下应回到式（5.227）离散，或完整展开所有变系数项。

### 5.20.3 变量存储与声压插值

采用格点型有限体积法时：

- 节点存储声压 $p'_k$；
- 单元存储 $\rho_{0e}$、$c_{0e}$ 和 $\mathbf M_e$；
- 单元内声压梯度由节点声压和形函数梯度计算；
- 节点控制体由相邻单元中心、面中心和边中点连接形成。

单元内声压插值为

$$
p'_e(\mathbf x,t)=\sum_{a=1}^{n_e}N_a(\mathbf x)p'_a(t),
\qquad
\nabla p'_e=\sum_{a=1}^{n_e}p'_a\nabla N_a,
\tag{5.231}
$$

形函数及其雅可比变换与第5.14节相同。三节点三角形和四节点四面体给出常声压梯度，四边形、六面体和三棱柱通常需要在积分点重新计算梯度。

单位分解条件 $\sum_aN_a=1$ 保证常声压场的梯度为零。对每个单元检查 $\sum_a\nabla N_a=0$，可以在声学计算前发现节点方向、雅可比转置和单元畸变等错误。

## 5.21 格点型有限体积空间离散

### 5.21.1 节点控制体与集中声学质量

设 $\Omega_k$ 为节点 $k$ 的对偶控制体，$\mathcal E_k$ 为共享该节点的单元集合，$V_{k,e}=\Omega_k\cap\Omega_e$ 为单元 $e$ 对节点控制体的体积贡献。若把节点二阶时间导数作为控制体代表值，则左端离散为

$$
\int_{\Omega_k}\frac1{\rho_0c_0^2}
\frac{\partial^2p'}{\partial t^2}\,\mathrm d\Omega
\approx g_k\ddot p'_k,
\qquad
g_k=\sum_{e\in\mathcal E_k}
\frac{V_{k,e}}{\rho_{0e}c_{0e}^2}.
\tag{5.232}
$$

$g_k$ 是节点的集中声学属性，组装后形成对角矩阵 $\mathsf M$。二维问题中体积替换为控制体面积，并按模型需要乘以厚度。

四面体中位控制体常给每个节点贡献 $V_e/4$，规则六面体贡献 $V_e/8$；三角形和规则四边形分别贡献 $A_e/3$ 和 $A_e/4$。一般扭曲单元应直接计算对偶子区域体积或面积，并检查所有节点贡献之和等于原单元测度。

### 5.21.2 积分平衡与三类空间项

在 $\Omega_k$ 上积分式（5.227），得到

$$
g_k\ddot p'_k
=D_k-Q_k-H_k.
\tag{5.233}
$$

其中 $D_k$ 为扩散通量，$Q_k$ 为混合时间—空间项，$H_k$ 为二阶马赫数项。下面分别给出其离散形式。

对扩散项应用高斯定理：

$$
D_k=
\int_{\Omega_k}\nabla\cdot\left(\frac1{\rho_0}\nabla p'\right)\mathrm d\Omega
=\sum_{e\in\mathcal E_k}\int_{\Gamma_{k,e}}
\frac1{\rho_{0e}}\nabla p'_e\cdot\mathbf n\,\mathrm d\Gamma.
\tag{5.234}
$$

$\Gamma_{k,e}$ 是单元 $e$ 内部属于 $\Omega_k$ 的控制体边界。边界节点还需单独加入真实物理边界上的通量。对每个面片使用有向面积向量 $\mathbf A_s=A_s\mathbf n_s$，可将面通量写成 $\rho_{0e}^{-1}\nabla p'_e\cdot\mathbf A_s$。闭合控制体必须满足 $\sum_s\mathbf A_s=0$。

### 5.21.3 混合时间—空间项

若 $\rho_0$、$c_0$ 和 $\mathbf M$ 在单元内为常数，则

$$
Q_k=
2\int_{\Omega_k}\frac1{\rho_0c_0}
\frac{\partial(\mathbf M\cdot\nabla p')}{\partial t}\,\mathrm d\Omega
\approx
2\sum_{e\in\mathcal E_k}
\frac{V_{k,e}}{\rho_{0e}c_{0e}}
\mathbf M_e\cdot\frac{\partial\nabla p'_e}{\partial t}.
\tag{5.235}
$$

将式（5.231）代入后，$Q_k$ 只包含节点 $\dot p'$，因而组装到 $\mathsf C\dot{\mathbf p}$ 中。这里的 $\mathsf C$ 主要由对流产生，**不是物理阻尼矩阵**；实际黏热损耗、阻抗边界或人工阻尼应另外加入并明确标识。

如果背景场随时间变化，严格展开还会出现 $\partial(\rho_0^{-1}c_0^{-1}\mathbf M)/\partial t$ 等附加项。课件的形式对应背景系数在当前时间层冻结的处理。

### 5.21.4 二阶马赫数项

在单元内 $\mathbf M$ 为常数的条件下，利用式（5.229）和高斯定理：

$$
H_k=
\int_{\Omega_k}\frac1{\rho_0}
\mathbf M\cdot\nabla(\mathbf M\cdot\nabla p')\,\mathrm d\Omega
=\sum_{e\in\mathcal E_k}\int_{\Gamma_{k,e}}
\frac1{\rho_{0e}}
[(\mathbf M_e\otimes\mathbf M_e)\nabla p'_e]\cdot\mathbf n\,\mathrm d\Gamma.
\tag{5.236}
$$

展开后包含

$$
M_x^2p'_{xx}+M_y^2p'_{yy}+M_z^2p'_{zz}
+2M_xM_yp'_{xy}+2M_xM_zp'_{xz}+2M_yM_zp'_{yz}.
$$

矩阵形式避免了二阶交叉导数逐项展开的符号错误，也说明为什么含平均流时需要单元内声压梯度和控制体面法向。若不对 $\mathbf M$ 作单元常数假设，则还应保留体积中的变系数项，不能只把式（5.236）当作普遍恒等式。

### 5.21.5 节点离散方程与整体矩阵

将式（5.234）—（5.236）及式（5.231）代入，可得节点 $k$ 的完整离散方程：

$$
\boxed{
\begin{aligned}
g_k\ddot p'_k={}&
\sum_{e\in\mathcal E_k}\int_{\Gamma_{k,e}}
\frac{\nabla p'_e\cdot\mathbf n}{\rho_{0e}}\,\mathrm d\Gamma\\
&-2\sum_{e\in\mathcal E_k}
\frac{V_{k,e}}{\rho_{0e}c_{0e}}
\mathbf M_e\cdot\frac{\partial\nabla p'_e}{\partial t}\\
&-\sum_{e\in\mathcal E_k}\int_{\Gamma_{k,e}}
\frac{[(\mathbf M_e\otimes\mathbf M_e)\nabla p'_e]\cdot\mathbf n}{\rho_{0e}}
\,\mathrm d\Gamma
+\mathcal S_k.
\end{aligned}}
\tag{5.237}
$$

其中 $\mathcal S_k$ 表示体积声源和物理边界的附加贡献。二维时，$V_{k,e}$ 换成面积 $S_{k,e}$，体积分换成面积分，面积分换成线积分。

把所有节点未知量排列成 $\mathbf p$，并将时间导数项移到左端，得到

$$
\boxed{
\mathsf M\ddot{\mathbf p}
+\mathsf C\dot{\mathbf p}
+\mathsf K\mathbf p
=\mathbf F(t)
}.
\tag{5.238}
$$

$\mathsf K$ 包含 $H_k-D_k$ 的系数，$\mathsf C$ 包含 $Q_k$ 的系数，$\mathbf F$ 包含声源和外边界输入。静止均匀介质且无物理阻尼时，$\mathsf C=0$，一致离散应恢复经典声波方程的对称空间算子。

## 5.22 边界条件与声源处理

### 5.22.1 绝对硬边界与振动壁面

静止绝对硬边界的法向质点速度为零，因此

$$
\nabla p'\cdot\mathbf n=\frac{\partial p'}{\partial n}=0.
\tag{5.239}
$$

在静止各向同性介质中，这对应齐次 Neumann 条件，边界声学通量为零。若壁面以法向加速度 $a_n$ 振动，则由线性 Euler 方程有

$$
\frac{\partial p'}{\partial n}=-\rho_0a_n.
\tag{5.240}
$$

该条件把结构振动输入声场，是结构声耦合的基础。对于含平均流的完整算子，法向通量还包含 $\mathbf M\otimes\mathbf M$ 和切向导数，不能在所有情况下都把两个边界积分同时简单置零。

### 5.22.2 阻抗边界与静止介质吸收边界

局部反应阻抗边界以 $p'=Z\mathbf u'\cdot\mathbf n$ 联系声压和法向质点速度。若 $Z$ 为与时间无关的常数，则

$$
\frac{\partial p'}{\partial n}
=-\frac{\rho_0}{Z}\frac{\partial p'}{\partial t}.
\tag{5.241}
$$

当 $Z=\rho_0c_0$ 时得到平面波吸收边界：

$$
\boxed{
\frac1{c_0}\frac{\partial p'}{\partial t}
+\frac{\partial p'}{\partial n}=0
}.
\tag{5.242}
$$

它只精确吸收垂直入射的平面波，对斜入射、近场和非平面波仍有反射。频变阻抗应在频域使用，或通过合适的时域卷积、辅助变量模型实现。

将式（5.242）代入静止介质边界扩散通量，得到

$$
\int_{\Gamma_b}\frac1{\rho_0}
\frac{\partial p'}{\partial n}\,\mathrm d\Gamma
=-\int_{\Gamma_b}\frac1{\rho_0c_0}
\frac{\partial p'}{\partial t}\,\mathrm d\Gamma.
\tag{5.243}
$$

因此吸收边界在 $\mathsf C$ 中增加真正具有耗散作用的正贡献。对于课件中坐标面与平均流方向对齐的简化情况，二阶马赫数边界项也可用 $M_n^2c_0^{-1}\partial p'/\partial n$ 表示；合并后边界系数包含 $1-M_n^2$。若存在切向平均流和非零切向声压梯度，则应保留完整的 $\mathsf A\nabla p'\cdot\mathbf n$，不能直接套用该标量系数。

### 5.22.3 平均流出流、压力边界与声源

平均流存在时，一维出射特征可近似写为

$$
\frac1{U_n+c_0}\frac{\partial p'}{\partial t}
+\frac{\partial p'}{\partial n}=0,
\tag{5.244}
$$

其中 $U_n=\mathbf u_0\cdot\mathbf n$。这里 $U_n+c_0$ 是沿外法向的出射声速，必须为正。多维出流还涉及涡波、熵波和不同入射角，宜采用第6章讨论的特征边界、Tam 条件或 PML。

指定声压 $p'=\overline p'_b$ 属于 Dirichlet 条件，可通过自由度消去施加；单极子、偶极子或分布声源则进入 $\mathbf F$。源项离散应检查总源强、空间归一化和时间带宽。把点源简单放入单个节点时，其幅值会随网格尺度改变，应按控制体体积或规定的源分布进行一致归一化。

边界条件必须进入与其物理含义相同的通量项。仅在求解完成后覆盖边界节点声压，不能等价实现吸收边界或振动壁面。

## 5.23 时间积分方法

### 5.23.1 显式中心差分

二阶时间导数采用中心差分：

$$
\ddot{\mathbf p}^{\,n}
\approx\frac{\mathbf p^{n+1}-2\mathbf p^n+\mathbf p^{n-1}}{\Delta t^2}.
\tag{5.245}
$$

一阶时间导数可使用二阶中心差分

$$
\dot{\mathbf p}^{\,n}
\approx\frac{\mathbf p^{n+1}-\mathbf p^{n-1}}{2\Delta t},
\tag{5.246}
$$

也可对混合项使用课件中的后差形式

$$
\left.\frac{\partial(\mathbf M\cdot\nabla p')}{\partial t}\right|^n
\approx\frac{(\mathbf M\cdot\nabla p')^n-(\mathbf M\cdot\nabla p')^{n-1}}{\Delta t}.
\tag{5.247}
$$

后差只有一阶时间精度，与二阶空间项和二阶时间导数组合后可能降低整体精度；若目标是二阶声学推进，宜优先采用相容的中心形式或经过验证的其他二阶格式。

将式（5.245）和式（5.246）代入式（5.238），得到

$$
\left(\frac{\mathsf M}{\Delta t^2}+\frac{\mathsf C}{2\Delta t}\right)\mathbf p^{n+1}
=\mathbf F^n-\mathsf K\mathbf p^n
+\frac{2\mathsf M}{\Delta t^2}\mathbf p^n
+\left(\frac{\mathsf C}{2\Delta t}-\frac{\mathsf M}{\Delta t^2}\right)\mathbf p^{n-1}.
\tag{5.248}
$$

当 $\mathsf C=0$ 且 $\mathsf M$ 为对角矩阵时，可直接逐节点更新，形成完全显式格式。若存在对流混合项、吸收边界或一致阻尼矩阵，即使时间导数采用中心差分，也可能需要求解一个线性系统。

初始步还需 $\mathbf p^{-1}$。可根据初始 $\dot{\mathbf p}^0$ 和初始加速度构造相容起步值，或用一阶起步后尽快切换到二阶格式；不应不加说明地令 $\mathbf p^{-1}=\mathbf p^0$，因为这相当于额外假设初始速度为零。

### 5.23.2 Newmark 方法

令 $\mathbf v_p=\dot{\mathbf p}$、$\mathbf a_p=\ddot{\mathbf p}$。Newmark 假设为

$$
\mathbf p_{n+1}=\mathbf p_n+\Delta t\mathbf v_{p,n}
+\Delta t^2\left[\left(\frac12-\beta\right)\mathbf a_{p,n}
+\beta\mathbf a_{p,n+1}\right],
\tag{5.249}
$$

$$
\mathbf v_{p,n+1}=\mathbf v_{p,n}
+\Delta t[(1-\gamma)\mathbf a_{p,n}+\gamma\mathbf a_{p,n+1}].
\tag{5.250}
$$

定义系数

$$
\begin{aligned}
&a_0=\frac1{\beta\Delta t^2},&&a_1=\frac\gamma{\beta\Delta t},\\
&a_2=\frac1{\beta\Delta t},&&a_3=\frac1{2\beta}-1,\\
&a_4=\frac\gamma\beta-1,&&a_5=\frac{\Delta t}{2}\left(\frac\gamma\beta-2\right).
\end{aligned}
\tag{5.251}
$$

则新时间层的二阶和一阶导数为

$$
\mathbf a_{p,n+1}
=a_0(\mathbf p_{n+1}-\mathbf p_n)-a_2\mathbf v_{p,n}-a_3\mathbf a_{p,n},
\tag{5.252}
$$

$$
\mathbf v_{p,n+1}
=a_1(\mathbf p_{n+1}-\mathbf p_n)-a_4\mathbf v_{p,n}-a_5\mathbf a_{p,n}.
\tag{5.253}
$$

代入 $\mathsf M\mathbf a_{p,n+1}+\mathsf C\mathbf v_{p,n+1}+\mathsf K\mathbf p_{n+1}=\mathbf F_{n+1}$，得到以 $\mathbf p_{n+1}$ 为未知量的系统

$$
\widehat{\mathsf K}\mathbf p_{n+1}=\widehat{\mathbf F}_{n+1},
\tag{5.254}
$$

其中

$$
\widehat{\mathsf K}=\mathsf K+a_0\mathsf M+a_1\mathsf C,
\tag{5.255}
$$

$$
\begin{aligned}
\widehat{\mathbf F}_{n+1}={}&\mathbf F_{n+1}
+\mathsf M(a_0\mathbf p_n+a_2\mathbf v_{p,n}+a_3\mathbf a_{p,n})\\
&+\mathsf C(a_1\mathbf p_n+a_4\mathbf v_{p,n}+a_5\mathbf a_{p,n}).
\end{aligned}
\tag{5.256}
$$

求解后，再用式（5.252）和式（5.253）更新两个时间导数。初始二阶导数应由初始平衡计算：

$$
\mathsf M\mathbf a_{p,0}
=\mathbf F_0-\mathsf C\mathbf v_{p,0}-\mathsf K\mathbf p_0.
\tag{5.257}
$$

$\gamma=1/2$、$\beta=1/4$ 为平均加速度法。对于对称正定的 $\mathsf M$、$\mathsf K$ 和非负阻尼矩阵的标准线性系统，它无条件稳定且二阶精度；含非对称对流矩阵时仍应检查完整离散算子的谱性质。$\gamma=1/2$、$\beta=0$ 时 Newmark 关系退化为显式中心差分形式，时间步受到条件稳定限制。

无条件稳定不代表时间步可以任意大。声周期内仍需足够采样点，平均流、边界阻抗和最高关注频率都应纳入时间精度评估。

### 5.23.3 频域 Helmholtz 形式

若所有系数不随时间变化，且声源和响应均为稳态简谐量，采用 $e^{\mathrm i\omega t}$ 约定可得

$$
(-\omega^2\mathsf M+\mathrm i\omega\mathsf C+\mathsf K)
\widehat{\mathbf p}=\widehat{\mathbf F}.
\tag{5.258}
$$

频域法不需要逐步推进，适合扫频、阻抗边界和模态分析，但每个频率都要求解复数稀疏系统。平均流和阻抗可能使矩阵非对称、非 Hermitian，不能直接套用仅适用于实对称正定矩阵的求解器。

## 5.24 稳定性、数值色散与验证

### 5.24.1 正三角形网格上的色散分析

对静止均匀介质中的二维正三角形网格，课件采用的中心差分离散可写为

$$
\frac{p_{i,j}^{n+1}-2p_{i,j}^n+p_{i,j}^{n-1}}{\Delta t^2}
=\frac{2c_0^2}{3\Delta x^2}
\left(\sum_{q=1}^{6}p_q^n-6p_{i,j}^n\right),
\tag{5.259}
$$

其中 $p_q$ 为六个相邻节点声压，$\Delta x$ 为正三角形边长。设平面波

$$
p_{i,j}^n=p_0\exp[\mathrm i(\omega n\Delta t-k_xx_i-k_yy_j)],
\qquad
(\omega/c_0)^2=k_x^2+k_y^2.
\tag{5.260}
$$

代入后得到数值色散关系

$$
\sin^2\left(\frac{\omega\Delta t}{2}\right)
=\frac{c_0^2\Delta t^2}{6\Delta x^2}
\left[
8-\left(2\cos\frac{k_x\Delta x}{2}
+\cos\frac{\sqrt3k_y\Delta x}{2}\right)^2
+\cos^2\frac{\sqrt3k_y\Delta x}{2}
\right].
\tag{5.261}
$$

当 $\Delta x\to0$、$\Delta t\to0$ 时，它恢复连续色散关系 $\omega^2=c_0^2(k_x^2+k_y^2)$。稳定性要求右端对所有可解析波数都不超过1；由括号项的最大值可得到该特定格式的条件

$$
\frac{c_0\Delta t}{\Delta x}\leq\sqrt{\frac23}.
\tag{5.262}
$$

这一推导同时说明：稳定时间步来自**整个空间—时间离散算子的谱范围**，不能仅根据波速和单元类型凭经验确定。

### 5.24.2 CFL 参考值与平均流修正

课件针对相应规则网格和中心差分格式给出以下参考值：

| 单元类型 | 稳定性参考条件 |
|---|---|
| 二维正三角形 | $c_0\Delta t/\Delta x\leq\sqrt{2/3}$ |
| 二维规则四边形 | $c_0\Delta t/\Delta x\leq1$ |
| 正四面体 | $c_0\Delta t/\Delta x\leq1/\sqrt2$ |
| 规则三棱柱 | $c_0\Delta t/\Delta x\leq1/\sqrt2$ |
| 规则六面体 | $c_0\Delta t/\Delta x\leq1$ |

一般可写成

$$
Co=\frac{c_{\mathrm{eff}}\Delta t}{h_{\mathrm{eff}}}
\leq\eta Co_{\max},
\tag{5.263}
$$

其中 $\eta<1$ 为安全裕量。含平均流时，最大特征传播速度至少应考虑 $|\mathbf u_0|+c_0$，而不能只使用 $c_0$。扭曲网格、非均匀介质、混合单元、PML 和耦合边界还会改变有效尺度及稳定上限。最可靠的做法是对实际离散矩阵进行谱估计，并通过数值试验验证。

### 5.24.3 时间精度、每波节点数与相位误差

CFL 只回答“是否稳定”，不回答“传播是否准确”。空间分辨率可用每波节点数

$$
N_\lambda=\frac{\lambda}{h}
=\frac{2\pi}{kh}
\tag{5.264}
$$

表示。低阶格式通常可先以 $N_\lambda\geq10$ 作为起点，再根据传播距离、允许相位误差和格式色散曲线增加分辨率。高阶低色散格式可以降低所需节点数，但对网格质量和边界实现更敏感。

时间分辨率可用 $\omega\Delta t$ 或每周期的采样点数衡量。若数值角频率为 $\omega_h$，相对相位速度误差可写为

$$
e_p=\frac{\omega_h}{kc_0}-1.
\tag{5.265}
$$

经过传播时间 $T$ 后，累计相位差约为 $(\omega_h-\omega)T$。因此，短脉冲局部波形正确并不能保证长时间、多波长传播后的远场相位正确。

### 5.24.4 完整求解流程

```text
读取节点、单元、边界、rho_0、c_0、u_0 和声源
构造节点对偶控制体，计算体积、面积向量和集中属性 g_k
构建形函数及整体坐标梯度，检查 sum(grad N_a)=0
逐单元计算扩散、混合时间项和二阶马赫数项系数
装配 M、C、K 和体积声源向量
施加硬壁、振动壁、阻抗、吸收或压力边界
选择显式中心差分、Newmark 或频域求解
每个时间步更新声压及时间导数，检查残差和能量
进行网格、时间步、传播距离和人工边界敏感性分析
输出声压、质点振速、频谱和关注位置的时间历程
```

若 $\rho_0$、$c_0$ 或 $\mathbf M$ 随背景流场更新，应在外迭代中重新组装受影响系数。不能在材料或平均流改变后仍使用旧的 $\mathsf M$、$\mathsf C$ 和 $\mathsf K$。

### 5.24.5 分层验证与常见错误

建议按以下顺序验证：

1. **几何与插值**：控制体体积剖分、闭合面积向量、形函数单位分解和常声压零梯度；
2. **静止波动解**：一维或二维平面波的波速、频率和幅值；
3. **硬壁反射**：刚性壁面反射系数和驻波节点位置；
4. **吸收边界**：不同入射角下的反射，以及增大计算域的参考解；
5. **平均流对流**：顺流和逆流声速、波长及群速度；
6. **时间收敛**：统一空间网格下减小时间步，区分时间误差与空间误差；
7. **长距离色散**：比较多个传播距离后的相位和包络；
8. **矩阵与能量**：检查 $\mathsf M$ 对角质量、边界阻尼符号和无耗散系统的能量变化。

归一化离散残差可写为

$$
\eta_R=
\frac{\|\mathbf F-\mathsf M\ddot{\mathbf p}
-\mathsf C\dot{\mathbf p}-\mathsf K\mathbf p\|_2}{F_{\mathrm{ref}}},
\tag{5.266}
$$

同时还应检查声压增量、界面通量和目标测点响应。静止、无源、无耗散声场的离散能量可用

$$
E_a=\frac12\int_\Omega
\left(\frac{p'^2}{\rho_0c_0^2}+\rho_0|\mathbf u'|^2\right)\mathrm d\Omega
\tag{5.267}
$$

构造参考指标。显式格式在稳定范围内仍可能出现周期性能量振荡，应区分正常的时间离散误差与失稳增长。

常见实现错误包括：

- 把马赫数矢量 $\mathbf M$ 与质量矩阵 $\mathsf M$ 混用；
- 在 $\nabla\cdot\mathbf M\ne0$ 时仍直接使用式（5.230）；
- 形函数梯度遗漏雅可比转置，或单元节点方向不一致；
- 将混合对流矩阵 $\mathsf C$ 当成必然耗散的物理阻尼；
- 吸收边界时间项符号错误，使耗散变为能量输入；
- 含切向平均流时遗漏 $\mathsf A\nabla p'$ 中的交叉导数；
- 二阶时间导数用中心差分，却用一阶后差处理其他时间项并仍声称整体二阶；
- $\beta=0$ 时忽略非对角 $\mathsf C$，误以为不需要求解线性系统；
- 只按静止声速选择含平均流问题的时间步；
- 只检查是否发散，不检查相位、反射和网格收敛。

## 本章小结

格点型声学有限体积法把声压存储在节点，把背景介质参数存储在单元，通过对偶控制体积分和单元形函数梯度构造离散方程。集中声学质量为 $g_k=\sum_eV_{k,e}/(\rho_{0e}c_{0e}^2)$；扩散项和二阶马赫数项转化为控制体面通量，混合时间—空间项组装为一阶时间导数矩阵。

静止硬壁对应零法向声压梯度，振动壁由法向加速度给定 Neumann 通量，阻抗和静止介质吸收边界产生时间导数边界项。含平均流和切向导数时，应保留完整张量通量，不能直接套用坐标对齐的简化公式。

显式中心差分简单高效，但受 CFL 条件限制；Newmark 方法可统一处理一阶和二阶时间导数，并转化为新时间层声压的有效刚度系统。正三角形网格的色散关系给出了特定稳定上限，其他单元参考值只适用于相应规则网格和格式。可靠声学计算还必须检查每波节点数、相位误差、边界反射和网格—时间步收敛。

## 习题

1. 从式（5.227）出发，说明 $\mathbf M=0$ 时如何恢复经典波动方程。
2. 在 $\nabla\cdot\mathbf M=0$ 条件下推导式（5.229），并列出变马赫数时可能遗漏的项。
3. 证明 $\mathsf A=\mathsf I-\mathbf M\otimes\mathbf M$ 的特征值为 $1$、$1$ 和 $1-|\mathbf M|^2$，讨论近声速时的数值困难。
4. 推导节点集中声学属性式（5.232），并验证其量纲与扩散通量项相容。
5. 将式（5.231）代入式（5.235），写出节点 $k$ 对任意邻节点 $j$ 的一阶时间导数系数。
6. 展开式（5.236），验证六个二阶导数项及三个交叉项的系数。
7. 推导整体矩阵 $\mathsf M\ddot{\mathbf p}+\mathsf C\dot{\mathbf p}+\mathsf K\mathbf p=\mathbf F$ 中各项与 $D_k$、$Q_k$、$H_k$ 的对应关系。
8. 从线性 Euler 方程推导振动壁条件式（5.240），说明法向符号。
9. 从 $p'=Z u_n'$ 推导阻抗条件式（5.241），并证明 $Z=\rho_0c_0$ 对应平面波吸收边界。
10. 推导式（5.248），讨论 $\mathsf C$ 非对角时显式格式为什么仍需线性求解。
11. 从 Newmark 假设推导式（5.252）—（5.256），并说明初始二阶导数如何确定。
12. 推导频域方程式（5.258），讨论对流和阻抗对矩阵性质的影响。
13. 将平面波式（5.260）代入正三角形离散式（5.259），推导色散关系式（5.261）。
14. 证明式（5.261）在时间和空间步长趋近于零时恢复 $\omega^2=c_0^2(k_x^2+k_y^2)$。
15. 说明为什么 $Co<1$ 不是所有网格和时间格式的普适稳定条件，并比较表中各单元的参考上限。
16. 对含平均流的介质，为什么稳定时间步应考虑 $|\mathbf u_0|+c_0$？这与频域绝对稳定性有何区别？
17. 设计一个平面波数值色散试验，分别估计空间步长和时间步长对相位误差的贡献。
18. 设计吸收边界验证方案，至少包含正入射、斜入射和更大计算域参考解。
19. 为什么静止无耗散系统中的能量振荡不一定意味着失稳？应如何区分误差与指数增长？
20. 编写声学求解器的实现检查表，覆盖几何、形函数、矩阵装配、边界、时间推进、稳定性和色散验证。

## 延伸阅读

1. Tam, C. K. W., and Webb, J. C. Dispersion-relation-preserving finite difference schemes for computational acoustics. *Journal of Computational Physics*, 1993, 107: 262–281.
2. Tam, C. K. W. Computational aeroacoustics: issues and methods. *AIAA Journal*, 1995, 33: 1788–1796.
3. Newmark, N. M. A method of computation for structural dynamics. *Journal of the Engineering Mechanics Division, ASCE*, 1959, 85(3): 67–94.
4. Hughes, T. J. R. *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*. Dover, 2000.
5. Eymard, R., Gallouët, T., and Herbin, R. Finite volume methods. In *Handbook of Numerical Analysis*, Vol. 7, North-Holland, 2000: 713–1020.
6. Berenger, J. P. A perfectly matched layer for the absorption of electromagnetic waves. *Journal of Computational Physics*, 1994, 114: 185–200.
7. Fahy, F. J., and Gardonio, P. *Sound and Structural Vibration: Radiation, Transmission and Response*. 2nd ed., Academic Press, 2007.
