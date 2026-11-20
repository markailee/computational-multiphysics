# 第5章（续）结构静力学与动力学：结构场的有限体积离散

结构场求解的基本未知量是位移。由位移可进一步计算应变、应力和节点反力，并与流体压力、温度场或其他物理场交换数据。与流场不同，结构方程中的未知量通常定义在网格节点上，而密度、弹性模量、泊松比和体积力等材料或载荷参数可定义在单元内。因此，结构离散需要同时处理**节点控制体上的力平衡**和**单元内部的位移插值**。

本章以课件采用的格点型有限体积法为主线：围绕节点构造控制体，对控制体积分动量平衡方程，利用高斯定理把应力散度转化为控制体边界上的面力通量，再用单元形函数计算应变和应力。内容按“控制方程—形函数与应变—二维离散—三维推广—边界条件与装配—静力及动力求解”的顺序展开，节号和公式编号承接本章前两篇。除特别说明外，假设结构为小变形、线弹性体，单元和节点编号均为局部编号。

## 5.13 结构控制方程与离散框架

### 5.13.1 运动学、本构关系与力平衡

设结构位移矢量为 $\mathbf u(\mathbf x,t)$。在小变形条件下，应变张量是位移梯度的对称部分：

$$
\boldsymbol\varepsilon
=\frac12\left(\nabla\mathbf u+\nabla\mathbf u^{\mathrm T}\right).
\tag{5.163}
$$

线弹性本构关系为

$$
\boldsymbol\sigma=\mathsf C:\boldsymbol\varepsilon,
\qquad
\sigma_{ij}=C_{ijkl}\varepsilon_{kl},
\tag{5.164}
$$

其中 $\mathsf C$ 为四阶弹性张量。对于各向同性材料，只需两个独立材料参数，例如弹性模量 $E$ 和泊松比 $\nu$，或 Lamé 参数 $\lambda$ 和 $\mu$。

结构的局部动量平衡方程为

$$
\rho\frac{\partial^2\mathbf u}{\partial t^2}
=\nabla\cdot\boldsymbol\sigma+\mathbf b
\qquad \text{在 }\Omega\times(0,T]\text{ 内},
\tag{5.165}
$$

其中 $\rho$ 为密度，$\mathbf b$ 为单位体积体力，例如重力、离心力或电磁力。惯性项 $\rho\ddot{\mathbf u}$ 来自结构自身的加速度；在流固耦合等移动结构问题中，还应确保位移、速度和加速度采用一致的时间参考系。

位移和力边界条件分别为

$$
\mathbf u=\overline{\mathbf u}
\quad\text{在 }\Gamma_u\text{ 上},
\qquad
\boldsymbol\sigma\mathbf n=\overline{\mathbf t}
\quad\text{在 }\Gamma_t\text{ 上},
\qquad
\Gamma_u\cup\Gamma_t=\partial\Omega,
\qquad
\Gamma_u\cap\Gamma_t=\varnothing,
\tag{5.166}
$$

其中 $\overline{\mathbf u}$ 为指定位移，$\overline{\mathbf t}$ 为指定面力，$\mathbf n$ 为结构边界单位外法线。动力问题还需给定初始位移和初始速度：

$$
\mathbf u(\mathbf x,0)=\mathbf u_0(\mathbf x),
\qquad
\dot{\mathbf u}(\mathbf x,0)=\mathbf v_0(\mathbf x).
\tag{5.167}
$$

静力问题是惯性可以忽略或结构已达到静力平衡时的特例：

$$
\nabla\cdot\boldsymbol\sigma+\mathbf b=\mathbf 0.
\tag{5.168}
$$

“稳态动力响应”与静力问题并不相同。若结构在周期载荷下达到稳定的简谐响应，惯性项仍然存在，只是各自由度的相位和幅值不再随平均时间变化，需用频域方程或足够长的瞬态计算获得。

### 5.13.2 待解变量与材料参数的存储位置

结构有限体积离散可采用两类节点布置。

- **单元中心法**：位移等未知量存储在单元中心，控制体就是网格单元。单元面两侧位移和材料不连续时，需要重构界面位移梯度及面力通量。
- **格点法或节点中心法**：位移存储在网格节点，围绕节点重新构造对偶控制体。材料参数和单元应力通常存储在原始单元内。该方法也称顶点中心有限体积法，是本章采用的方法。

格点法中，一个原始网格单元会被其边中点、面中心和单元中心分割，并分别贡献给相邻节点的控制体。于是出现两套几何对象：

1. **原始单元**：用于定义形函数、材料参数、应变和应力；
2. **节点控制体**：用于列写惯性、应力和外载荷的积分平衡。

位移在节点上天然连续，而应力只要求在单元内部按本构关系计算。跨过原始单元边界时，数值应力可以发生跳变，这与低阶位移单元的精度特征一致。

### 5.13.3 与有限元方法的联系和区别

本方法借用有限元的形函数计算单元内位移和应变，但控制方程的离散方式仍是有限体积积分平衡，而不是 Galerkin 弱形式。两者不能直接混称。

| 项目 | 格点型有限体积法 | 位移有限元法 |
|---|---|---|
| 基本未知量 | 节点位移 | 节点位移 |
| 方程来源 | 节点控制体积分力平衡 | 加权余量或虚功原理 |
| 应力项 | 对偶控制体边界的面力通量 | 单元域内的 $\mathbf B^{\mathrm T}\boldsymbol\sigma$ 积分 |
| 质量矩阵 | 常按控制体质量集中为对角矩阵 | 可采用一致质量或集中质量 |
| 局部守恒 | 在节点控制体上直接表达 | 通过弱形式及后处理分析 |

对某些规则网格、特殊积分和材料假设，两种离散可能得到相同或相近的系数；一般情况下，**控制体面通量矩阵不能未经推导就替换为有限元的 $\int\mathbf B^{\mathrm T}\mathbf D\mathbf B\,\mathrm d\Omega$**。程序验证时应以所采用的离散公式为准。

## 5.14 形函数、几何映射与应变计算

### 5.14.1 位移插值与形函数的基本性质

设单元有 $n_e$ 个节点，节点 $a$ 的形函数为 $N_a$。单元内位移插值为

$$
\mathbf u^e(\mathbf x,t)
=\sum_{a=1}^{n_e}N_a(\mathbf x)\mathbf u_a(t),
\qquad
u_i^e=\sum_{a=1}^{n_e}N_a u_{a,i}.
\tag{5.169}
$$

形函数应满足节点插值条件和单位分解性：

$$
N_a(\mathbf x_b)=\delta_{ab},
\qquad
\sum_{a=1}^{n_e}N_a(\mathbf x)=1,
\qquad
\sum_{a=1}^{n_e}\nabla N_a(\mathbf x)=\mathbf 0.
\tag{5.170}
$$

前两条保证位移在节点处等于节点值，并能精确表示刚体平移和常应变状态；第三条是计算梯度时的重要检查条件。等参单元还采用同一组形函数插值几何坐标，因此也能精确表示单元的刚体运动。

### 5.14.2 三节点三角形单元

对逆时针编号的三角形节点 $(x_1,y_1)$、$(x_2,y_2)$、$(x_3,y_3)$，其有符号面积为

$$
A_e=\frac12
\begin{vmatrix}
1&x_1&y_1\\
1&x_2&y_2\\
1&x_3&y_3
\end{vmatrix}
=\frac12[(x_2-x_1)(y_3-y_1)-(x_3-x_1)(y_2-y_1)].
\tag{5.171}
$$

假设 $u=a_1+a_2x+a_3y$，将三个节点的位移代入并求解系数，可得线性形函数

$$
N_a=\frac{\alpha_a+\beta_ax+\gamma_ay}{2A_e},
\qquad a=1,2,3,
\tag{5.172}
$$

其中

$$
\begin{aligned}
&\alpha_1=x_2y_3-x_3y_2,&&\beta_1=y_2-y_3,&&\gamma_1=x_3-x_2,\\
&\alpha_2=x_3y_1-x_1y_3,&&\beta_2=y_3-y_1,&&\gamma_2=x_1-x_3,\\
&\alpha_3=x_1y_2-x_2y_1,&&\beta_3=y_1-y_2,&&\gamma_3=x_2-x_1.
\end{aligned}
$$

因此

$$
\frac{\partial N_a}{\partial x}=\frac{\beta_a}{2A_e},
\qquad
\frac{\partial N_a}{\partial y}=\frac{\gamma_a}{2A_e}.
\tag{5.173}
$$

形函数在整个单元内为线性函数，其梯度为常数，所以三节点三角形是常应变三角形（CST）。若节点顺时针编号，则式（5.171）给出负面积。程序可以统一交换两个节点使其变为逆时针，而不应简单地在梯度公式中混用绝对面积和原始有向面积。

### 5.14.3 四边形、四面体和六面体单元

四节点双线性四边形在参考坐标 $(\xi,\eta)\in[-1,1]^2$ 中的形函数为

$$
N_a(\xi,\eta)=\frac14(1+\xi_a\xi)(1+\eta_a\eta),
\qquad a=1,\ldots,4,
\tag{5.174}
$$

其中 $(\xi_a,\eta_a)$ 依次为 $(-1,-1)$、$(1,-1)$、$(1,1)$、$(-1,1)$。对于非矩形四边形，形函数梯度随位置变化，需在积分点计算。

四节点四面体可采用参考坐标中的体积坐标

$$
N_1=1-\xi-\eta-\zeta,
\qquad N_2=\xi,
\qquad N_3=\eta,
\qquad N_4=\zeta.
\tag{5.175}
$$

其整体坐标梯度为常数，因而也是常应变单元。实际实现可构造矩阵 $[1\ x_a\ y_a\ z_a]$ 并求逆，一次获得四个形函数的常数项和梯度，避免逐项抄写冗长的行列式公式。

八节点三线性六面体的形函数为

$$
N_a(\xi,\eta,\zeta)
=\frac18(1+\xi_a\xi)(1+\eta_a\eta)(1+\zeta_a\zeta),
\qquad a=1,\ldots,8,
\tag{5.176}
$$

其中 $(\xi_a,\eta_a,\zeta_a)$ 是参考立方体的八个顶点。三棱柱等单元可分别采用三角形面积坐标与一维形函数的乘积。不同文献的节点顺序可能不同，程序必须先固定节点编号约定，再配套使用形函数、面定义和控制体拓扑。

### 5.14.4 等参映射与整体坐标梯度

等参单元的几何映射写为

$$
\mathbf x(\boldsymbol\xi)=\sum_{a=1}^{n_e}N_a(\boldsymbol\xi)\mathbf x_a,
\qquad
\boldsymbol\xi=(\xi,\eta)\ \text{或}\ (\xi,\eta,\zeta).
\tag{5.177}
$$

定义雅可比矩阵

$$
\mathbf J=\frac{\partial\mathbf x}{\partial\boldsymbol\xi},
\qquad
J_{ij}=\frac{\partial x_i}{\partial\xi_j}
=\sum_{a=1}^{n_e}x_{a,i}\frac{\partial N_a}{\partial\xi_j}.
\tag{5.178}
$$

由链式法则 $\nabla_{\boldsymbol\xi}N_a=\mathbf J^{\mathrm T}\nabla_{\mathbf x}N_a$，得到

$$
\nabla_{\mathbf x}N_a
=\mathbf J^{-\mathrm T}\nabla_{\boldsymbol\xi}N_a.
\tag{5.179}
$$

式（5.179）是形函数导数计算的核心。实现时应解一个小的线性系统，而不是显式写出矩阵逆。若 $\det\mathbf J=0$，映射退化；若在单元内变号，则单元翻转或严重畸变，不能继续沿用正常积分。

参考域积分转换为整体坐标积分：

$$
\int_{\Omega_e}\phi(\mathbf x)\,\mathrm d\Omega
=\int_{\widehat\Omega}\phi(\boldsymbol\xi)
|\det\mathbf J|\,\mathrm d\widehat\Omega
\approx\sum_{q=1}^{n_q}w_q|\det\mathbf J_q|\phi_q.
\tag{5.180}
$$

二维平面问题还需乘以厚度 $h$；轴对称问题还需包含径向权重并增加环向应变，不能直接套用平面公式。

### 5.14.5 应变—位移矩阵

二维工程记法取

$$
\boldsymbol\varepsilon_v=
\begin{bmatrix}\varepsilon_{xx}&\varepsilon_{yy}&\gamma_{xy}\end{bmatrix}^{\mathrm T},
\qquad
\gamma_{xy}=2\varepsilon_{xy},
$$

以及

$$
\boldsymbol\sigma_v=
\begin{bmatrix}\sigma_{xx}&\sigma_{yy}&\sigma_{xy}\end{bmatrix}^{\mathrm T}.
$$

对每个节点定义

$$
\mathbf B_a=
\begin{bmatrix}
N_{a,x}&0\\
0&N_{a,y}\\
N_{a,y}&N_{a,x}
\end{bmatrix},
\qquad
\mathbf B_e=
\begin{bmatrix}\mathbf B_1&\cdots&\mathbf B_{n_e}\end{bmatrix},
\qquad
\boldsymbol\varepsilon_v=\mathbf B_e\mathbf d_e,
\tag{5.181}
$$

其中 $\mathbf d_e=[u_{1x},u_{1y},\ldots,u_{n_ex},u_{n_ey}]^{\mathrm T}$，$N_{a,x}=\partial N_a/\partial x$，其余记号类推。工程剪应变与弹性剪应力配套使用，不能再额外乘二。

各向同性材料的 Lamé 参数为

$$
\mu=\frac{E}{2(1+\nu)},
\qquad
\lambda=\frac{E\nu}{(1+\nu)(1-2\nu)}.
\tag{5.182}
$$

平面应变时采用 $\lambda$；平面应力时使用等效参数

$$
\lambda^*=\frac{2\mu\lambda}{\lambda+2\mu}
=\frac{E\nu}{1-\nu^2}.
\tag{5.183}
$$

二维弹性矩阵可统一写为

$$
\mathbf D=
\begin{bmatrix}
2\mu+\bar\lambda&\bar\lambda&0\\
\bar\lambda&2\mu+\bar\lambda&0\\
0&0&\mu
\end{bmatrix},
\qquad
\boldsymbol\sigma_v=\mathbf D\boldsymbol\varepsilon_v,
\tag{5.184}
$$

其中平面应变取 $\bar\lambda=\lambda$，平面应力取 $\bar\lambda=\lambda^*$。因此，平面应力与平面应变的差别不仅是厚度，还包括弹性矩阵。

三维工程应变和应力取

$$
\boldsymbol\varepsilon_v=
[\varepsilon_{xx},\varepsilon_{yy},\varepsilon_{zz},
\gamma_{xy},\gamma_{yz},\gamma_{zx}]^{\mathrm T},
$$

$$
\boldsymbol\sigma_v=
[\sigma_{xx},\sigma_{yy},\sigma_{zz},
\sigma_{xy},\sigma_{yz},\sigma_{zx}]^{\mathrm T}.
$$

节点 $a$ 对应的三维 $\mathbf B_a$ 为

$$
\mathbf B_a=
\begin{bmatrix}
N_{a,x}&0&0\\
0&N_{a,y}&0\\
0&0&N_{a,z}\\
N_{a,y}&N_{a,x}&0\\
0&N_{a,z}&N_{a,y}\\
N_{a,z}&0&N_{a,x}
\end{bmatrix}.
\tag{5.185}
$$

各向同性三维弹性矩阵为

$$
\mathbf D=
\begin{bmatrix}
\lambda+2\mu&\lambda&\lambda&0&0&0\\
\lambda&\lambda+2\mu&\lambda&0&0&0\\
\lambda&\lambda&\lambda+2\mu&0&0&0\\
0&0&0&\mu&0&0\\
0&0&0&0&\mu&0\\
0&0&0&0&0&\mu
\end{bmatrix}.
\tag{5.186}
$$

若存在热应变、初始应变或本征应变 $\boldsymbol\varepsilon_0$，则应使用 $\boldsymbol\sigma=\mathbf D(\boldsymbol\varepsilon-\boldsymbol\varepsilon_0)$，而不是把热应力直接叠加到已经平衡的节点力上。各向异性材料只需替换 $\mathbf D$，后续面力离散形式不变。

## 5.15 二维节点控制体离散

### 5.15.1 对偶控制体的构造

围绕节点 $a$，依次连接相邻原始单元的边中点和单元中心，形成控制体 $\Omega_a$。在一个三角形单元内，节点控制体边界由“边中点—单元中心—相邻边中点”两段组成；在四边形单元内，也可采用同样的中位线构造。边界节点控制体还包含一部分真实结构边界。

记 $\mathcal E_a$ 为共享节点 $a$ 的单元集合，$\Omega_{a,e}=\Omega_a\cap\Omega_e$ 为单元 $e$ 对节点控制体的贡献，则

$$
\Omega_a=\bigcup_{e\in\mathcal E_a}\Omega_{a,e},
\qquad
\sum_{a\in\mathcal N_e}|\Omega_{a,e}|=|\Omega_e|.
\tag{5.187}
$$

第二个等式是控制体剖分必须满足的几何一致性。三角形中位控制体有 $|\Omega_{a,e}|=A_e/3$；矩形和平行四边形有 $|\Omega_{a,e}|=A_e/4$。一般畸变四边形应直接计算相应多边形面积，不能盲目使用 $A_e/4$。

### 5.15.2 惯性项与集中质量

在 $\Omega_a$ 上积分式（5.165）：

$$
\int_{\Omega_a}\rho\ddot{\mathbf u}\,\mathrm d\Omega
=\int_{\Omega_a}\nabla\cdot\boldsymbol\sigma\,\mathrm d\Omega
+\int_{\Omega_a}\mathbf b\,\mathrm d\Omega.
\tag{5.188}
$$

假设节点加速度代表整个节点控制体的平均加速度，得到集中质量近似

$$
\int_{\Omega_a}\rho\ddot{\mathbf u}\,\mathrm d\Omega
\approx m_a\ddot{\mathbf u}_a,
\qquad
m_a=\sum_{e\in\mathcal E_a}\rho_e|\Omega_{a,e}|.
\tag{5.189}
$$

二维时 $m_a$ 还包含厚度。该近似产生的质量矩阵为对角矩阵，存储和显式动力计算都很方便，并保证各单元质量按控制体体积分配。其代价是不能像一致质量矩阵那样完整表示单元内惯性分布；高频波动和振动问题应检查质量集中对频率及波速的影响。

### 5.15.3 应力散度的面力通量形式

利用高斯定理，应力散度转化为控制体边界上的面力：

$$
\int_{\Omega_a}\nabla\cdot\boldsymbol\sigma\,\mathrm d\Omega
=\int_{\partial\Omega_a}\boldsymbol\sigma\mathbf n\,\mathrm d\Gamma
=\sum_{e\in\mathcal E_a}
\int_{\Gamma_{a,e}}\boldsymbol\sigma_e\mathbf n\,\mathrm d\Gamma
+\int_{\Gamma_{a,t}}\overline{\mathbf t}\,\mathrm d\Gamma,
\tag{5.190}
$$

其中 $\Gamma_{a,e}$ 是位于单元 $e$ 内部的对偶控制体边界，$\Gamma_{a,t}$ 是真实力边界。内部节点没有第二项。面力矢量的分量为

$$
\begin{bmatrix}t_x\\t_y\end{bmatrix}
=
\begin{bmatrix}
\sigma_{xx}n_x+\sigma_{xy}n_y\\
\sigma_{xy}n_x+\sigma_{yy}n_y
\end{bmatrix}.
\tag{5.191}
$$

把 $\Gamma_{a,e}$ 分成直线段 $s$，可定义二维有向面积向量

$$
\mathbf A_s=A_s\mathbf n_s=(\Delta y_s,-\Delta x_s),
\tag{5.192}
$$

这里假设控制体边界按逆时针方向前进，内部在行进方向左侧，故右侧外法线为 $(\Delta y,-\Delta x)/A_s$。若程序采用其他绕向，必须同步改变法线公式。一个基本几何检查是：任一闭合控制体的有向面积向量之和应为零。

为了把面力和工程应力矢量联系起来，定义

$$
\mathsf N(\mathbf n)=
\begin{bmatrix}
n_x&0&n_y\\
0&n_y&n_x
\end{bmatrix},
\qquad
\mathbf t=\mathsf N(\mathbf n)\boldsymbol\sigma_v.
\tag{5.193}
$$

该矩阵将在刚度系数推导中反复使用。

### 5.15.4 由形函数导出控制体刚度系数

在单元内，应力由 $\boldsymbol\sigma_v=\mathbf D\mathbf B_e\mathbf d_e$ 计算。将内部控制体边界上的面力移到左端，定义正的内力刚度贡献：

$$
\mathbf f_{\mathrm{int},a}
=-\int_{\Gamma_{a,\mathrm{int}}}\mathbf t\,\mathrm d\Gamma
=\sum_{e\in\mathcal E_a}\sum_{b\in\mathcal N_e}
\mathbf K_{ab}^e\mathbf u_b.
\tag{5.194}
$$

于是每个单元对节点自由度 $a$、$b$ 的二维刚度块可写为

$$
\boxed{
\mathbf K_{ab}^e
=-\int_{\Gamma_{a,e}}
\mathsf N(\mathbf n)\mathbf D\mathbf B_b\,\mathrm d\Gamma
}
\tag{5.195}
$$

若采用数值积分，则

$$
\mathbf K_{ab}^e
\approx-\sum_{q\in\Gamma_{a,e}}
w_q\mathsf N(\mathbf n_q)\mathbf D_q\mathbf B_{b,q}.
\tag{5.196}
$$

式（5.195）就是课件中各 $k_{xx}$、$k_{xy}$、$k_{yx}$ 和 $k_{yy}$ 系数的紧凑矩阵形式。展开 $\mathbf D\mathbf B_b$ 后，即可得到由弹性常数、形函数偏导数和边界法向积分组成的标量公式。对于各向异性材料，只需使用完整的 $\mathbf D$，不必重新推导平衡方程。

对于 CST，$\mathbf B_b$ 和单元应力在单元内为常数，只需计算 $\int_{\Gamma_{a,e}}n_x\,\mathrm d\Gamma$ 和 $\int_{\Gamma_{a,e}}n_y\,\mathrm d\Gamma$。对双线性四边形，$\mathbf B_b$ 随位置变化，应沿“边中点—单元中心—边中点”各段选取积分点，并由式（5.179）计算局部梯度。积分点数量不足会降低系数精度，过多积分则不会弥补单元本身的插值误差。

### 5.15.5 节点平衡与整体矩阵

定义节点外载荷

$$
\mathbf F_a=
\int_{\Omega_a}\mathbf b\,\mathrm d\Omega
+\int_{\Gamma_{a,t}}\overline{\mathbf t}\,\mathrm d\Gamma.
\tag{5.197}
$$

节点离散平衡方程为

$$
m_a\ddot{\mathbf u}_a
+\sum_{e\in\mathcal E_a}\sum_{b\in\mathcal N_e}
\mathbf K_{ab}^e\mathbf u_b
=\mathbf F_a.
\tag{5.198}
$$

将所有节点的二维自由度按 $[u_{1x},u_{1y},u_{2x},u_{2y},\ldots]^{\mathrm T}$ 排列，得到

$$
\boxed{\mathbf M\ddot{\mathbf d}+\mathbf K\mathbf d=\mathbf F(t)}.
\tag{5.199}
$$

其中 $\mathbf M$ 为整体质量矩阵，$\mathbf K$ 为整体刚度矩阵，$\mathbf F$ 为等效节点载荷。刚度装配采用有限元中常见的“单元局部自由度到整体自由度”映射，但单元系数来自式（5.195）的控制体边界积分。

课件采用 $\mathbf M\ddot{\mathbf U}=\mathbf K_{\mathrm{slide}}\mathbf U+\mathbf F$ 的写法。与本章约定比较可得

$$
\mathbf K_{\mathrm{slide}}=-\mathbf K.
\tag{5.200}
$$

本章把内力项移到左端，因此在材料稳定且刚体位移被合理约束后，$\mathbf K$ 为对称正定型；未约束刚体模态时则为半正定型。若直接把课件右端矩阵当成通常的正刚度矩阵，会造成内力符号错误。

静力问题忽略加速度，得到

$$
\mathbf K\mathbf d=\mathbf F.
\tag{5.201}
$$

若体力和面力随时间缓慢变化，可对每个载荷时刻分别求解式（5.201），这称为准静态分析；它不包含惯性和阻尼，不能代替真正的动力计算。

## 5.16 三维空间离散与数值积分

### 5.16.1 三维控制体和质量分配

三维节点控制体由单元中心、面中心和边中点连接形成。四面体中，节点周围的控制体边界通常由若干三角形面片组成；六面体中则通常由四边形面片组成。三棱柱及其他多面体单元可按相同的几何剖分思想处理。

集中质量仍按式（5.189）计算，只需把 $|\Omega_{a,e}|$ 理解为体积。四面体中位剖分给每个节点贡献 $V_e/4$；规则六面体给每个节点贡献 $V_e/8$。对一般扭曲六面体或非标准剖分，应直接计算多面体体积，并检查所有节点贡献之和等于单元体积。

多面体可分解为若干四面体后计算体积。例如以参考点 $\mathbf x_0$ 和定向三角形面片 $(\mathbf x_1,\mathbf x_2,\mathbf x_3)$ 构造四面体，其有符号体积为

$$
V=\frac16(\mathbf x_1-\mathbf x_0)\cdot
[(\mathbf x_2-\mathbf x_0)\times(\mathbf x_3-\mathbf x_0)].
\tag{5.202}
$$

若采用有符号体积求和，所有面片必须具有一致的外法向；若仅计算凸控制体的正体积，也应避免不同面片方向混乱导致相互抵消。

### 5.16.2 三维面力通量与刚度块

三维高斯定理给出

$$
\int_{\Omega_a}\nabla\cdot\boldsymbol\sigma\,\mathrm dV
=\int_{\partial\Omega_a}\boldsymbol\sigma\mathbf n\,\mathrm dS.
\tag{5.203}
$$

采用式（5.185）的应力排列，三维面力可写为 $\mathbf t=\mathsf N_3(\mathbf n)\boldsymbol\sigma_v$，其中

$$
\mathsf N_3(\mathbf n)=
\begin{bmatrix}
n_x&0&0&n_y&0&n_z\\
0&n_y&0&n_x&n_z&0\\
0&0&n_z&0&n_y&n_x
\end{bmatrix}.
\tag{5.204}
$$

于是三维刚度块与二维形式完全相似：

$$
\boxed{
\mathbf K_{ab}^e
=-\int_{S_{a,e}}
\mathsf N_3(\mathbf n)\mathbf D\mathbf B_b\,\mathrm dS
}
\tag{5.205}
$$

这里 $S_{a,e}$ 是单元 $e$ 内属于节点 $a$ 控制体的内部边界面。该式展开后即可得到课件中沿 $x$、$y$、$z$ 三个方向的全部分量系数。矩阵形式能减少下标错误，并使二维和三维程序共享相同的装配逻辑。

### 5.16.3 面积向量与曲面积分

对每个平面多边形面片，可直接构造有向面积向量 $\mathbf A_f=A_f\mathbf n_f$。若面片顶点 $\mathbf x_i$ 按外法向满足右手规则的顺序排列，则

$$
\mathbf A_f
=\frac12\sum_{i=1}^{n_f}\mathbf x_i\times\mathbf x_{i+1},
\qquad
\mathbf x_{n_f+1}=\mathbf x_1.
\tag{5.206}
$$

该公式适用于顶点位于同一平面的多边形。对于非平面四边形，可将其三角化，或用两条对角线构造近似面积向量；必须在整个程序中采用一致的面片定义。所有闭合控制体仍应满足 $\sum_f\mathbf A_f=\mathbf 0$。

曲面积分可采用

$$
\int_{S_{a,e}}\phi\,\mathrm dS
\approx\sum_{q=1}^{n_q}w_q\phi_q,
\tag{5.207}
$$

其中 $w_q$ 已包含面片面积和参考域映射权重。四面体单元的 $\mathbf B_b$ 为常数，因此应力通量积分可退化为“常应力乘面积向量”；六面体和三棱柱的 $\mathbf B_b$ 通常随积分点变化，需要在每个面积分点计算雅可比矩阵及形函数梯度。

对于曲面六面体，单元面和对偶控制体面可能不是平面。此时应从参数映射计算切向量和法向量，而不能只把顶点坐标代入平面多边形公式。网格畸变会使雅可比、应力梯度和控制体几何同时恶化，网格质量检查是结构计算的一部分。

### 5.16.4 低阶单元的精度限制

三节点三角形和四节点四面体只能表示单元内常应变，往往需要较密网格才能准确描述弯曲和应力集中。双线性四边形和三线性六面体表达能力较强，但在某些弯曲问题中仍可能出现剪切锁定；采用减缩积分又可能引入零能沙漏模式。

因此，刚度矩阵具有刚体零空间只是必要条件，并不保证离散正确。还应检查：

- 是否存在除刚体模态之外的非物理零能模式；
- 常应变分片试验能否复现均匀应力；
- 厚板、弯曲和近不可压材料是否出现锁定；
- $\nu\rightarrow0.5$ 时平面应变或三维公式是否出现严重体积锁定；
- 网格加密后位移、反力和应力是否趋于稳定。

近不可压问题需要混合插值、选择性积分或其他专门格式，不能仅靠提高线性方程求解精度解决。

## 5.17 边界条件、载荷与整体装配

### 5.17.1 位移边界条件

指定位移属于 Dirichlet 条件。若节点 $a$ 的某自由度已知为 $\bar d_a$，最直接的处理是从未知量系统中消去该自由度。对任意自由行 $i$，原方程为

$$
\sum_jK_{ij}d_j=F_i.
$$

将已知位移移至右端：

$$
\sum_{j\in\mathcal F}K_{ij}d_j
=F_i-\sum_{j\in\mathcal D}K_{ij}\bar d_j,
\tag{5.208}
$$

其中 $\mathcal F$ 和 $\mathcal D$ 分别表示自由和受约束自由度。随后可直接求解缩减系统。与简单把对角元改成大数相比，消去法条件数更好，也便于准确计算反力。

对称边界通常要求法向位移为零、切向面力为零；反对称边界则应按具体载荷和位移场设置。三维无约束结构具有三个平移和三个转动刚体模态，二维平面结构具有两个平移和一个转动刚体模态。静力计算必须施加足够且不过度的约束，否则矩阵奇异或产生非物理附加反力。

### 5.17.2 面力、压力和体力

力边界不单独修改刚度矩阵，而是积分进入载荷向量：

$$
\mathbf F_{t,a}=\int_{\Gamma_{a,t}}\overline{\mathbf t}\,\mathrm d\Gamma.
\tag{5.209}
$$

若边界面承受流体压力 $p$，结构外法向为 $\mathbf n$，则通常有 $\overline{\mathbf t}=-p\mathbf n$。符号必须按“压力作用在结构表面并指向结构内部”检查，而不能只根据程序中法线数组的方向猜测。黏性面力、接触力或随变形变化的随动载荷也应按各自定义加入。

线性分布面力不宜在所有情况下都简单平均到两个节点。应沿真实边界段积分，并确保同一物理边界面只被计入一次。体力可按控制体积分；若 $\mathbf b_e$ 在单元内近似常数，则单元对节点 $a$ 的贡献为 $\mathbf b_e|\Omega_{a,e}|$。

热载荷可有两种等价实现方式：在本构关系中扣除热应变后形成热应力通量，或将其整理为等效节点热载荷。两种方式必须保持符号和边界处理一致，不能重复计入。

### 5.17.3 反力与整体平衡检查

得到位移后，未修正前的整体离散方程可用于计算约束反力。定义残差形式

$$
\mathbf R=\mathbf F-\mathbf M\ddot{\mathbf d}
-\mathbf C\dot{\mathbf d}-\mathbf K\mathbf d,
\tag{5.210}
$$

其中阻尼矩阵 $\mathbf C$ 在无阻尼时取零。受约束自由度上的外部反力与 $\mathbf R$ 的关系取决于其定义方向，后处理时应明确输出约定。静力状态下，所有外部载荷和约束反力之和应接近零，总力矩也应平衡。

整体平衡检查能够发现遗漏边界、错误载荷方向、重复装配和错误约束等问题。仅观察位移云图平滑并不能证明力平衡正确。

## 5.18 静力、瞬态动力与稳态响应求解

### 5.18.1 线性静力求解与收敛

线弹性静力问题只需在施加位移边界后求解式（5.201）。小规模或需要多次求解相同矩阵的问题可采用稀疏直接分解；大规模问题可采用共轭梯度等迭代法，并配合适当的预条件。合理约束后 $\mathbf K$ 对称正定，是使用共轭梯度法的重要依据。

若材料非线性、接触或几何非线性使内力依赖位移，则需求解

$$
\mathbf R(\mathbf d)=\mathbf f_{\mathrm{int}}(\mathbf d)-\mathbf F=0.
$$

Newton 迭代为

$$
\mathbf K_T^{(k)}\Delta\mathbf d^{(k)}
=-\mathbf R^{(k)},
\qquad
\mathbf d^{(k+1)}=\mathbf d^{(k)}+\alpha_k\Delta\mathbf d^{(k)},
\tag{5.211}
$$

其中 $\mathbf K_T$ 为切线刚度，$\alpha_k$ 可由线搜索确定。大变形问题还需更新构形、使用适当的应力度量并考虑几何刚度，不能继续把本章小变形矩阵当作常数重复使用。

### 5.18.2 阻尼与动力方程

实际结构常有材料、连接和周围介质引起的能量耗散。在线性模型中可加入速度项：

$$
\mathbf M\ddot{\mathbf d}
+\mathbf C\dot{\mathbf d}
+\mathbf K\mathbf d=\mathbf F(t).
\tag{5.212}
$$

一种常用模型是 Rayleigh 阻尼

$$
\mathbf C=\alpha_M\mathbf M+\beta_K\mathbf K,
\tag{5.213}
$$

其中 $\alpha_M$ 和 $\beta_K$ 分别控制质量比例和刚度比例阻尼。系数通常根据两个目标频率处的阻尼比确定。过高刚度比例阻尼会使高频模态衰减很快，并显著限制显式时间步，因此不能仅把它当作改善数值稳定性的随意参数。

初始加速度不一定为零。若初始位移和速度已给定，应由初始动力平衡求得

$$
\mathbf M\mathbf a_0
=\mathbf F(0)-\mathbf C\mathbf v_0-\mathbf K\mathbf d_0.
\tag{5.214}
$$

直接令 $\mathbf a_0=0$ 相当于额外施加了一个初始冲击，可能污染后续响应。

### 5.18.3 Newmark 时间积分

Newmark 方法假设新时间层位移和速度满足

$$
\mathbf d_{n+1}=\mathbf d_n+\Delta t\mathbf v_n
+\Delta t^2\left[
\left(\frac12-\beta\right)\mathbf a_n
+\beta\mathbf a_{n+1}\right],
\tag{5.215}
$$

$$
\mathbf v_{n+1}=\mathbf v_n+\Delta t
\left[(1-\gamma)\mathbf a_n+\gamma\mathbf a_{n+1}\right].
\tag{5.216}
$$

其中 $\mathbf v=\dot{\mathbf d}$，$\mathbf a=\ddot{\mathbf d}$。定义

$$
\begin{aligned}
&a_0=\frac1{\beta\Delta t^2},&&a_1=\frac\gamma{\beta\Delta t},\\
&a_2=\frac1{\beta\Delta t},&&a_3=\frac1{2\beta}-1,\\
&a_4=\frac\gamma\beta-1,&&a_5=\frac{\Delta t}{2}\left(\frac\gamma\beta-2\right).
\end{aligned}
\tag{5.217}
$$

由式（5.215）和式（5.216）可得

$$
\mathbf a_{n+1}
=a_0(\mathbf d_{n+1}-\mathbf d_n)-a_2\mathbf v_n-a_3\mathbf a_n,
\tag{5.218}
$$

$$
\mathbf v_{n+1}
=a_1(\mathbf d_{n+1}-\mathbf d_n)
-a_4\mathbf v_n-a_5\mathbf a_n.
\tag{5.219}
$$

代入 $t_{n+1}$ 时刻的动力平衡，整理为关于 $\mathbf d_{n+1}$ 的线性系统：

$$
\widehat{\mathbf K}\mathbf d_{n+1}=\widehat{\mathbf F}_{n+1},
\tag{5.220}
$$

其中

$$
\widehat{\mathbf K}=\mathbf K+a_0\mathbf M+a_1\mathbf C,
\tag{5.221}
$$

$$
\begin{aligned}
\widehat{\mathbf F}_{n+1}={}&\mathbf F_{n+1}
+\mathbf M(a_0\mathbf d_n+a_2\mathbf v_n+a_3\mathbf a_n)\\
&+\mathbf C(a_1\mathbf d_n+a_4\mathbf v_n+a_5\mathbf a_n).
\end{aligned}
\tag{5.222}
$$

求出位移后，再用式（5.218）和式（5.219）更新加速度与速度。线性问题中若 $\Delta t$、矩阵和边界约束不变，$\widehat{\mathbf K}$ 只需分解一次，可显著降低多时间步成本。

平均加速度法取 $\beta=1/4$、$\gamma=1/2$，在线性无阻系统中无条件稳定且二阶精度。更一般地，线性 Newmark 方法无条件稳定的常用参数条件为

$$
\gamma\geq\frac12,
\qquad
\beta\geq\frac14\left(\gamma+\frac12\right)^2.
\tag{5.223}
$$

“无条件稳定”不表示“任意时间步都准确”。时间步仍应足以解析载荷变化和所关注的最高结构频率。$\gamma>1/2$ 会引入算法阻尼，可用于衰减无法控制的高频离散模态，但也会改变物理响应。

一个时间步的实施顺序如下：

```text
已知 d_n、v_n、a_n 和外部载荷历史
由时间步和 beta、gamma 计算 a_0 至 a_5
形成 K_hat 和 F_hat，并施加新时刻的位移边界
求解 K_hat d_(n+1) = F_hat
由位移更新 a_(n+1) 和 v_(n+1)
计算单元应变、应力、节点反力及耦合界面量
检查动力平衡、能量和时间步敏感性
进入下一时间步
```

若位移边界随时间变化，不仅新时刻位移应满足约束，速度和加速度更新也应与约束时间导数一致。简单地在每步结束后覆盖位移而不处理速度，会引入虚假的动力误差。

### 5.18.4 稳态简谐响应与固有频率

若外载荷为 $\mathbf F(t)=\Re(\widehat{\mathbf F}e^{\mathrm i\omega t})$，稳态响应设为 $\mathbf d(t)=\Re(\widehat{\mathbf d}e^{\mathrm i\omega t})$，代入式（5.212）得

$$
\left(-\omega^2\mathbf M+\mathrm i\omega\mathbf C+\mathbf K\right)
\widehat{\mathbf d}=\widehat{\mathbf F}.
\tag{5.224}
$$

每个激励频率需求解一次复代数系统。无阻尼且 $\omega$ 接近固有频率时，响应会显著放大；有阻尼时幅值有限，但相位仍发生明显变化。

固有频率和模态由广义特征值问题确定：

$$
\mathbf K\boldsymbol\phi_r
=\omega_r^2\mathbf M\boldsymbol\phi_r.
\tag{5.225}
$$

模态可用于判断共振风险、选择时间步和构造降阶模型。集中质量与一致质量对高阶频率的影响通常比低阶频率更明显，因此验证时应同时比较若干阶频率，而不是只比较基频。

### 5.18.5 多物理场耦合中的结构更新

在流固耦合中，流体提供结构表面压力和黏性面力，结构向流体提供界面位移和速度。弱耦合可在每个时间步内只交换一次数据，强耦合则需反复执行“流体求解—结构求解—界面与网格更新”，直至界面残差收敛。

需要特别注意：

- 界面两侧面力应满足作用与反作用，位移和速度应采用一致的映射；
- 保守插值应保持总力，必要时还应检查总力矩和虚功；
- 结构集中质量、流体网格运动和时间积分必须使用同一时间层定义；
- 轻质结构与稠密流体耦合时可能出现显著附加质量效应，需要更强耦合或亚迭代；
- 热应力、接触和材料非线性应进入同一残差与收敛检查，而不能只检查位移变化。

## 5.19 完整求解流程与数值验证

### 5.19.1 结构求解流程

综合上述推导，格点型有限体积结构求解可按以下步骤实现：

1. **读取网格与材料**：确定节点、单元、边界、厚度、密度和弹性参数，检查单元方向与畸变。
2. **构造节点控制体**：连接边中点、面中心和单元中心，计算控制体体积、内部面片及其外法向。
3. **构建形函数**：按单元类型计算 $N_a$、参考坐标梯度、雅可比矩阵和整体坐标梯度。
4. **形成单元系数**：由 $\mathbf B$ 和 $\mathbf D$ 计算应力算子，按式（5.195）或式（5.205）积分刚度块，同时分配控制体质量。
5. **装配整体系统**：将局部自由度映射到整体自由度，累加 $\mathbf M$、$\mathbf K$、体力载荷和面力载荷。
6. **施加边界条件**：消去或约束指定位移自由度，保留用于反力恢复的原始矩阵和载荷。
7. **选择分析类型**：静力问题求解 $\mathbf K\mathbf d=\mathbf F$；动力问题采用 Newmark 等方法逐时间步求解；简谐问题求解频域方程。
8. **恢复物理量**：计算应变、应力、反力、应变能和界面运动，必要时进行应力平滑但保留原始单元值用于平衡检查。
9. **执行验证**：检查几何闭合、力平衡、网格与时间步收敛，以及解析解或实验对照。

非线性问题在步骤4至步骤7之间增加残差、切线刚度和外层迭代。耦合问题还需在外层交换其他物理场的数据。

### 5.19.2 分层验证策略

结构程序不宜直接用复杂工程算例作为第一次测试。建议按以下层级逐步验证。

**几何与插值检查。**

- 每个单元的控制体子区域体积之和等于单元体积；
- 闭合控制体的有向面积向量之和为零；
- $\sum_aN_a=1$，$\sum_a\nabla N_a=0$；
- 随机线性位移场在单元内能被精确插值。

**单元与分片试验。**

- 刚体平移和小转动不产生应变与应力；
- 常应变分片试验复现解析应力和节点平衡；
- 单轴拉伸复现 $E$、$\nu$ 对应的轴向及横向变形；
- 纯弯、纯剪算例用于识别锁定和异常零能模式。

**整体平衡与收敛。**

- 外载荷与约束反力满足合力和合力矩平衡；
- 对称问题采用全模型与半模型得到一致结果；
- 至少使用三套逐渐加密的网格检查位移、反力和应力；
- 动力问题同时检查时间步收敛，简谐问题检查频率间隔是否足以解析共振峰。

可定义归一化离散残差和位移增量指标

$$
\eta_R=\frac{\|\mathbf R\|_2}{F_{\mathrm{ref}}},
\qquad
\eta_d=\frac{\|\Delta\mathbf d\|_2}{d_{\mathrm{ref}}},
\tag{5.226}
$$

其中 $F_{\mathrm{ref}}$ 和 $d_{\mathrm{ref}}$ 是与问题尺度相适应的参考值。只使用绝对残差会在载荷尺度改变时失去可比性。应力是位移梯度量，通常比位移收敛慢，在单元边界还可能不连续，因此应分别设置位移和应力的网格收敛判断。

### 5.19.3 常见实现错误

以下错误在结构有限体积实现中尤其常见：

1. 三角形面积取绝对值，却没有统一节点方向，导致部分单元梯度符号错误；
2. 把 $\nabla_{\mathbf x}N$ 错写为 $\mathbf J^{-1}\nabla_{\boldsymbol\xi}N$，遗漏式（5.179）中的转置约定；
3. 将工程剪应变与张量剪应变混用，使剪切刚度相差两倍；
4. 平面应力问题误用平面应变弹性矩阵；
5. 控制体面法线未统一朝外，或只积分部分面片；
6. 把课件右端刚度矩阵直接当作正刚度，遗漏式（5.200）的符号转换；
7. 面力在相邻边界上重复积分，或压力法向采用流体侧而非结构侧方向；
8. 静力分析未消除全部刚体模态，却用迭代求解器的不收敛误判为网格问题；
9. 动力分析随意令初始加速度为零，或改变时间步后仍沿用旧的有效矩阵；
10. 只看位移云图，不检查反力、能量、网格和时间步敏感性。

## 本章小结

结构场格点型有限体积离散包含两套相互配合的几何与插值结构：原始单元用于定义形函数、材料参数和应力，节点对偶控制体用于建立惯性、面力和外载荷的积分平衡。形函数满足节点插值和单位分解条件，整体坐标梯度通过雅可比矩阵由参考坐标梯度求得，应变则由 $\mathbf B$ 矩阵组装。

对节点控制体应用高斯定理后，应力散度转化为控制体边界面力通量。二维和三维可统一写成 $-\int\mathsf N\mathbf D\mathbf B\,\mathrm d\Gamma$ 的刚度块形式，区别仅在工程应力维数、面力矩阵和边界测度。集中质量由各单元贡献的控制体体积计算。课件把刚度项保留在右端，本章则采用 $\mathbf M\ddot{\mathbf d}+\mathbf K\mathbf d=\mathbf F$ 的正刚度约定，两者相差一个负号。

线弹性静力问题最终求解 $\mathbf K\mathbf d=\mathbf F$；动力问题还需处理阻尼、初始加速度和时间积分。Newmark 方法可整理为以新时刻位移为未知量的有效刚度系统，稳态简谐响应则对应复频域方程。可靠的结构计算还必须通过几何闭合、刚体与常应变分片试验、反力平衡、网格收敛和时间步收敛检查，不能只以位移云图作为正确性依据。

## 习题

1. 从牛顿第二定律推导式（5.165），说明静力、准静态和稳态简谐响应三种分析的区别。
2. 证明形函数的单位分解性是精确表示刚体平移的必要条件，并推导 $\sum_a\nabla N_a=0$。
3. 从 $u=a_1+a_2x+a_3y$ 出发，完整推导三节点三角形形函数及其梯度。
4. 对四节点四边形推导式（5.174），并计算参考单元中心和四个节点处的形函数梯度。
5. 由链式法则推导式（5.179）。说明为什么畸变单元中不能直接使用参考坐标梯度代替整体坐标梯度。
6. 写出二维平面应力和平面应变的 $\mathbf D$，并验证式（5.183）的等效 Lamé 参数。
7. 对任意二维节点控制体证明闭合边界满足 $\sum_s\mathbf A_s=\mathbf 0$，讨论该恒等式如何用于程序调试。
8. 将式（5.195）展开为二维 $2\times2$ 刚度块，分别写出 $k_{xx}$、$k_{xy}$、$k_{yx}$ 和 $k_{yy}$。
9. 对常应变三角形说明为什么应力通量积分只需边界法向积分；对双线性四边形为何还需要面积分点？
10. 推导式（5.198），并将课件写法 $\mathbf M\ddot{\mathbf U}=\mathbf K_{\mathrm{slide}}\mathbf U+\mathbf F$ 转换为本章的正刚度形式。
11. 说明格点型有限体积刚度块与有限元 $\int\mathbf B^{\mathrm T}\mathbf D\mathbf B\,\mathrm d\Omega$ 的方程来源有何不同。
12. 推导三维面力矩阵式（5.204），验证其每一行分别对应 $t_x$、$t_y$ 和 $t_z$。
13. 说明四面体中位控制体为何向每个节点分配 $V_e/4$，并讨论扭曲六面体为什么不能机械套用 $V_e/8$。
14. 设计一个常应变分片试验，列出网格、材料、边界位移、解析应力及反力检查项目。
15. 为什么未约束刚体模态时静力刚度矩阵奇异？分别给出二维和三维合理约束的例子。
16. 推导 Newmark 有效刚度式（5.221）和有效载荷式（5.222）。
17. 证明平均加速度法的参数 $\beta=1/4$、$\gamma=1/2$ 满足式（5.223），并解释无条件稳定与时间精度的区别。
18. 从式（5.212）推导稳态简谐方程式（5.224），说明共振时阻尼的作用。
19. 在流固耦合中，为什么只交换力而不检查界面位移、速度和时间层一致性仍可能产生错误？
20. 编写一份结构求解器验证清单，至少包含几何、插值、分片试验、整体平衡、网格收敛和时间步收敛六类检查。

## 延伸阅读

1. Zienkiewicz, O. C., Taylor, R. L., and Zhu, J. Z. *The Finite Element Method: Its Basis and Fundamentals*. 7th ed., Butterworth-Heinemann, 2013.
2. Bathe, K. J. *Finite Element Procedures*. 2nd ed., Prentice Hall, 1996.
3. Hughes, T. J. R. *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*. Dover, 2000.
4. Belytschko, T., Liu, W. K., Moran, B., and Elkhodary, K. *Nonlinear Finite Elements for Continua and Structures*. 2nd ed., Wiley, 2014.
5. Newmark, N. M. A method of computation for structural dynamics. *Journal of the Engineering Mechanics Division, ASCE*, 1959, 85(3): 67–94.
6. Eymard, R., Gallouët, T., and Herbin, R. Finite volume methods. In *Handbook of Numerical Analysis*, Vol. 7, North-Holland, 2000: 713–1020.
