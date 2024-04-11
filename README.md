<div align="center">

# 一维有限元算例(README还需要修改, 请勿阅读)

</div>

考虑如下的特征值问题, $k$ 为参数, $\lambda$ 为特征值:

$$
\begin{cases}
-u^{''}(x)+cos(kx)u(x)={\lambda}u(x),x\in[-1,1]\\
u(-1)=u(1)=0
\end{cases}
$$

记 $q(x)=cos(kx)$ , 则原问题弱形式为: 求 $u{\in}H^1_0(-1,1)$ 满足

$$
\int^1_{-1}u^{'}(x)v^{'}(x)dx+\int^1_{-1}q(x)u(x)v(x)dx=\lambda\int^1_{-1}u(x)v(x)dx,\,\,{\forall}v(x){\in}H^1_0(-1,1)
$$

将 $[-1,1]$ 均匀划分为 $N$ 个子区间, 子区间长为 $h=\frac{2}{N}$ , 区间端点分别为 $x_0,x_1,...,x_N$ 且 $x_i=-1+ih$ , 构造对应于 $x_i$ 的节点基函数:


$$
\phi_i(x)=
\begin{cases}
\frac{1}{h}(x-x_{i-1}),x\in[x_{i-1},x_i)\\
-\frac{1}{h}(x-x_{i+1}),x\in[x_i,x_{i+1}]\\
0,其它\\
\end{cases}
$$


故可求节点基函数的导数为：

$$
\phi_i^{'}(x)=
\begin{cases}
\frac{1}{h},x\in(x_{i-1},x_i)\\
-\frac{1}{h},x\in(x_i,x_{i+1})\\
0,x\in[-1,x_{i-1})\bigcup(x_{i+1},1]\\
不存在,其它\\
\end{cases}
$$

据此可求 $\phi_i(x)\phi_j(x)$ 与 $\phi_i^{'}(x)\phi_j^{'}(x)$ : $(i,j=1,2,...,N-1)$ <br/>
<br/>
(1) 当 $i=j$ 时:

$$
\phi_i(x)\phi_i(x)=
\begin{cases}
\frac{1}{h^2}(x-x_{i-1})^2,x\in[x_{i-1},x_i)\\
\frac{1}{h^2}(x-x_{i+1})^2,x\in[x_i,x_{i+1}]\\
0,其它\\
\end{cases}
$$

及

$$
\phi_i^{'}(x)\phi_i^{'}(x)=
\begin{cases}
\frac{1}{h^2},x\in(x_{i-1},x_i)\bigcup(x_i,x_{i+1})\\
0,x\in[-1,x_{i-1})\bigcup(x_{i+1},1]\\
不存在,其它\\
\end{cases}
$$

(2) 当 $j=i+1$ 时( $i=j+1$ 的情形类似 ):

$$
\phi_i(x)\phi_{i+1}(x)=
\begin{cases}
-\frac{1}{h^2}(x-x_{i+1})(x-x_i),x\in[x_i,x_{i+1}]\\
0,其它\\
\end{cases}
$$

及

$$
\phi_i^{'}(x)\phi_{i+1}^{'}(x)=
\begin{cases}
-\frac{1}{h^2},x\in(x_i,x_{i+1})\\
0,x\in[-1,x_{i-1})\bigcup(x_{i-1},x_i)\bigcup(x_{i+1},1]\\
不存在,其它\\
\end{cases}
$$

(3) 当 $|i-j|>1$ 时:

$$
\phi_i(x)\phi_j(x)=0,x\in[-1,1]\\
$$

及

$$
\phi_i^{'}(x)\phi_j^{'}(x)=
\begin{cases}
0,x\in[-1,x_{i-1})\bigcup(x_{i-1},x_i)\bigcup(x_i,x_{i+1})\bigcup(x_{i+1},1]\\
不存在,其它\\
\end{cases}
$$

为方便讨论, 将基函数导数无定义的点处的基函数导数值设为 $0$

弱问题的有限元估计为: 求 $u_h(x){\in}V_h=span[\phi_1,\phi_2,...,\phi_{N-1}]$ 满足:

$$
\int^1_{-1}u_h^{'}(x)v_h^{'}(x)dx+\int^1_{-1}q(x)u_h(x)v_h(x)dx=\lambda\int^1_{-1}u_h(x)v_h(x)dx,{\forall}v_h{\in}V_h
$$

设 $u_h(x)={\sum}_{i=1}^{N-1}U_i\phi_i(x)$ , 则上述问题可写为: 求 $U_i{\in}\mathbb{R}$ 满足: 

$$
\sum_{i=1}^{N-1}{U_i\int^1_{-1}[\phi_i^{'}(x)\phi_j^{'}(x)+q(x)\phi_i(x)\phi_j(x)]dx}=\sum_{i=1}^{N-1}{{\lambda}U_i\int^1_{-1}\phi_i(x)\phi_j(x)dx,{\forall}j\in\{1,2,...,N-1\}}
$$

记

$$
a_{ji}=\int^1_{-1}[\phi_i^{'}(x)\phi_j^{'}(x)+q(x)\phi_i(x)\phi_j(x)]dx
$$

及

$$
b_{ji}=\int^1_{-1}\phi_i(x)\phi_j(x)dx
$$

设 $A=(a_{ji}){\in}\mathbb{R}^{(N-1)\times(N-1)}$ , $B=(b_{ji}){\in}\mathbb{R}^{(N-1)\times(N-1)}$ , $U=(U_{i}){\in}\mathbb{R}^{N-1}$, 则易知 $A$ 与 $B$ 均为对称三对角矩阵, 且问题转化为特征值问题:

$$
AU={\lambda}BU
$$

下面求解 $A$ 与 $B$ 的主对角元与次对角元:

(1) $A$ 的主对角元, 对应代码中 **a_diag** :

$$
a_{ii}=\frac{1}{h}\left[2+\frac{1}{h}\left[\int_{-1+(i-1)h}^{-1+ih}q(x)(x+1-(i-1)h)^2dx+\int_{-1+ih}^{-1+(i+1)h}q(x)(x+1-(i+1)h)^2dx\right]\right]
$$

(2) $A$ 的上次对角元, 对应代码中 **a_subdiag**:

$$
a_{i,i+1}=-\frac{1}{h}\left[1+\frac{1}{h}\int^{-1+(i+1)h}_{-1+ih}q(x)(x+1-(i+1)h)(x+1-ih)dx\right]
$$

(3) $B$ 的主对角元, 对应代码中 **b_diag**:

$$
b_{ii}=\frac{1}{h^2}\left[\int_{-1+(i-1)h}^{-1+ih}(x+1-(i-1)h)^2dx+\int_{-1+ih}^{-1+(i+1)h}(x+1-(i+1)h)^2dx\right]
$$

(4) $B$ 的上次对角元, 对应代码中 **b_subdiag**:

$$
b_{i,i+1}=-\frac{1}{h^2}\int_{-1+ih}^{-1+(i+1)h}(x+1-(i+1)h)(x+1-ih)dx
$$

本说明更新时间为2024年3月23日 <br/>

