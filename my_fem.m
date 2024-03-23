function [solution, u, lambda] = my_fem(k, N, draw_plot, n)
% 基于有限元方法求解特征值问题-u''+cos(kx)*u=lambda*u的最小特征值解
% u(x)定义域为(-1,1), 边值条件为u(-1)=u(1)=0
% k为参数, lambda为特征值, N为划分区间个数
% 设g_i为对应分割点x_i的分段线性基函数
% 设u=U_1*g_1+...+U_(N-1)*g_(N-1)
% 记a_ji=integral(g_i'*g_j'+cos(kx)*g_i*g_j, a, b)
% 记b_ji=integral(g_i*g_j, a, b)
% 记U=(U_1,...,U_(N-1))'
% 记A=(a_ji), B=(b_ji)为(N-1)*(N-1)的对称三对角矩阵
% 则问题化为广义特征值问题A*U=lambda*B*U
% 返回值解释:
% solution 为解函数
% u 为解函数在区间上的 n 个等分点上的值, 用于绘图
% 部分参数解释:
% draw_plot 为1将返回解函数图象(基于fplot函数)
% n 与 u 的计算有关

    syms x ii
    a = -1;
    b = 1;
    h = (b-a)/N; % 区间长
    q = cos(k*x);

    % 构建矩阵A
    a_diag = 1:N-1; % A对角元a_i_i
    a_diag_formula = ((int(q*(x+1-(ii-1)*h)^2, x, -1+(ii-1)*h, -1+ii*h)+ ...
        int(q*(x+1-(ii+1)*h)^2, x, -1+ii*h, -1+(ii+1)*h))/h+2)/h;
    for jj = 1:N-1
        a_diag(jj) = subs(a_diag_formula, ii, jj);
    end
    a_subdiag = 1:N-2; % A下次对角元a_{i+1}_i = a_i_{i+1}
    a_subdiag_formula = -(int(q*(x+1-(ii+1)*h)*(x+1-ii*h), ...
        -1+ii*h, -1+(ii+1)*h)/h+1)/h;
    for jj = 1:N-2
        a_subdiag(jj) = subs(a_subdiag_formula, ii, jj);
    end
    A = sparse([1:N-1 1:N-2 2:N-1], [1:N-1 2:N-1 1:N-2], ...
        [a_diag a_subdiag a_subdiag], N-1, N-1);
    
    % 构建矩阵B
    b_diag = 1:N-1; % B对角元
    b_diag_formula = (int((x+1-(ii-1)*h)^2, -1+(ii-1)*h, -1+ii*h)+ ...
        int((x+1-(ii+1)*h)^2, -1+ii*h, -1+(ii+1)*h))/h^2;
    for jj = 1:N-1
        b_diag(jj) = subs(b_diag_formula, ii, jj);
    end
    b_subdiag = 1:N-2; % A下次对角元a_i+1_i = a_i_i+1
    b_subdiag_formula = -int((x+1-(ii+1)*h)*(x+1-ii*h), ...
        -1+ii*h, -1+(ii+1)*h)/h^2;
    for jj = 1:N-2
        b_subdiag(jj) = subs(b_subdiag_formula, ii, jj);
    end
    B = sparse([1:N-1 1:N-2 2:N-1], [1:N-1 2:N-1 1:N-2], ...
        [b_diag b_subdiag b_subdiag], N-1, N-1);

    % 计算广义特征对
    [U, lambda] = eigs(A, B, 1, "smallestreal");
    if norm(U) ~= 1
        U = U/norm(U);
    end

    % 返回有限元解
    solution = 0; % 初始化解函数
    for jj = 1:N-1
        solution = solution+U(jj)*piecewise(x <= -1+(jj-1)*h, ...
            0, -1+(jj-1)*h < x & x < -1+jj*h, (x+1-(jj-1)*h)/h, ...
            -1+jj*h <= x & x <= -1+(jj+1)*h, -(x+1-(jj+1)*h)/h, ...
            -1+(jj+1)*h < x, 0);
    end

    % 绘制图象
    if draw_plot
        fplot(solution, [-1 1], MeshDensity = n)
    end
    
    % 将区间n等分并计算解函数在n等分点上的值
    u = double(subs(solution, x, linspace(-1, 1, n)));
    for jj = 1:n
        if u(jj) < 0
            u = -u;
            solution = -solution;
            break
        end
    end

    % 节点处误差分析(待写)
end