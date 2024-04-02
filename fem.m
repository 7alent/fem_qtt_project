% 基于有限元方法求解一维二阶椭圆型微分方程特征值问题(狄利克雷边界)
% -u''(x)+q(x)u(x)=lambda*u(x), a<=x<=b
% u(a)=0; u(b)=0
% 求(按模)最小特征值对应的解函数
% x 为自变量(符号变量)
% q 为系数函数(关于 x 的单变量函数句柄)
% N 为元个数(整数)
% a, b 为求解域边界点(标量)
% (暂不支持的参数: P, T, P_b, T_b, 积分仅涉及一阶导数)
% (暂不支持采用位于不同空间的试探与测试函数, 不支持其它边界条件)

function [U, lambda, L, R, A, B, Ud1, Ud2, qU, V, D, Uerr] = fem(x, q, N, a, b, normalize_by)
    
    % 初始化变量及信息矩阵
    h = (b-a)/N;
    P = linspace(a, b, N+1); % 第n列为第n个节点的坐标
    N_m = size(P, 2); % 节点总数
    T = [1:N_m-1; 2:N_m]; % 第n列为第n个单元上所有节点的全局编号
    P_b = P; % 第n列为第n个有限元基函数的节点坐标
    T_b_trial = T; % 第n列为第n个单元上所有试探基函数节点的全局编号
    T_b_test = T; % 第n列为第n个单元上所有试探基函数节点的全局编号
    N_lb_trial = size(T_b_trial, 1); % 局部节点试探基函数个数
    N_lb_test = size(T_b_test, 1); % 局部节点测试基函数个数
    A_1 = sparse(N_m, N_m); % 等号左侧总刚度矩阵(拉普拉斯阵)
    A_2 = sparse(N_m, N_m); % 等号左侧总刚度矩阵(势能项)
    B = sparse(N_m, N_m); % 等号右侧总刚度矩阵
    
    % 合成刚度矩阵
    for n = 1:N
        nf = cell2struct(cell(4, 1), ...
            {'f1', 'f1d', 'f2', 'f2d'}, 1); % 局部节点基函数及其导数
        if n>1
            nf.f1 = (P(n+1)-x)/h;
            nf.f1d = -1/h;
        else % 第一个单元
            nf.f1 = 0*x;
            nf.f1d = 0*x;
        end
        if n<N
            nf.f2 = (x-P(n))/h;
            nf.f2d = 1/h;
        else % 最后一个单元
            nf.f2 = 0*x;
            nf.f2d = 0*x;
        end
        S_A_1 = zeros(N_lb_trial, N_lb_test); % 等号左侧单元刚度矩阵(拉普拉斯阵)
        S_A_2 = zeros(N_lb_trial, N_lb_test); % 等号左侧单元刚度矩阵(势能项)
        S_B = zeros(N_lb_trial, N_lb_test); % 等号右侧单元刚度矩阵
        for alpha = 1:N_lb_trial
            for beta = 1:N_lb_test % 局部节点测试基函数个数
                phi = nf.(['f', num2str(alpha)]); % 局部试探基函数
                phi_d = nf.(['f', num2str(alpha), 'd']); % 局部试探基函数导数
                psi = nf.(['f', num2str(beta)]); % 局部测试基函数
                psi_d = nf.(['f', num2str(beta), 'd']); % 局部测试基函数导数
                integrand_A_1 = phi_d*psi_d;
                integrand_A_2 = q*phi*psi;
                integrand_B = phi*psi;
                S_A_1(alpha, beta) = gauss_quad(x, integrand_A_1, ...
                    P(T(1, n)), P(T(2, n)));
                S_A_2(alpha, beta) = gauss_quad(x, integrand_A_2, ...
                    P(T(1, n)), P(T(2, n)));
                S_B(alpha, beta) = gauss_quad(x, integrand_B, ...
                    P(T(1, n)), P(T(2, n)));
            end
        end
        % 单元刚度矩阵合成总刚度矩阵
        A_1(T_b_trial(:, n), T_b_test(:, n)) = ...
            A_1(T_b_trial(:, n), T_b_test(:, n))+S_A_1;
        A_2(T_b_trial(:, n), T_b_test(:, n)) = ...
            A_2(T_b_trial(:, n), T_b_test(:, n))+S_A_2;
        B(T_b_trial(:, n), T_b_test(:, n)) = ...
            B(T_b_trial(:, n), T_b_test(:, n))+S_B;
    end
    A_1 = A_1(2:N_m-1, 2:N_m-1);
    A_2 = A_2(2:N_m-1, 2:N_m-1);
    A = A_1+A_2;
    B = B(2:N_m-1, 2:N_m-1);
    
    % 求解特征值问题
    %[U, lambda] = eigs(A, B, 1, 'smallestreal');
    [V, D] = eig(full(B)\full(A));
    [lambda, I] = min(diag(D));
    U = V(:, I);
    
    % 确定特征向量符号
    % -U与U均可作为特征向量
    % 为使得绘图后曲线整体位于x轴上方, 调整U的符号使得U具有更多非负元素
    if sum(U >= 0) < size(U, 1)/2
		U = -U;
    end

    % 计算部分输出变量
    L = A*U;
    R = lambda*B*U;
    Ud2 = A_1*U; % u的二阶导数
    Ud2 = [0; Ud2; 0];
    qU = A_2*U; % u与q的乘积
    qU = [0; qU; 0];
    U = [0; U; 0];
    Ud1 = U; % (近似)一阶导数
    Ud1(2:N_m-1) = (Ud1(3:N_m)-Ud1(1:N_m-2))/(2*h); % 中心差分
    Ud1(1) = Ud1(2);
    Ud1(N_m) = Ud1(N_m-1);

    % 归一化返回解函数
    U = U/fem_norm(U, 0.*U, normalize_by);
    Uerr = U-cos(pi*linspace(a, b, N+1)/2)'; % 解函数与cos(pi*x/2)的绝对误差

    % 绘制图象
    plot(P, U);

end
