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

function [U, lambda] = fem(x, q, N, a, b, normalize_by)
    
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
    A = sparse(N_m, N_m); % 等号左侧总刚度矩阵
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
        S_A = zeros(N_lb_test, N_lb_trial); % 等号左侧单元刚度矩阵
        S_B = zeros(N_lb_test, N_lb_trial); % 等号右侧单元刚度矩阵
        for alpha = 1:N_lb_trial
            for beta = 1:N_lb_test % 局部节点测试基函数个数
                phi = nf.(['f', num2str(alpha)]); % 局部试探基函数
                phi_d = nf.(['f', num2str(alpha), 'd']); % 局部试探基函数导数
                psi = nf.(['f', num2str(beta)]); % 局部测试基函数
                psi_d = nf.(['f', num2str(beta), 'd']); % 局部测试基函数导数
                integrand_A = phi_d*psi_d+q*phi*psi;
                integrand_B = phi*psi;
                if size(char(integrand_A), 2) > 1
                    S_A(beta, alpha) = quadgk(matlabFunction(integrand_A), ...
                        P(T(1, n)), P(T(2, n)));
                end
                if size(char(integrand_B), 2) > 1
                    S_B(beta, alpha) = quadgk(matlabFunction(integrand_B), ...
                        P(T(1, n)), P(T(2, n)));
                end
            end
        end
        % 单元刚度矩阵合成总刚度矩阵
        A(T_b_test(:, n), T_b_trial(:, n)) = ...
            A(T_b_test(:, n), T_b_trial(:, n))+S_A;
        B(T_b_test(:, n), T_b_trial(:, n)) = ...
            B(T_b_test(:, n), T_b_trial(:, n))+S_B;
    end
    A = A(2:N_m-1, 2:N_m-1);
    B = B(2:N_m-1, 2:N_m-1);

    % 求解特征值问题
    [U, lambda] = eigs(A, B, 1, 'smallestabs');
    U = [0; U; 0];

    % 确定特征向量符号
    % -U与U均可作为特征向量
    % 为使得绘图后曲线整体位于x轴上方, 调整U的符号使得U具有更多非负元素
    if sum(U >= 0) < size(U, 1)/2
		U = -U;
    end

    % 归一化返回解函数
    U = U/fem_norm(U, 0.*U, normalize_by);

    % 绘制图象
    plot(P, U);

end