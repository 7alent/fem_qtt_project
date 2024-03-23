% 基于有限元方法求解特征值问题-u''+cos(kx)*u=lambda*u
% u(x)定义域为(-1,1), 边值条件为u(-1)=u(1)=0

%--------------------------------------------------------------------------
% 固定k, 比较不同N对求解结果的影响

% 初始化结果数组
n = 1000;
k = 1;
m = 5;
N_vec = linspace(100, 200, m);
lambda_vec = 1:m;
time_vec = zeros(1, m);
u_mat = zeros(n, m);

% 求解
for ii = 1:m
    fprintf(sprintf("正在进行 k = %d , N = %d 的求解...", [k N_vec(ii)]))
    tic;
    [~, u_mat(:, ii), lambda_vec(ii)] = my_fem(k, N_vec(ii), 0, n);
    time_vec(ii) = time_vec(ii)+toc;
    fprintf(sprintf("已完成求解, 耗时 %.6f 秒\n", time_vec(ii)));
end

% 绘图
hold on
legend_cell = cell(m, 1);
for ii = 1:m
    legend_cell(ii, 1) = {sprintf("k = %d, N = %d, t = %d, λ = %.6g", ...
        [k, N_vec(ii), round(time_vec(ii)), lambda_vec(ii)])};
    plot(linspace(-1, 1, n), u_mat(:, ii), LineWidth = 1.5);
end
legend(legend_cell, "FontSize", 12);
title(sprintf("Smallest Eigenvalue Solution when k = %d", k), ...
    "FontSize", 20);
xlabel("x","FontSize", 20);
ylabel("u(x)","FontSize", 20);
hold off

%--------------------------------------------------------------------------
% 固定N, 比较不同k对求解结果的影响

% 初始化结果数组
n = 1000;
N = 100;
m = 5;
k_vec = [1 10 100 1000 10000];
lambda_vec = 1:m;
time_vec = zeros(1, m);
u_mat = zeros(n, m);

% 求解
for ii = 1:m
    fprintf(sprintf("正在进行 k = %d , N = %d 的求解...", [k_vec(ii) N]))
    tic;
    [~, u_mat(:, ii), lambda_vec(ii)] = my_fem(k_vec(ii), N, sigma, 0, n);
    time_vec(ii) = time_vec(ii)+toc;
    fprintf(sprintf("已完成求解, 耗时 %.6f 秒\n", time_vec(ii)));
end

% 绘图
hold on
legend_cell = cell(m, 1);
for ii = 1:m
    legend_cell(ii, 1) = {sprintf("k = %d, N = %d, t = %d, λ = %.6g", ...
        [k_vec(ii), N, round(time_vec(ii)), lambda_vec(ii)])};
    plot(linspace(-1, 1, n), u_mat(:, ii), LineWidth = 1.5);
end
legend(legend_cell, "FontSize", 12);
title(sprintf("Smallest Eigenvalue Solution when N = %d", N), ...
    "FontSize", 20);
xlabel("x","FontSize", 20);
ylabel("u(x)","FontSize", 20);
hold off
