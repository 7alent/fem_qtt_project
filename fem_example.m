% 基于有限元方法求解一维二阶椭圆型微分方程特征值问题(狄利克雷边界)
% -u''(x)+p(x)u'(x)+q(x)u(x)+r(x)=lambda*u(x), a<=x<=b
% u(a)=0; u(b)=0
% 求(按模)最小特征值对应的解函数

% ------------------探究 U_err 与 k 的关系(k = 2^4 ~ 2^8)-------------------
% 参数初始化
syms x
N = 2^11;
a= -1;
b = 1;
normalize_by = "L2";

% 计算(固定 N 调整 k)
m = 5;
U_mat_1 = zeros(N+1, m);
Uerr_mat_1 = zeros(N+1, m);
Ud1_mat_1 = zeros(N+1, m);
Ud2_mat_1 = zeros(N+1, m);
qU_mat_1 = zeros(N+1, m);
lambda_vec_1 = 1:m;
for k = 1:m
    [U_mat_1(:, k), lambda_vec_1(k), ~, ~, ~, ~, Ud1_mat_1(:, k), ...
        Ud2_mat_1(:, k), qU_mat_1(:, k), ~, ~, Uerr_mat_1(:, k)] = ...
        fem(x, cos(2^(k+3)*x), N, a, b, normalize_by);
end

% 解函数与 cos(pi*x/2) 的图象
for k = 1:5
    plot(linspace(a, b, N+1), U_mat_1(:, k), "LineWidth", 1.5)
    hold on
end
plot(linspace(a, b, N+1), cos(pi*linspace(a, b, N+1)/2), "LineWidth", 1.5)
set(gca,"FontSize", 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of $cos(\frac{{\pi}}{2}x)$ and $u(x)$ when $k = 2^4\sim2^8, N = 2^{11}$", "FontSize", 24, "Interpreter", "latex")
legend("$k = 16$", "$k = 32$", "$k = 64$", "$k = 128$", "$k = 256$", "$cos(\frac{{\pi}}{2}x)$", "FontSize", 20, "Interpreter", "latex")
hold off

% 比较不同 k 下 U_err 的波动性
plot(linspace(a, b, N+1), Uerr_mat_1(:, 1), "LineWidth", 1.5)
set(gca,'FontSize', 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of $u(x)-cos(\frac{{\pi}}{2}x)$ when $k = 16, N = 2^{11}$", "FontSize", 24, "Interpreter", "latex")

plot(linspace(a, b, N+1), Uerr_mat_1(:, 2), "LineWidth", 1.5)
set(gca,'FontSize', 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of $u(x)-cos(\frac{{\pi}}{2}x)$ when $k = 32, N = 2^{11}$", "FontSize", 24, "Interpreter", "latex")

plot(linspace(a, b, N+1), Uerr_mat_1(:, 3), "LineWidth", 1.5)
set(gca,'FontSize', 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of $u(x)-cos(\frac{{\pi}}{2}x)$ when $k = 64, N = 2^{11}$", "FontSize", 24, "Interpreter", "latex")

plot(linspace(a, b, N+1), Uerr_mat_1(:, 4), "LineWidth", 1.5)
set(gca,'FontSize', 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of $u(x)-cos(\frac{{\pi}}{2}x)$ when $k = 128, N = 2^{11}$", "FontSize", 24, "Interpreter", "latex")

plot(linspace(a, b, N+1), Uerr_mat_1(:, 5), "LineWidth", 1.5)
set(gca,'FontSize', 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of $u(x)-cos(\frac{{\pi}}{2}x)$ when $k = 256, N = 2^{11}$", "FontSize", 24, "Interpreter", "latex")

% 解函数一阶导数的图象
for k = 1:5
    plot(linspace(a, b, N+1), Ud1_mat_1(:, k), "LineWidth", 1.5)
    hold on
end
set(gca,'FontSize', 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of and $u^\prime(x)$ when $k = 2^4\sim2^8, N = 2^{11}$", "FontSize", 24, "Interpreter", "latex")
legend("$k = 16$", "$k = 32$", "$k = 64$", "$k = 128$", "$k = 256$", "FontSize", 20, "Interpreter", "latex")
hold off

% 解函数二阶导数的图象
for k = 1:5
    plot(linspace(a, b, N+1), Ud2_mat_1(:, k), "LineWidth", 1.5)
    hold on
end
set(gca,'FontSize', 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of and $u^{\prime\prime}(x)$ when $k = 2^4\sim2^8, N = 2^{11}$", "FontSize", 24, "Interpreter", "latex")
legend("$k = 16$", "$k = 32$", "$k = 64$", "$k = 128$", "$k = 256$", "FontSize", 20, "Interpreter", "latex")
hold off

% cos(kx)u 的图象
for k = 1:5
    plot(linspace(a, b, N+1), qU_mat_1(:, k), "LineWidth", 1.5)
    hold on
end
set(gca,'FontSize', 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of and $cos(kx)u(x)$ when $k = 2^4\sim2^8, N = 2^{11}$", "FontSize", 24, "Interpreter", "latex")
legend("$k = 16$", "$k = 32$", "$k = 64$", "$k = 128$", "$k = 256$", "FontSize", 20, "Interpreter", "latex")
hold off

% ------------------探究 N 对 U_err 的影响(N = 2^5 ~ 2^11)------------------

% 参数初始化
syms x
k = 16;
a= -1;
b = 1;
normalize_by = "L2";

% 计算(固定 k 调整 N)
m = 5;
U_mat_2 = zeros(2^11+1, m);
Uerr_mat_2 = zeros(2^11+1, m);
Ud1_mat_2 = zeros(2^11+1, m);
Ud2_mat_2 = zeros(2^11+1, m);
qU_mat_2 = zeros(2^11+1, m);
lambda_vec_2 = 1:m;
for j = 1:m
    n_j = 1:2^(j+6)+1;
    [U_mat_2(n_j, j), lambda_vec_2(j), ~, ~, ~, ~, Ud2_mat_2(n_j, j), ...
        Ud2_mat_2(n_j, j), qU_mat_2(n_j, j), ~, ~, Uerr_mat_2(n_j, j)] = ...
        fem(x, cos(k*x), 2^(j+6), a, b, normalize_by);
end

% 解函数与 cos(pi*x/2) 的图象
for j = 1:m
    n_j = 1:2^(j+6)+1;
    plot(linspace(a, b, 2^(j+6)+1), U_mat_2(n_j, j), "LineWidth", 1.5)
    hold on
end
plot(linspace(a, b, 2^11+1), cos(pi*linspace(a, b, 2^11+1)/2), "LineWidth", 1.5)
set(gca,"FontSize", 18)
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of $cos(\frac{{\pi}}{2}x)$ and $u(x)$ when $N = 2^7\sim2^{11}, k = 16$", "FontSize", 24, "Interpreter", "latex")
legend("$N = 128$", "$N = 256$", "$N = 512$", "$N = 1024$", "$N = 2048$", "$cos(\frac{{\pi}}{2}x)$", "FontSize", 20, "Interpreter", "latex")
hold off

% 比较不同 N 下 U_err 的波动性
for j = 1:m
    n_j = 1:2^(j+6)+1;
    plot(linspace(a, b, 2^(j+6)+1), Uerr_mat_2(n_j, j), "LineWidth", 1.5)
    hold on
end
set(gca,'FontSize', 18);
xlabel("X", "FontSize", 20)
ylabel("Y", "FontSize", 20)
title("Plot of $u(x)-cos(\frac{{\pi}}{2}x)$ when $N = 2^7\sim2^{11}, k = 16$", "FontSize", 24, "Interpreter", "latex")
legend("$N = 128$", "$N = 256$", "$N = 512$", "$N = 1024$", "$N = 2048$", "FontSize", 20, "Interpreter", "latex")
hold off
