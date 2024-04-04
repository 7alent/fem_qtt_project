% 基于有限元方法求解一维二阶椭圆型微分方程特征值问题(狄利克雷边界)
% -u''(x)+q*cos(kx)u(x)=lambda*u(x), a<=x<=b
% u(a)=0; u(b)=0
% q,k > 0 为参数, M 取较大的值以放大解函数振荡程度
% 求(按模)最小特征值对应的解函数


% ------------------探究 u_err 与 k 的关系(k = 2^4 ~ 2^8)-------------------
% 参数初始化
syms x
N = 2^11;
a= -1;
b = 1;
normalize_by = 'L2';
q = 100;

% 计算(固定 N 调整 k)
m = 5;
U_mat_1 = zeros(N+1, m);
Uerr_mat_1 = zeros(N+1, m);
A1U_mat_1 = zeros(N+1, m);
A2U_mat_1 = zeros(N+1, m);
L_mat_1 = zeros(N+1, m);
R_mat_1 = zeros(N+1, m);
lambda_vec_1 = 1:m;
for k = 1:m
    [U_mat_1(:, k), lambda_vec_1(k), L_mat_1(:, k), R_mat_1(:, k), ...
        A1U_mat_1(:, k), A2U_mat_1(:, k), Uerr_mat_1(:, k)] = ...
        fem(x, q*cos(2^(k+3)*x), N, a, b, normalize_by);
end

% 解函数与 cos(pi*x/2) 的图象
for k = 1:5
    plot(linspace(a, b, N+1), U_mat_1(:, k), 'LineWidth', 1.5)
    hold on
end
plot(linspace(a, b, N+1), cos(pi*linspace(a, b, N+1)/2), 'LineWidth', 1.5)
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $cos(\frac{{\pi}}{2}x)$ and $u_h(x)$ when $k = 2^4\sim2^8, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')
legend('$k = 16$', '$k = 32$', '$k = 64$', '$k = 128$', '$k = 256$', '$cos(\frac{{\pi}}{2}x)$', 'FontSize', 20, 'Interpreter', 'latex')
hold off

% 比较不同 k 下 u_err 的波动性
plot(linspace(a, b, N+1), Uerr_mat_1(:, 1), 'LineWidth', 1.5)
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $u_h(x)-cos(\frac{{\pi}}{2}x)$ when $k = 16, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')

plot(linspace(a, b, N+1), Uerr_mat_1(:, 2), 'LineWidth', 1.5)
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $u_h(x)-cos(\frac{{\pi}}{2}x)$ when $k = 32, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')

plot(linspace(a, b, N+1), Uerr_mat_1(:, 3), 'LineWidth', 1.5)
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $u_h(x)-cos(\frac{{\pi}}{2}x)$ when $k = 64, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')

plot(linspace(a, b, N+1), Uerr_mat_1(:, 4), 'LineWidth', 1.5)
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $u_h(x)-cos(\frac{{\pi}}{2}x)$ when $k = 128, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')

plot(linspace(a, b, N+1), Uerr_mat_1(:, 5), 'LineWidth', 1.5)
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $u_h(x)-cos(\frac{{\pi}}{2}x)$ when $k = 256, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')

% A_1*U 的图象
for k = 1:5
    plot(linspace(a, b, N+1), A1U_mat_1(:, k), 'LineWidth', 1.5)
    hold on
end
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $A_1U$ when $k = 2^4\sim2^8, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')
legend('$k = 16$', '$k = 32$', '$k = 64$', '$k = 128$', '$k = 256$', 'FontSize', 20, 'Interpreter', 'latex')
hold off

% A_2*U 的图象
for k = 1:5
    plot(linspace(a, b, N+1), A2U_mat_1(:, k), 'LineWidth', 1.5)
    hold on
end
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $A_2U$ when $k = 2^4\sim2^8, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')
legend('$k = 16$', '$k = 32$', '$k = 64$', '$k = 128$', '$k = 256$', 'FontSize', 20, 'Interpreter', 'latex')
hold off

% A1*U, A2*U, A*U, lambda*B*U 的图象比较(k = 16, N = 2^11)
plot(linspace(a, b, N+1), A1U_mat_1(:, 1), 'LineWidth', 1.5)
hold on
plot(linspace(a, b, N+1), A2U_mat_1(:, 1), 'LineWidth', 1.5)
hold on
plot(linspace(a, b, N+1), L_mat_1(:, 1), 'LineWidth', 1.5)
hold on
plot(linspace(a, b, N+1), R_mat_1(:, 1), 'LineWidth', 1.5)
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $A_1U, A_2U, AU, {\lambda}BU$ when $k = 16, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')
legend('$A_1U$', '$A_2U$', '$AU$', '${\lambda}BU$', 'FontSize', 20, 'Interpreter', 'latex')
hold off

% A1*U, A2*U, A*U, lambda*B*U 的图象比较(k = 32, N = 2^11)
plot(linspace(a, b, N+1), A1U_mat_1(:, 2), 'LineWidth', 1.5)
hold on
plot(linspace(a, b, N+1), A2U_mat_1(:, 2), 'LineWidth', 1.5)
hold on
plot(linspace(a, b, N+1), L_mat_1(:, 2), 'LineWidth', 1.5)
hold on
plot(linspace(a, b, N+1), R_mat_1(:, 2), 'LineWidth', 1.5)
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title('Plot of $A_1U, A_2U, AU, {\lambda}BU$ when $k = 32, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')
legend('$A_1U$', '$A_2U$', '$AU$', '${\lambda}BU$', 'FontSize', 20, 'Interpreter', 'latex')
hold off

%  A*U 与 lambda*B*U 的误差绝对值(k = 16, N = 2^11)
plot(linspace(a, b, 2^11+1), L_mat_1(:, 1)-R_mat_1(:, 1))
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('$AU-{\lambda}BU$', 'FontSize', 20, 'Interpreter', 'latex')
title('Plot of $AU-{\lambda}BU$ when $k = 16, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')

%  A*U 与 lambda*B*U 的误差绝对值(k = 32, N = 2^11)
plot(linspace(a, b, 2^11+1), L_mat_1(:, 2)-R_mat_1(:, 2))
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('$AU-{\lambda}BU$', 'FontSize', 20, 'Interpreter', 'latex')
title('Plot of $AU-{\lambda}BU$ when $k = 32, N = 2^{11}$', 'FontSize', 24, 'Interpreter', 'latex')


% ------------------探究 N 对 u_err 的影响(N = 2^3 ~ 2^11)------------------

% 参数初始化
syms x
k = 16;
a= -1;
b = 1;
normalize_by = 'L2';

% 计算(固定 k 调整 N)
m = 9;
n = 2;
U_mat_2 = zeros(2^(m+n)+1, m);
Uerr_mat_2 = zeros(2^(m+n)+1, m);
A1U_mat_2 = zeros(2^(m+n)+1, m);
A2U_mat_2 = zeros(2^(m+n)+1, m);
lambda_vec_2 = 1:m;
for j = 1:m
    n_j = 1:2^(j+n)+1;
    [U_mat_2(n_j, j), lambda_vec_2(j), ~, ~, A1U_mat_2(n_j, j), ...
        A2U_mat_2(n_j, j), Uerr_mat_2(n_j, j)] = ...
        fem(x, q*cos(k*x), 2^(j+n), a, b, normalize_by);
end

% 解函数与 cos(pi*x/2) 的图象
for j = 1:m
    n_j = 1:2^(j+n)+1;
    plot(linspace(a, b, 2^(j+n)+1), U_mat_2(n_j, j), 'LineWidth', 1.5)
    hold on
end
plot(linspace(a, b, 2^(m+n)+1), cos(pi*linspace(a, b, 2^(m+n)+1)/2), 'LineWidth', 1.5)
set(gca, 'FontSize', 18)
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title(['Plot of $cos(\frac{{\pi}}{2}x)$ and $u_h(x)$ when $N = 2^{', num2str(1+n), '}\sim2^{', num2str(m+n), '}, k = ', num2str(k), '$'], "FontSize", 24, "Interpreter", "latex")
legend_cell = cell(1, m+1);
for j = 1:m
    legend_cell(1, j) = {['$N = 2^{', num2str(j+n), '}$']};
end
legend_cell(1, m+1) = {'$cos(\frac{{\pi}}{2}x)$'};
legend(legend_cell, 'FontSize', 20, 'Interpreter', 'latex')
hold off

% 比较不同 N 下 u_err 的波动性
for j = 1:m
    n_j = 1:2^(j+n)+1;
    plot(linspace(a, b, 2^(j+n)+1), Uerr_mat_2(n_j, j), 'LineWidth', 1.5)
    hold on
end
set(gca, 'FontSize', 18);
xlabel('X', 'FontSize', 20)
ylabel('Y', 'FontSize', 20)
title(['Plot of $u_h(x)-cos(\frac{{\pi}}{2}x)$ when $N = 2^{', num2str(1+n), '}\sim2^{', num2str(m+n), '}, k = ', num2str(k),'$'], "FontSize", 24, "Interpreter", "latex")
legend_cell = cell(1, m);
for j = 1:m
    legend_cell(1, j) = {['$N = 2^{', num2str(j+n), '}$']};
end
legend(legend_cell, 'FontSize', 20, 'Interpreter', 'latex')
hold off

% log(||u||_L2) 与 log(h) 的折线图(以 N = 2^11 的求解结果作为参照)
u_L2_err_2 = 1:m-1;
for j = 1:m-1
    n_j = 1:2^(j+n)+1;
    u_L2_err_2(j) = fem_norm(U_mat_2(n_j, j), U_mat_2(1:2^(m+n)+1, m), 'L2');
end
plot(log((b-a)./2.^((1:m-1)+n)), log(u_L2_err_2), 'LineWidth', 1.5, 'Marker', '*', 'MarkerSize', 16)
hold on
plot(log((b-a)./2.^((1:m-1)+n)), 2*log((b-a)./2.^((1:m-1)+n))-5, 'LineWidth', 1.5)
set(gca, 'FontSize', 18);
xlabel('$log(h)$', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$log({\Vert}u(x)-u(x)\Vert_{L_2})$', 'FontSize', 20, 'Interpreter', 'latex')
title('Plot of $log({\Vert}u_h(x)-\widetilde{u}(x)\Vert_{L_2})$ when $log(h)$ changes ($k = 16$)', 'FontSize', 24, 'Interpreter', 'latex')
text(log((b-a)/2^(m-1+n)), log(u_L2_err_2(1))-.5, ['$\widetilde{u}(x)$ is the $u_h(x)$ when $N = 2^{11}, k = ', num2str(k), '$'], 'FontSize', 20, 'Interpreter', 'latex')
text(log((b-a)/2^(m-1+n)), log(u_L2_err_2(1))-2, 'Slope of the straight line is $2$', 'FontSize', 20, 'Interpreter', 'latex')
hold off

% log(|lambda|) 与 log(h) 的折线图(以 N = 2^11 的求解结果作为参照)
lambda_err_2 = 1:m-1;
for j = 1:m-1
    lambda_err_2(j) = abs(lambda_vec_2(j)-lambda_vec_2(m));
end
plot(log((b-a)./2.^((1:m-1)+n)), log(lambda_err_2), 'LineWidth', 1.5, 'Marker', '*', 'MarkerSize', 16)
hold on
plot(log((b-a)./2.^((1:m-1)+n)), 2*log((b-a)./2.^((1:m-1)+n)), "LineWidth", 1.5)
set(gca,'FontSize', 18);
xlabel('$log(h)$', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$log(\vert\lambda-\widetilde{\lambda}\vert)$', 'FontSize', 20, 'Interpreter', 'latex')
title(['Plot of $log(\vert\lambda-\widetilde{\lambda}\vert)$ when $log(h)$ changes ($k = ', num2str(k), '$)'], 'FontSize', 24, 'Interpreter', 'latex')
text(log((b-a)/2^(m-1+n)), log(u_L2_err_2(1))+1, ['$\widetilde{\lambda}$ is the $\lambda$ when $N = 2^{11}, k = ', num2str(k), '$'], 'FontSize', 20, 'Interpreter', 'latex')
text(log((b-a)/2^(m-1+n)), log(u_L2_err_2(1)), 'Slope of the straight line is $2$', 'FontSize', 20, 'Interpreter', 'latex')
hold off

