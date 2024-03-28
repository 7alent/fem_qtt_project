% 基于有限元方法求解一维二阶椭圆型微分方程特征值问题(狄利克雷边界)
% -u''(x)+p(x)u'(x)+q(x)u(x)+r(x)=lambda*u(x), a<=x<=b
% u(a)=g_a; u(b)=g_b
% 求(按模)最小特征值对应的解函数

% 求 q=cos(10x) 时 N=2^5~2^9 的解
syms x
k = 10;
q = cos(k*x);
a = -1;
b = 1;
normalize_by = "H1";
lambda_vec = 1:5;
U_mat = zeros(2^9+1, 5);
time_vec = zeros(5, 1);

% 求解并绘制解函数图象
for n = 1:5
    tic;
    [U_mat(1:2^(n+4)+1, n), lambda_vec(n)] = fem(x, q, 2^(n+4), a, b, normalize_by);
    time_vec(n) = time_vec(n)+toc;
    hold on
end
hold off
legend({'2^5', '2^6', '2^7', '2^8', '2^9'})

% 时间-N折线图
plot(2.^(5:9), time_vec)

% 特征值误差-h折线图
plot(2./2.^(5:9), abs(lambda_vec-lambda_vec(5)))
xlabel("h")
ylabel("Error of Lambda")

% 解函数H1模误差-h折线图
solution_err_vec = zeros(5, 1);
for n = 1:4
    solution_err_vec(n) = fem_norm(U_mat(1:2^(n+4)+1, n), U_mat(1:2^9+1, 5), "H1");
end
plot((2./2.^(5:9)).^2, solution_err_vec)
xlabel("h^2")
ylabel("Error of Solution")
