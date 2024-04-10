% Solve the case of q(x) = 10000*cos(2^8*x)

% Initialization
syms x
a = -1;
b = 1;
q = 10000*cos(2^8*x);
l = 10;
N = 2^l+1;
tol = 1e-6;

% Assemble Stiffness Matrices
[A, B] = fem_mat(x, q, N, a, b);

% Solve Eigenvalue Problem
[U, lambda, t] = qtt_eig(A, B, 738, tol^2, tol);

% Normalization
U = [[0, 0]; U; [0, 0]]; % The first & last entries of the solution is 0 due to the boundary condition
U(:, 1) = U(:, 1)/fem_norm(U(:, 1), 0.*U(:, 1), 'L2');
U(:, 2) = U(:, 2)/fem_norm(U(:, 2), 0.*U(:, 2), 'L2');

% Error Analysis
lambda_err = abs((lambda(1)-lambda(2))/lambda(2));
U_err = fem_norm(U(:, 1), U(:, 2), 'L2');

% Plot
plot(linspace(a, b, N+1), U(:, 2), 'LineWidth', 1.5)
hold on
plot(linspace(a, b, N+1), U(:, 1), 'LineWidth', 1.5)
set(gca, 'FontSize', 18);
xlabel('$X$', 'FontSize', 20, 'Interpreter', 'latex')
ylabel('$Y$', 'FontSize', 20, 'Interpreter', 'latex')
title('Plot of $u_h(x)$ when $q(x) = 10000cos(2^8x)$', 'FontSize', 24, 'Interpreter', 'latex')
legend('Sparse Matrix Eigen Solver', 'DMRG Eigen Solver', 'FontSize', 20, 'Interpreter', 'latex')
hold off







