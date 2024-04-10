% Solve Generalized Eigenvalues Problem in QTT format
%     [x, lambda] = qtt_eig(A, B, eps, tol)
%
%     Input Arguments:
%         A - Symmetric matrix
%         B - Symmetric positive definite matrix
%         eps - Accuracy of transform a matrix into QTT format
%         tol - Stopping tolerance of eigenvalue solver.
%
%     Output Arguments:
%         x - Eigenvector. 2^d * 2 matrix, first column is the DMRG
%             solution, second column is the eigs() solution.
%         lambda - Eigenvalue. 2 * 1 vector, first entry is the DMRG
%             solution, second entry is the eigs() solution.
%         t - Running time. 2 * 1 vector, first entry is the running
%             time of DMRG method, second entry is the running time of
%             eigs().
%
%     Consider the generalized eigenvalues problem as
%                         A * x = lambda * B * x
%     where A, B are matrices of size 2^d * 2^d , d is a positive integer
%     and A, B are symmetric, B is postive definite;
%     lambda is the smallest eigenvalue, x is the eigenvector respectly.
%     Throgh the Cholesky decomposition of B = LL', the problem can be 
%     rewrote as
%                           C * y = lambda * y
%     with C = inv(L) * A * inv(L'), y = L' * x .
%     If C is not positive definite, we modify C by subtracting m * I 
%     to C and solve the eigenvalue problem
%                  ( C + m * I ) * y = ( lambda + m ) * y
%     here m can be seen as an estimate of -lambda.
%     Noticed that y is the eigenvector w.r.t. eigenvalue ( lambda + m )
%     of matrix ( C + m * I ) .
%     This problem is solved by transforming the matrix into Quantitized
%     Tensor Train (QTT) format and using Density Matrix Renormalization 
%     Group (DMRG) method, the solution of sparse matrix solver (eigs()) 
%     will also be computed.
%     Tips: lambda obtained by DMRG can be use to tune the parameter m,
%           a simple way is to set m = ceil(-lambda).


function [x, lambda, t] = qtt_eig(A, B, m, eps, tol)
    % Initialization
    d = log(size(A, 1))/log(2);
    x = zeros(2^d, 2);
    lambda = zeros(2, 1);
    t = zeros(2, 1);

    % Problem Transformation
    L = chol(B);
    if istriu(L)
        L = L';
    end
    C = (L\A)/(L');
    C = triu(C)+triu(C, 1)';
    
    % Precondition
    C = C+m*eye(2^d);
    try
        chol(C);
    catch
        warning('C is not positive definite, please modify parameter m.');
    end
    
    % DMRG
    C_ttm = tt_matrix(tt_qfromfull(full(C), 2, d, eps, 2));
    tic;
    [x_d, lambda(1), ~] = dmrg_eig(C_ttm, tol);
    t(1) = t(1)+toc;
    x_d = full(x_d);
    x_d = (L')\x_d;
    if sum(x_d >= 0) < size(x_d, 1)/2 % We need a solution with more positive entries(if u is an eigenvector, so is -u)
        x_d = -x_d;
    end
    x(:, 1) = x_d;
    lambda(1) = lambda(1)-m;
    

    % eigs
    tic;
    [x(:, 2), lambda(2)] = eigs(C, 1, 'smallestabs');
    t(2) = t(2)+toc;
    x(:, 2) = (L')\x(:, 2);
    if sum(x(:, 2) >= 0) < size(x(:, 2), 1)/2 % We need a solution with more positive entries(if u is an eigenvector, so is -u)
        x(:, 2) = -x(:, 2);
    end
    lambda(2) = lambda(2)-m;


end