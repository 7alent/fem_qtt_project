% 计算FEM求得的两个解向量对应解函数的误差范数
% V, W 为解向量(必须为列向量)
% P 为节点坐标信息矩阵, 第n列为第n个节点的坐标
% norm_type 为范数类型: "L2" 或 "H1"
% 若要求V(W)的范数, 设定W(V)为0向量即可

function output = fem_norm(W, V, norm_type)
    W_sz = size(W, 1);
    V_sz = size(V, 1);
    if W_sz < V_sz
        U = interp1(linspace(-1, 1, W_sz), W, linspace(-1, 1, V_sz), "linear");
        U = U'-V;
    elseif W_sz > V_sz
        U = interp1(linspace(-1, 1, V_sz), V, linspace(-1, 1, W_sz), "linear");
        U = U'-W;
    else
        U = W-V;
    end
    N_m = size(U, 1);
    P = linspace(-1, 1, N_m);
    h = 2/(N_m-1);
    L2_norm = 0;
    H1_norm = 0;
    syms x
    for n = 1:N_m-1
        k_n = (U(n+1)-U(n))/(P(n+1)-P(n)); % 解函数在某单元上的斜率
        L2_norm = L2_norm+quadgk(matlabFunction( ...
            k_n^2*(x-P(n))^2+2*k_n*U(n)*(x-P(n))+U(n)^2), P(n), P(n+1));
        H1_norm = L2_norm+k_n^2*h;
    end
    if norm_type == "L2"
        output = sqrt(L2_norm);
    else
        output = sqrt(H1_norm);
    end
end