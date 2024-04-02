% 计算FEM求得的两个解向量对应解函数的误差范数
% V, W 为解向量(二维数组)
% P 为节点坐标信息矩阵, 第n列为第n个节点的坐标
% norm_type 为范数类型: "L2" 或 "H1"
% 若要求V(W)的范数, 设定W(V)为0向量即可
% 暂时只支持[-1, 1]上的函数

function output = fem_norm(W, V, norm_type)
    % 将输入向量转为列向量
    if size(W, 1) < 2
        W = W';
    end
    if size(V, 1) < 2
        V = V';
    end
    
    % 基于线性插值补全较短的向量
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

    % 变量初始化
    N_m = size(U, 1);
    h = 2/(N_m-1);
    
    % 基于高斯积分的方法(H1模代码待修改: 先将基函数化为标准区间[-1,1]再对新自变量求导并乘上jacobi阵)
    gauss_point = [-sqrt(5+2*sqrt(10/7))/3; -sqrt(5-2*sqrt(10/7))/3; 0; ...
        sqrt(5-2*sqrt(10/7))/3; sqrt(5+2*sqrt(10/7))/3];
    weight = [(322-13*sqrt(70))/900; (322+13*sqrt(70))/900; 128/225; ...
        (322+13*sqrt(70))/900; (322-13*sqrt(70))/900];
    index = fix((gauss_point+1)/h)+1; % 距离每个高斯点最近的左侧格点编号
    left_node = -1+(index-1)*h; % 距离每个高斯点最近的左侧格点坐标
    slope = (U(index+1)-U(index))/h; % 高斯点所在单元的解函数斜率
    U_gauss = (gauss_point-left_node).*slope+U(index); % 解函数在高斯点处的值
    output = dot(weight, U_gauss.^2);
    if norm_type == "H1"
        U_d = U; % (近似)一阶导数
        U_d(2:N_m-1) = (U_d(3:N_m)-U_d(1:N_m-2))/(2*h);
        U_d(1) = U_d(2);
        U_d(N_m) = U_d(N_m-1);
        output = sqrt(output+fem_norm(U_d, 0.*U_d, "L2")^2);
    else
        output = sqrt(output);
    end
    
end
