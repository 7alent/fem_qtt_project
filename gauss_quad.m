% 基于符号表达式计算高斯积分(5个高斯点)
% x 为积分变量(符号变量)
% integrand 为积分式(关于 x 的符号表达式)
% a, b 为积分下, 上限(标量)

function output = gauss_quad(x, integrand, a, b)
    gauss_point = [-sqrt(5+2*sqrt(10/7))/3; -sqrt(5-2*sqrt(10/7))/3; 0; ...
        sqrt(5-2*sqrt(10/7))/3; sqrt(5+2*sqrt(10/7))/3];
    weight = [(322-13*sqrt(70))/900; (322+13*sqrt(70))/900; 128/225; ...
        (322+13*sqrt(70))/900; (322-13*sqrt(70))/900];
    output = double(dot(weight, ...
        subs(integrand, x, (b-a)/2*gauss_point+(b+a)/2))*(b-a)/2);
end
