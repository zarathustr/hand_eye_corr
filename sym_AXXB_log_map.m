% Correspondence Matching for Hand-eye Calibration
% Author: Jin Wu
% E-mail: jin_wu_uestc@hotmail.com
% Website: www.jinwu.science
%
% Citation: Wu, J. and Liu, M. (2019) 
%                 Correspondence Matching for Hand-eye Calibration. 
%                 IEEE Trans. Instrum. Meas. (Submitted)
%
% Copyright (c) 2019 Jin Wu

clear all
close all
clc

format long g

% global dim len AA A

dim = 4;
len = dim * (dim - 1) / 2;

str1 = 'syms';
str2 = 'syms';
for i = 1 : len
    str1 = strcat(str1, sprintf(' a%d', i));
    str2 = strcat(str2, sprintf(' b%d', i));
end
str1 = strcat(str1, ' real');
str2 = strcat(str2, ' real');
eval(str1);
eval(str2);

for i = 1 : dim
    str1 = 'syms';
    for j = 1 : dim
        str1 = strcat(str1, sprintf(' x%d%d', i, j));
    end
    str1 = strcat(str1, ' real');
    eval(str1);
end

str1 = 'X = [';
for i = 1 : dim
    for j = 1 : dim
        if(j ~= dim)
            str1 = strcat(str1, sprintf('x%d%d,', i, j));
        else
            str1 = strcat(str1, sprintf('x%d%d;', i, j));
        end
    end
end
str1 = strcat(str1, '];');
eval(str1);

str1 = 'a = [';
str2 = 'b = [';
for i = 1 : len
    str1 = strcat(str1, sprintf(' a%d;', i));
    str2 = strcat(str2, sprintf(' b%d;', i));
end
str1 = strcat(str1, '];');
str2 = strcat(str2, '];');
eval(str1);
eval(str2);

alpha = times_(a, dim);
beta = times_(b, dim);

[A, B] = equationsToMatrix(wedge(alpha, dim) == wedge(X * beta * X', dim), b);

res = wedge(alpha, dim) - wedge(X * beta * X', dim);
simplify(A * b - B - res)'

xx = reshape(X, [dim * dim, 1]);
x = xx;

num = 1;

for k = 1 : num
    XX = randn(dim, dim);
    XX = orthonormalize(XX);

    for i = 1 : dim
        for j = 1 : dim
            eval(sprintf('x%d%d = XX(%d, %d);', i, j, i, j));
        end
    end

    AA = eval(A);
end
sgn = sign(det(AA));
noise = 1e-1;
eq = sgn * A - orthonormalize(sgn * AA + randn(len, len) * noise);
eq_vec = reshape(eq, [len * len, 1]);
Jacob = jacobian(eq_vec, x);
residual = eq_vec - Jacob * x;
mu = 0.1;
iter = 1000;

xx = randn(dim * dim, 1);
for k = 1 : iter
    K = reshape(xx, [dim, dim]);
    for i = 1 : dim
        for j = 1 : dim
            eval(sprintf('x%d%d = K(%d, %d);', i, j, i, j));
        end
    end
    
    JJ = eval(Jacob);
    eq_vec_num = eval(eq_vec);
    xx = xx - inv(JJ' * JJ + mu * eye(dim * dim)) * JJ' * eq_vec_num;
    
    if(k > 1)
        if(norm(xx - last_xx) < 1e-16)
            k
            break;
        end
    end
    K = orthonormalize(sign(det(K)) * K);
    last_xx = xx;
end

XXX = orthonormalize(sign(det(K)) * K)
XX


function f = J(x)
global dim len AA A
X = reshape(x, [dim, dim]);

for i = 1 : dim
    for j = 1 : dim
        eval(sprintf('x%d%d = X(%d, %d);', i, j, i, j));
    end
end

AAA = eval(A);

f = 0;
for i = 1 : len
    for j = 1 : len
        f = f + (AAA(i, j) - AA(i, j));
    end
end
end
