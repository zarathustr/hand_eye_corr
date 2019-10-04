function R = orthonormalize(A)
[u, ~, v] = svd(A);
s = size(A);
S = eye(s(1));
S(s(1), s(1)) = det(u * v);
R = u * S * v';
end