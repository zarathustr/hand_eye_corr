function A = times_(r, dim)

len = dim * (dim - 1) / 2;

if(isnumeric(r(1)))
    A = zeros(dim, dim);
    S_idx = skewdec(dim, 1);
    for j = 1 : dim
        for i = j + 1 : dim
            idx = abs(S_idx(i, j)) - 1;
            A(j, i) = (-1)^(i - j) * r(len - idx + 1);
            A(i, j) = - (-1)^(i - j) * r(len - idx + 1);
        end
    end
else
    A = sym(zeros(dim, dim));
    S_idx = skewdec(dim, 1);
    for j = 1 : dim
        for i = j + 1 : dim
            idx = abs(S_idx(i, j)) - 1;
            A(j, i) = sym((-1)^(i - j)) * r(len - idx + 1);
            A(i, j) = sym(- (-1)^(i - j)) * r(len - idx + 1);
        end
    end
end

end