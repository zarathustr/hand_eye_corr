function vA = wedge(A, dim)

    len = dim * (dim - 1) / 2;
    if(isnumeric(A(1, 1)))
        S = zeros(dim, dim);
        S_idx = skewdec(dim, 1);
        S_sign = zeros(dim, dim);
        for j = 1 : dim
            for i = j + 1 : dim
                idx = abs(S_idx(i, j)) - 1;
                S(j, i) = len - idx + 1;
                S(i, j) = len - idx + 1;
                S_sign(j, i) = (-1)^(i - j);
                S_sign(i, j) = - (-1)^(i - j);
            end
        end

        vA = zeros(len, 1);
        for j = 1 : dim
            for i = j + 1 : dim
                vA(S(j, i)) = A(j, i) * S_sign(j, i);
            end
        end
    else
        S = sym(zeros(dim, dim));
        S_idx = skewdec(dim, 1);
        S_sign = sym(zeros(dim, dim));
        for j = 1 : dim
            for i = j + 1 : dim
                idx = abs(S_idx(i, j)) - 1;
                S(j, i) = sym(len - idx + 1);
                S(i, j) = sym(len - idx + 1);
                S_sign(j, i) = sym((-1)^(i - j));
                S_sign(i, j) = sym(- (-1)^(i - j));
            end
        end

        vA = sym(zeros(len, 1));
        for j = 1 : dim
            for i = j + 1 : dim
                vA(S(j, i)) = A(j, i) * S_sign(j, i);
            end
        end
    end

end