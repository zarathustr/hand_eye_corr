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

function R = orthonormalize(A)
[u, ~, v] = svd(A);
s = size(A);
S = eye(s(1));
S(s(1), s(1)) = det(u * v);
R = u * S * v';
end
