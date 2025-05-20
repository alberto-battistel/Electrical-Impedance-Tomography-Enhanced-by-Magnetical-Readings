function [x] = TSVD(A, b, trunction_idx)

% Input: A (m x n), b (m x 1)
[U, S, V] = svd(A, 'econ');
s = diag(S);

s = s(1:trunction_idx);

Ut_b = U' * b;

% Tikhonov filter factors
filt = s ./ (s.^2);
x = V(:,1:length(s)) * (filt .* Ut_b(1:length(s)));


end