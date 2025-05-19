function [res_norms, x_norms, x_lambdas, solutions] = calc_L_curve(A, b, lambdas)

[U, S, V] = svd(A, 'econ');
s = diag(S);
Ut_b = U' * b;

% Define range of lambda values
% lambdas = logspace(-5, 10, 1000);  % adjust as needed
gcv_vals = zeros(size(lambdas));
x_norms = zeros(size(lambdas));
res_norms = zeros(size(lambdas));
x_lambdas = zeros(size(A,2), length(lambdas));
solutions = zeros(length(b), length(lambdas));

for i = 1:length(lambdas)
    lambda = lambdas(i);
    
    % Tikhonov filter factors
    filt = s ./ (s.^2 + lambda^2);
    x_lambda = V * (filt .* Ut_b(1:length(s)));
    x_lambdas(:,i) = x_lambda;
    
    % Compute residual and solution norms
    solution = A * x_lambda;
    solutions(:,i) = solution;
    res = solution - b;
    res_norms(i) = norm(res);
    x_norms(i) = norm(x_lambda);
    
    % GCV computation
    trace_term = sum(s.^2 ./ (s.^2 + lambda^2));
    gcv_vals(i) = norm(res)^2 / (length(b) - trace_term)^2;
end

end