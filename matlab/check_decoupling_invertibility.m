% Check invertibility of the decoupling matrix A(x) over a sampled region, to
% assess whether a control-affine system is (input-output) feedback
% linearizable there. A(x) is the matrix in y^(r) = A(x) u + b(x); the feedback
% law u = A(x)^{-1} b(x) is only well-posed where A(x) is nonsingular, so a
% near-singular A is exactly where the synthesized controller blows up.
%
% Reuses reach_avoid_controller to obtain A_matrix (it only depends on f, g, h).
%
% INPUTS
%   name        : label for the system (string), used in the printout
%   fx,gx,hx    : symbolic f(x) (n x 1), g(x) (n x m), h(x) (p x 1)
%   x_vars      : symbolic state vector (n x 1)
%   y_vars      : symbolic output vector (p x 1)
%   safe_set    : symbolic safe-set polynomial in y (>=0 inside); used only to
%                 flag which samples lie in the safe set (region of interest)
%   bound_min/max: sampling box for the state (n x 1)
%   n_samples   : number of random samples
%   out_csv     : path to write the per-sample CSV (state cols + det,min_sv,cond,in_safe)
%
% OUTPUT
%   stats : struct with summary fields (also printed)
function stats = check_decoupling_invertibility(name, fx, gx, hx, x_vars, y_vars, ...
        safe_set, bound_min, bound_max, n_samples, out_csv)

    % A(x) only depends on (f, g, h); reach_avoid_controller computes it for us.
    [~, ~, ~, ~, ~, ~, ~, A_matrix, ~, ~, p, ~] = ...
        reach_avoid_controller(fx, gx, hx, x_vars, y_vars, safe_set);

    n = numel(x_vars);
    m = size(gx, 2);

    % Vectorized evaluator for each entry of A (handles trig terms via matlabFunction).
    entry_fun = cell(p, m);
    for i = 1:p
        for j = 1:m
            entry_fun{i, j} = matlabFunction(A_matrix(i, j), 'Vars', x_vars);
        end
    end
    safe_fun = matlabFunction(subs(safe_set, y_vars, hx), 'Vars', x_vars);

    rng(0);
    X = bound_min(:)' + (bound_max(:)' - bound_min(:)') .* rand(n_samples, n);
    Xc = num2cell(X, 1);  % one column-vector per state variable (vectorized args)

    % Evaluate every A entry over all samples at once.
    Aent = cell(p, m);
    for i = 1:p
        for j = 1:m
            v = entry_fun{i, j}(Xc{:});
            if isscalar(v); v = v * ones(n_samples, 1); end
            Aent{i, j} = v(:);
        end
    end
    sv_safe = safe_fun(Xc{:});
    if isscalar(sv_safe); sv_safe = sv_safe * ones(n_samples, 1); end
    in_safe = sv_safe(:) >= 0;

    detA = nan(n_samples, 1);
    min_sv = nan(n_samples, 1);
    condA = nan(n_samples, 1);
    for k = 1:n_samples
        A = zeros(p, m);
        for i = 1:p
            for j = 1:m
                A(i, j) = Aent{i, j}(k);
            end
        end
        s = svd(A);
        min_sv(k) = min(s);
        condA(k) = max(s) / max(min(s), eps);
        if p == m
            detA(k) = det(A);
        end
    end

    % Write per-sample CSV.
    T = array2table([X, detA, min_sv, condA, double(in_safe)]);
    state_names = arrayfun(@(k) sprintf('x%d', k), 1:n, 'UniformOutput', false);
    T.Properties.VariableNames = [state_names, {'det', 'min_sv', 'cond', 'in_safe'}];
    writetable(T, out_csv);

    % Summary over the full box and restricted to the safe set.
    tol_sv = 1e-3;          % "near-singular" threshold on the smallest singular value
    summarize = @(mask) struct( ...
        'count',        sum(mask), ...
        'min_min_sv',   min(min_sv(mask)), ...
        'max_cond',     max(condA(mask)), ...
        'frac_singular', mean(min_sv(mask) < tol_sv));

    s_box  = summarize(true(n_samples, 1));
    s_safe = summarize(in_safe);

    fprintf('\n=== Decoupling matrix A(x) invertibility: %s (A is %dx%d) ===\n', name, p, m);
    fprintf('  samples=%d, in-safe-set=%d\n', n_samples, sum(in_safe));
    fprintf('  [full box ]  min(min_sv)=%.3e  max(cond)=%.3e  near-singular(min_sv<%.0e)=%.2f%%\n', ...
        s_box.min_min_sv, s_box.max_cond, tol_sv, 100 * s_box.frac_singular);
    if sum(in_safe) > 0
        fprintf('  [safe set ]  min(min_sv)=%.3e  max(cond)=%.3e  near-singular=%.2f%%\n', ...
            s_safe.min_min_sv, s_safe.max_cond, 100 * s_safe.frac_singular);
        verdict = 'YES (well-conditioned over safe set)';
        if s_safe.min_min_sv < tol_sv
            verdict = sprintf('NO (A is near-singular on %.1f%% of the safe set)', 100 * s_safe.frac_singular);
        end
        fprintf('  feedback linearizable over safe set?  %s\n', verdict);
    end
    fprintf('  wrote %s\n', out_csv);

    stats = struct('name', name, 'box', s_box, 'safe', s_safe, 'csv', out_csv);
end
