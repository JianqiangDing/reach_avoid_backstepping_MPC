% Iteratively sample the FL box until exactly `n_valid` points lie inside the
% reach-avoid region (safe >= 0, outside target, V_vanilla >= 0). This is the
% set of points that will be added as per-sample SOP constraints, so this is the
% knob a user should set directly (rather than a pre-filter initial-sample count).
%
%   x_samples = sample_n_valid(n_valid, x_vars, y_vars, hx, safe_set, target_set,
%                              certificate, bound_min, bound_max[, opts])
%
% opts (struct, all optional):
%   max_iter            cap on retries (default 30)
%   initial_batch_factor first batch size = n_valid * factor (default 3.0)
%   min_batch           minimum batch size per iteration (default 200)
%   safety_factor       batch oversize factor in retry iterations (default 1.5)
%   verbose             print per-iteration progress (default false)
%
% Returns x_samples as an n_x × n_valid matrix (each column is one sample), to
% match what solve_k1_controller_sop_*.m expects.
function x_samples = sample_n_valid(n_valid, x_vars, y_vars, hx, safe_set, target_set, certificate, bound_min, bound_max, opts)
    if nargin < 10 || isempty(opts); opts = struct(); end
    if ~isfield(opts, 'max_iter');             opts.max_iter = 30;             end
    if ~isfield(opts, 'initial_batch_factor'); opts.initial_batch_factor = 3.0; end
    if ~isfield(opts, 'min_batch');            opts.min_batch = 200;           end
    if ~isfield(opts, 'safety_factor');        opts.safety_factor = 1.5;       end
    if ~isfield(opts, 'verbose');              opts.verbose = false;           end

    % Build the three test functions once (matlabFunction is the slow step).
    safe_set_x    = subs(safe_set,    y_vars, hx);
    target_set_x  = subs(target_set,  y_vars, hx);
    certificate_x = subs(certificate, y_vars, hx);
    safe_func   = matlabFunction(safe_set_x,    'Vars', x_vars);
    target_func = matlabFunction(target_set_x,  'Vars', x_vars);
    cert_func   = matlabFunction(certificate_x, 'Vars', x_vars);

    n = length(x_vars);
    accumulated = zeros(0, n);
    total_tried = 0;
    iter = 0;

    while size(accumulated, 1) < n_valid && iter < opts.max_iter
        iter = iter + 1;
        if iter == 1
            batch = max(ceil(n_valid * opts.initial_batch_factor), opts.min_batch);
        else
            yield = size(accumulated, 1) / max(total_tried, 1);
            remaining = n_valid - size(accumulated, 1);
            est = ceil(remaining / max(yield, 1e-3) * opts.safety_factor);
            batch = max(est, opts.min_batch);
        end

        candidates = bound_min' + (bound_max - bound_min)' .* rand(batch, n);
        args = num2cell(candidates', 2);
        sv = safe_func(args{:});
        tv = target_func(args{:});
        cv = cert_func(args{:});
        if isscalar(sv); sv = sv * ones(batch, 1); end
        if isscalar(tv); tv = tv * ones(batch, 1); end
        if isscalar(cv); cv = cv * ones(batch, 1); end
        valid_mask = (sv(:) >= 0) & (tv(:) >= 0) & (cv(:) >= 0);
        accumulated = [accumulated; candidates(valid_mask, :)]; %#ok<AGROW>
        total_tried = total_tried + batch;

        if opts.verbose
            fprintf('  sample_n_valid iter=%d  batch=%d  valid_this=%d  valid_total=%d/%d  total_tried=%d\n', ...
                iter, batch, sum(valid_mask), size(accumulated, 1), n_valid, total_tried);
        end
    end

    if size(accumulated, 1) < n_valid
        warning('sample_n_valid:notEnough', ...
            'Only %d of %d requested valid samples after %d iterations (tried %d total). Reach-avoid region may be very small at this (mu, xi0).', ...
            size(accumulated, 1), n_valid, opts.max_iter, total_tried);
    else
        accumulated = accumulated(1:n_valid, :);
    end

    overall_yield = size(accumulated, 1) / max(total_tried, 1);
    fprintf('sample_n_valid: collected %d valid (target %d) in %d iter, total tried %d, yield %.1f%%\n', ...
        size(accumulated, 1), n_valid, iter, total_tried, 100*overall_yield);

    x_samples = accumulated';   % n_x × n_valid, matches the SOP solver's expectation
end
