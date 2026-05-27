% using sampling method to estimate the lower bound & upper bound of given symbolic expression over the zero superlevel set of another symbolic expression
function [estimated_lower_bound, estimated_upper_bound] = compute_poly_bounds_sampling(vars, poly, superlevel_set_poly, samples_num, bound_min, bound_max)
    % sampling samples from the state space within the specified bounds
    num_vars = length(vars);
    bound_min = bound_min(:)'; % ensure row vector (1 x num_vars)
    bound_max = bound_max(:)'; % ensure row vector (1 x num_vars)
    samples = rand(samples_num, num_vars) .* (bound_max - bound_min) + bound_min; % uniformly sample within the bounds
    % split samples into one column-vector cell per variable (alphabetical order, matching matlabFunction)
    sample_cols = num2cell(samples, 1);
    % construct the matlab function for the polynomial and the superlevel set polynomial for evaluation
    poly_func = matlabFunction(poly, 'vars', vars);
    superlevel_set_func = matlabFunction(superlevel_set_poly, 'vars', vars);
    % evaluate the polynomial and the superlevel set polynomial at the sampled points
    poly_values = poly_func(sample_cols{:});
    superlevel_set_values = superlevel_set_func(sample_cols{:});

    % filter the samples that satisfy the superlevel set constraint (superlevel_set_poly >= 0)
    valid_indices = superlevel_set_values >= 0;
    valid_poly_values = poly_values(valid_indices);
    % compute the estimated lower and upper bounds from the valid samples
    estimated_lower_bound = min(valid_poly_values);
    estimated_upper_bound = max(valid_poly_values);
end
