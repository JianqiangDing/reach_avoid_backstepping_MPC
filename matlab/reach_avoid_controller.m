% function for reach-avoid controller synthesis using backstepping

function [u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = reach_avoid_controller(fx, gx, hx, x_vars, y_vars, safe_set)
    % documentation for the function
    % This function synthesizes a reach-avoid controller using backstepping technique.
    % It takes system dynamics, safe set definitions, and target sets as inputs,
    % and outputs a control policy that ensures the system reaches the target while avoiding unsafe states.

    % INPUTS:
    % fx: system dynamics function for state evolution, symbolic matrix expression for f(x) of size (n x 1), n is the state dimension
    % gx: control input function, symbolic matrix expression for g(x) of size (n x m), m is the control input dimension
    % hx: output function, symbolic matrix expression for h(x) of size (p x 1), p is the output dimension
    % x_vars: state variables, symbolic vector of size (n x 1)
    % safe_set: polynomial defining the safe set, such that safe_set >= 0 defines the safe region

    % OUTPUTS:
    % u: synthesized control policy as a symbolic function of the state variables and undefined parameters (mu, lambda, k1 controler polynomial)
    % k1: first layer control vector (symbolic, p x 1)
    % J_k1: Jacobian of k1 w.r.t. y (symbolic, p x p)
    % mu: cell array of mu parameter vectors, mu{i} contains [mu_i_1; mu_i_2; ...] for output i
    % lambda: lambda parameter for backstepping
    % certificate: reach-avoid certificate function
    % A_matrix: A(x) matrix for backstepping
    % b_vector: b(x) vector for backstepping
    % ks: cell array of auxiliary controllers for each output
    % cert_term_dict: dictionary of terms in the reach-avoid certificate
    % p: number of outputs
    % r_deg: vector of relative degrees for each output

    % compute the vector relative degree of the system and cache Lie derivatives
    % Lfh{i} = {h_i, Lf h_i, ..., Lf^{r_i} h_i}
    % LgLfh{i} = {Lg h_i, Lg Lf h_i, ..., Lg Lf^{r_i-1} h_i}
    [r_deg, Lfh, LgLfh] = vector_relative_degree(fx, gx, hx, x_vars);

    % declare symbolic variables for k1 controller as [k1_1, k1_2, ..., k1_p]
    p = length(hx); % number of outputs
    k1 = sym('k1', [p, 1], 'real');
    J_k1 = sym('J_k1', [p, p], 'real'); % Symbolic Jacobian of k1 w.r.t y

    % declare symbolic variables for undefined parameters mu values as a cell array
    % mu{i} contains mu parameters for output i, with r_deg(i)-1 elements
    mu = cell(p, 1);

    for i = 1:p

        if r_deg(i) > 1
            mu{i} = sym(['mu' num2str(i) '_'], [r_deg(i) - 1, 1], 'real');
        else
            mu{i} = [];
        end

    end

    % declare the scale variable lambda for backstepping
    lambda = sym('lambda', 'real');

    % initialize the Dy_psi matrix, which is the partial derivative of the safety contraint function psi(y) w.r.t. y
    Dy_psi = jacobian(safe_set, y_vars);
    Dy_psi = subs(Dy_psi, y_vars, hx); % substitute y with h(x) to get Dy_psi as a function of x

    % compute the A(x) matrix for backstepping controller design
    % A(x) = [L_g L_f^{r1-1} h1; L_g L_f^{r2-1} h2; ... ]
    A_matrix = A_matrix_for_backstepping(r_deg, LgLfh, p);

    % DEBUG
    % print the A_matrix for debugging
    % disp(A_matrix);
    % DEBUG

    % compute the b(x) vector for backstepping controller design
    % Also returns ks: cell array of auxiliary controllers for each output
    [b_vector, ks] = b_vector_for_backstepping(r_deg, Lfh, k1, J_k1, mu, lambda, Dy_psi, p);
    % DEBUG
    % print the b_vector for debugging
    % disp(b_vector);
    % DEBUG

    % compute the control input u = A^{-1} * b
    u = A_matrix \ b_vector;

    % Compute reach-avoid certificate
    [certificate, cert_term_dict] = compute_reach_avoid_certificate(safe_set, ks, mu, r_deg, Lfh, p);

    % Debug output
    % disp('=== A(x) matrix ===');
    % disp(A_matrix);
    % disp('=== b(x) vector ===');
    % disp(b_vector);
    % disp('=== Control u(x) ===');
    % disp(u);
end

% =====================================================
% Compute reach-avoid certificate
% =====================================================
function [certificate, cert_term_dict] = compute_reach_avoid_certificate(safe_set, ks, mu, r_deg, Lfh, p)
    % Synthesize the reach-avoid certificate function for backstepping controller design.
    %
    % INPUTS:
    % safe_set: safe set polynomial S(y) >= 0
    % ks: cell array where ks{i} is a containers.Map of k_l values for output i
    % mu: cell array of mu parameters for each output
    % r_deg: vector relative degree
    % Lfh: cell array of Lie derivatives for all outputs
    % p: number of outputs
    %
    % OUTPUTS:
    % certificate: reach-avoid certificate function
    %              For any state x, if certificate(x) >= 0, then starting from x,
    %              there exists a control to reach the target while staying safe
    % cert_term_dict: dictionary of the terms corresponding to each output in the certificate
    % (for constructing Schur complement conditions in the SOS program later)

    % Initialize the certificate function with the safe set condition S(y) >= 0
    certificate = safe_set;
    cert_term_dict = dictionary("Psi(y(x))", {safe_set}); % left top element of the Schur complement matrix

    % Add the terms for each output
    for i = 1:p
        rd = r_deg(i);

        % Compute the term for this output
        % term = sum from j=1 to rd-1 of: (1 / (2 * mu_{i,j})) * (Lf^j h_i - k_j)^2
        term = sym(0);

        for j = 1:(rd - 1)
            % Lf^j h_i = Lfh{i}{j+1} (1-indexed: Lf^j h_i is at index j+1)
            Lfj_hi = Lfh{i}{j + 1};

            % k_j for output i
            kj = ks{i}(j);

            % mu_{i,j} is the j-th element of mu{i}
            mu_ij = mu{i}(j);

            % Add the term: (1 / (2 * mu_ij)) * (Lfj_hi - kj)^2
            term = term + (1 / (2 * mu_ij)) * (Lfj_hi - kj) ^ 2;

            % save the term information in the term_dict for later reference
            cert_term_dict(sprintf("output_%d_k%d", i, j)) = {[2 * mu_ij, Lfj_hi - kj]};
        end

        % Subtract the term from the certificate
        certificate = certificate - term;
    end

end

% =====================================================
% A(x) matrix for backstepping controller design
% =====================================================
function A = A_matrix_for_backstepping(r_deg, LgLfh, p)
    % Synthesize the A(x) matrix for backstepping controller design.
    %
    % A(x) = [L_g L_f^{r1-1} h1; L_g L_f^{r2-1} h2; ... ]
    %
    % INPUTS:
    % r_deg: vector relative degree, vector of size (p x 1)
    % LgLfh: cell array of Lie derivatives Lg Lf^k h for each output
    % p: number of outputs
    %
    % OUTPUTS:
    % A: the A(x) matrix of size (p x m) where m is control dimension

    A = [];

    for i = 1:p
        rd = r_deg(i);
        % LgLfh{i}{rd} = Lg Lf^{rd-1} h_i (1-indexed, so index rd)
        LgLfh_rd = LgLfh{i}{rd};
        % Append as a new row
        A = [A; LgLfh_rd];
    end

end

% =====================================================
% b(x) vector for backstepping controller design
% =====================================================
function [b, ks] = b_vector_for_backstepping(r_deg, Lfh, k1, J_k1, mu, lambda, Dy_psi, p)
    % Synthesize the b(x) vector for backstepping controller design.
    %
    % INPUTS:
    % r_deg: vector relative degree, vector of size (p x 1)
    % Lfh: cell array of Lie derivatives Lf^k h for each output
    % k1: first layer control vector (p x 1)
    % mu: cell array of mu parameters for each output
    % lambda: the lambda parameter for backstepping
    % Dy_psi: partial derivative of safety constraint w.r.t. output y (1 x p)
    % p: number of outputs
    %
    % OUTPUTS:
    % b: the b(x) vector of size (p x 1)
    % ks: cell array where ks{i} is a containers.Map of k_l values for output i

    b = sym(zeros(p, 1));
    ks = cell(p, 1);

    for i = 1:p
        [bi, ks{i}] = bi_backstepping(r_deg, Lfh, k1, J_k1, mu, lambda, Dy_psi, p, i);
        b(i) = bi;
    end

end

% =====================================================
% Compute b_i for a single output using backstepping
% =====================================================
function [bi, ks] = bi_backstepping(r_deg, Lfh, k1, J_k1, mu, lambda, Dy_psi, p, output_idx)
    % Compute the b_i term for the output_idx-th output.
    %
    % INPUTS:
    % r_deg: vector relative degree
    % Lfh: cell array of Lie derivatives for all outputs
    % k1: first layer control vector (p x 1)
    % mu: cell array of mu parameters for each output
    % lambda: the lambda parameter
    % Dy_psi: jacobian of safety constraint w.r.t output (1 x p)
    % p: number of outputs
    % output_idx: index of the current output (1-indexed)
    %
    % OUTPUTS:
    % bi: the b_i scalar/expression
    % ks: containers.Map of k_l values for this output

    rd = r_deg(output_idx);
    ks = containers.Map('KeyType', 'int32', 'ValueType', 'any');

    % k_1 = k1(output_idx)
    ks(1) = k1(output_idx);

    % If relative degree is 1, then bi = k1_i - Lf^1 h_i
    if rd == 1
        % Lfh{output_idx}{2} = Lf^1 h_i
        bi = k1(output_idx) - Lfh{output_idx}{2};
        return;
    end

    % For rd >= 2, compute auxiliary controllers k_2, k_3, ..., k_rd
    mu_i = mu{output_idx};

    % Compute k_2:
    % k_2 = mu_{i,1} * Dy_psi_i + sum_j (dk1_i/dh_j * Lf h_j) + 0.5*lambda*(Lf h_i - k1_i)
    jacobian_term = sym(0);

    for j = 1:p
        % Lfh{j}{1} = h_j, Lfh{j}{2} = Lf h_j
        h_j = Lfh{j}{1};
        Lf_h_j = Lfh{j}{2};
        jacobian_term = jacobian_term + J_k1(output_idx, j) * Lf_h_j;
    end

    % Lf h_i = Lfh{output_idx}{2}
    Lf_hi = Lfh{output_idx}{2};

    ks(2) = mu_i(1) * Dy_psi(output_idx) + jacobian_term + 0.5 * lambda * (Lf_hi - k1(output_idx));

    % For l = 3 to rd, compute k_l recursively
    for l = 3:rd
        % term_0 = -mu_{i,l-1} / mu_{i,l-2} * (Lf^{l-2} h_i - k_{l-2})
        Lf_l2_hi = Lfh{output_idx}{l - 1}; % Lf^{l-2} h_i (1-indexed: l-2+1 = l-1)
        term_0 = -mu_i(l - 1) / mu_i(l - 2) * (Lf_l2_hi - ks(l - 2));

        % term_1 = sum over s=1 to l-2, sum over j=1 to p of dk_{l-1}/d(Lf^{s-1} h_j) * Lf^s h_j
        term_1 = sym(0);

        for s = 1:(l - 2)

            for j = 1:p
                % Lf^{s-1} h_j = Lfh{j}{s}
                Lf_s1_hj = Lfh{j}{s};
                % Lf^s h_j = Lfh{j}{s+1}
                Lf_s_hj = Lfh{j}{s + 1};
                this_diff = diff(ks(l - 1), Lf_s1_hj);
                term_1 = term_1 + this_diff * Lf_s_hj;
            end

        end

        % term_2 = 0.5 * lambda * (Lf^{l-1} h_i - k_{l-1})
        Lf_l1_hi = Lfh{output_idx}{l}; % Lf^{l-1} h_i (1-indexed: l-1+1 = l)
        term_2 = 0.5 * lambda * (Lf_l1_hi - ks(l - 1));

        ks(l) = term_0 + term_1 + term_2;
    end

    % Finally, bi = k_rd - Lf^rd h_i
    % Lf^rd h_i = Lfh{output_idx}{rd+1}
    Lf_rd_hi = Lfh{output_idx}{rd + 1};
    bi = ks(rd) - Lf_rd_hi;
end

% =====================================================

% function to find the vector relative degree of a control affine system with outputs
function [r_deg, Lfh, LgLfh] = vector_relative_degree(fx, gx, hx, vars)

    % INPUTS:
    % fx: system dynamics function for state evolution, symbolic matrix expression for f(x) of size (n x 1), n is the state dimension
    % gx: control input function, symbolic matrix expression for g(x) of size (n x m), m is the control input dimension
    % hx: output function, symbolic matrix expression for h(x) of size (p x 1), p is the output dimension
    % vars: state variables, symbolic vector of size (n x 1)

    % OUTPUTS:
    % r_deg: vector relative degree, vector of size (p x 1)
    % Lfh: cell array of size (p x 1), where Lfh{i} contains {h_i, Lf h_i, Lf^2 h_i, ..., Lf^{r_i} h_i}
    % LgLfh: cell array of size (p x 1), where LgLfh{i} contains {Lg h_i, Lg Lf h_i, ..., Lg Lf^{r_i-1} h_i}

    p = length(hx); % number of outputs
    r_deg = zeros(p, 1);
    Lfh = cell(p, 1);
    LgLfh = cell(p, 1);

    for i = 1:p
        [r_deg(i), Lfh{i}, LgLfh{i}] = relative_degree_with_cache(fx, gx, hx(i), vars);
    end

end

% function to find the relative degree of a single output with cached Lie derivatives
function [r_deg, Lfh_cache, LgLfh_cache] = relative_degree_with_cache(fx, gx, hi, vars)
    % INPUTS:
    % fx: system dynamics function for state evolution, symbolic matrix expression for f(x) of size (n x 1), n is the state dimension
    % gx: control input function, symbolic matrix expression for g(x) of size (n x m), m is the control input dimension
    % hi: single output function, symbolic expression for h_i(x)
    % vars: state variables, symbolic vector of size (n x 1)

    % OUTPUTS:
    % r_deg: relative degree of the output hi
    % Lfh_cache: cell array containing [h_i, Lf h_i, Lf^2 h_i, ..., Lf^{r_deg} h_i]
    %            Lfh_cache{k} = Lf^{k-1} h_i (1-indexed, so Lfh_cache{1} = h_i)
    % LgLfh_cache: cell array containing [Lg h_i, Lg Lf h_i, ..., Lg Lf^{r_deg-1} h_i]
    %              LgLfh_cache{k} = Lg Lf^{k-1} h_i (1-indexed)

    max_iteration = length(vars) + 1;

    % Initialize caches
    Lfh_cache = {};
    LgLfh_cache = {};

    r = 0;
    Lfh = hi; % start with h_i (corresponds to Lf^0 h_i = h_i)
    Lfh_cache{1} = Lfh; % cache Lf^0 h_i = h_i

    while r < max_iteration
        % Compute Lg Lf^r h_i = jacobian(Lf^r h_i) * g
        LgLfh = jacobian(Lfh, vars) * gx;
        LgLfh_cache{r + 1} = LgLfh; % cache Lg Lf^r h_i (1-indexed)

        % Check if any component of LgLfh is non-zero
        if ~all(LgLfh == 0)
            % Relative degree found: r + 1
            r_deg = r + 1;
            % Compute and cache the final Lf^{r+1} h_i for completeness
            Lfh_next = jacobian(Lfh, vars) * fx;
            Lfh_cache{r + 2} = Lfh_next;
            return;
        end

        % Continue: compute Lf^{r+1} h_i = jacobian(Lf^r h_i) * f
        Lfh = jacobian(Lfh, vars) * fx;
        r = r + 1;
        Lfh_cache{r + 1} = Lfh; % cache the new Lf^r h_i

        if r >= max_iteration
            error('The system does not have a well-defined relative degree for the given output.');
        end

    end

    error('Relative degree not found within maximum iterations.');
end
