% helper function to convert symbolic polynomial to pvar polynomial without dpvar
function pvar_poly = sym2pvar(sym_expr, sym_vars, pvar_vars)
    % Convert a symbolic expression to a pvar polynomial
    %
    % INPUTS:
    % sym_expr: symbolic expression (must be a polynomial)
    % sym_vars: symbolic variables used in the expression (column vector)
    % pvar_vars: pvar variables to use in the result (column vector)
    %
    % OUTPUTS:
    % pvar_poly: pvar polynomial object

    % Expand the expression to ensure it's in polynomial form
    sym_expr = expand(sym_expr);

    % Initialize the pvar polynomial
    pvar_poly = polynomial(0);

    % Get all terms by converting to char and parsing, or use children
    % Get the terms of the polynomial (each additive term)
    terms = children(sym_expr);

    % If children returns empty or the expression itself, it's a single term
    if isempty(terms) || (length(terms) == 1 && isequal(terms{1}, sym_expr))
        terms = {sym_expr};
    end

    n_vars = length(sym_vars);

    for t = 1:length(terms)
        term_sym = terms{t};

        % Extract coefficient and powers for this term
        % Get the coefficient (the numeric part)
        [coeff_sym, mono_sym] = coeffs(term_sym, sym_vars);

        if isempty(coeff_sym)
            % It's a constant
            pvar_poly = pvar_poly + double(term_sym);
        else

            for k = 1:length(coeff_sym)
                coeff_val = double(coeff_sym(k));
                pvar_mono = coeff_val;

                % Get degree of each variable in this monomial
                for j = 1:n_vars
                    deg = polynomialDegree(mono_sym(k), sym_vars(j));

                    if deg > 0
                        pvar_mono = pvar_mono * pvar_vars(j) ^ deg;
                    end

                end

                pvar_poly = pvar_poly + pvar_mono;
            end

        end

    end

end
