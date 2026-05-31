% helper function to convert polynomial object to symbolic expression
function sym_out = poly2sym(poly_in, pvar_vars, sym_vars)
    % Convert polynomial object to symbolic expression
    %
    % INPUTS:
    % poly_in: polynomial class object (pvar)
    % pvar_vars: vector of pvar variables corresponding to the symbolic variables
    % sym_vars: vector of symbolic variables corresponding to the pvar variables
    %
    % OUTPUTS:
    % sym_out: symbolic expression equivalent of poly_in

    if ~isa(poly_in, 'polynomial')
        error('Input must be of type polynomial');
    end

    % Convert polynomial to string representation
    poly_str = char(poly_in);

    % Replace pvar variable names with symbolic variable names in the string
    for i = 1:length(pvar_vars)
        poly_str = strrep(poly_str, char(pvar_vars(i)), char(sym_vars(i)));
    end

    % Convert the modified string back to a symbolic expression
    sym_out = str2sym(poly_str);
end
