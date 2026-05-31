% function to detect all trigonometric terms in the given symbolic expression
function trig_terms = detect_trigonometric_terms(expr)
    % expr: symbolic expression to be analyzed
    % trig_terms: cell array of detected trigonometric terms in the expression

    % get all the terms in the expression using the findSymType function
    % Search for various trigonometric function types
    trig_types = {'sin', 'cos', 'tan', 'cot', 'sec', 'csc', ...
                      'asin', 'acos', 'atan', 'acot', 'asec', 'acsc', ...
                      'sinh', 'cosh', 'tanh', 'coth', 'sech', 'csch', ...
                      'asinh', 'acosh', 'atanh', 'acoth', 'asech', 'acsch'};

    % Pre-allocate a cell array to hold results from each trig type
    n_types = length(trig_types);
    results = cell(n_types, 1);

    % Iterate through each trigonometric function type
    for i = 1:n_types
        terms = findSymType(expr, trig_types{i});

        if ~isempty(terms)
            results{i} = terms;
        end

    end

    % Concatenate all non-empty results and remove duplicates
    trig_terms = vertcat(results{~cellfun(@isempty, results)});

    if ~isempty(trig_terms)
        trig_terms = unique(trig_terms);
    end

end
