function plotPolyL(polyExpr, vars, x_range, y_range, color, level_values)
    % plotSymbolicPolynomialLevelSet - Plot the level set contour of a symbolic polynomial.
    %
    % Parameters:
    %   polyExpr (symbolic) - Symbolic polynomial expression
    %   x_range (vector) - Range of x values [xmin, xmax]
    %   y_range (vector) - Range of y values [ymin, ymax]
    %   level_values (vector) - Contour levels to plot, default is [0]

    % Handle default level_values parameter
    if nargin < 5
        level_values = [0 0]; % Default to plotting the contour f(x, y) = 0
    end

    % Define symbolic variables
    % syms x y

    % Create mesh grid for x and y
    x_vals = linspace(x_range(1), x_range(2), 500);
    y_vals = linspace(y_range(1), y_range(2), 500);
    [xGrid, yGrid] = meshgrid(x_vals, y_vals);
    zGrid = reshape(double(subs(polyExpr, vars, [xGrid(:), yGrid(:)].')), size(xGrid));

    % Plot contour
    % figure; % Create a new figure window
    % hold on;
    contour(xGrid, yGrid, zGrid, level_values, 'LineWidth', 2, 'EdgeColor', color, 'FaceColor', color, 'FaceAlpha', 0.2);
    xlabel('x');
    ylabel('y');
    title('Level Set of the Symbolic Polynomial');
    grid on;
    axis equal; % Ensure equal scaling for x and y axes
end
