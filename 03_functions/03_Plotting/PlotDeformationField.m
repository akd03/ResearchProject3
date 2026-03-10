function PlotDeformationField(Mx, DMx, skip_step)
%PLOTDEFORMATIONFIELD Creates a quiver plot showing mesh vertex displacements.
%   Mx: Original mesh points (Nm x 2)
%   DMx: Deformed absolute mesh points (Nm x 2)
%   skip_step: Integer to subsample points for readability (default is 50)

arguments
    Mx (:,2) double
    DMx (:,2) double
    skip_step (1,1) double = 50
end

% 1. Calculate the displacement vectors (U, V)
% Subtracting original positions from deformed positions yields the delta
U = DMx(:, 1) - Mx(:, 1);
V = DMx(:, 2) - Mx(:, 2);

% 2. Subsample the data to prevent visual overcrowding
idx = 1:skip_step:size(Mx, 1);

% 3. Generate the plot
figure;
hold on;
axis equal;
title('Mesh Vertex Deformation Vectors');
xlabel('X');
ylabel('Y');

% Plot a subset of the original mesh points lightly for origin context
plot(Mx(idx, 1), Mx(idx, 2), '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 3);

% Plot the actual displacement vectors
% The '0' argument forces the arrow lengths to equal the exact displacement magnitude
quiver(Mx(idx, 1), Mx(idx, 2), U(idx), V(idx), 0, 'r', 'LineWidth', 1);

hold off;
end