function PlotOMesh(Mx, DMx, Ni, Nj)
%PLOTOMESH Plots the structured O-mesh connectivity before and after deformation.

arguments
    Mx  (:,2) double   % Original volume mesh points
    DMx (:,2) double   % Deformed volume mesh points
    Ni  (1,1) double   % Number of nodes around the aerofoil
    Nj  (1,1) double   % Number of nodes in the radial direction
end

%% FUNCTION BODY
% 1. Reshape the 1D coordinate arrays into Ni x Nj grids
X_orig = reshape(Mx(:, 1), [Ni, Nj]);
Y_orig = reshape(Mx(:, 2), [Ni, Nj]);

X_def = reshape(DMx(:, 1), [Ni, Nj]);
Y_def = reshape(DMx(:, 2), [Ni, Nj]);

% 2. Create a dummy Z-axis matrix for the mesh function
Z_grid = zeros(Ni, Nj);

% 3. Generate the plot
figure;
tiledlayout(1, 2, "TileSpacing", "tight");

% Left Tile: Initial Mesh
nexttile;
hold on;
title('Initial O-Mesh');
xlabel('X');
ylabel('Y');
% Plot the grid lines in gray
mesh(X_orig, Y_orig, Z_grid, 'FaceColor', 'none', 'EdgeColor', [0.6 0.6 0.6]);
% Highlight the aerofoil surface (first column of the reshaped grid)
plot(X_orig(:, 1), Y_orig(:, 1), 'r-', 'LineWidth', 1.5); 
view(2); % Force 2D view directly from above
hold off;

% Right Tile: Deformed Mesh
nexttile;
hold on;
title('Deformed O-Mesh');
xlabel('X');
ylabel('Y');
% Plot the deformed grid lines in gray
mesh(X_def, Y_def, Z_grid, 'FaceColor', 'none', 'EdgeColor', [0.6 0.6 0.6]);
% Highlight the deformed aerofoil surface
plot(X_def(:, 1), Y_def(:, 1), 'r-', 'LineWidth', 1.5); 
view(2);
hold off;

end