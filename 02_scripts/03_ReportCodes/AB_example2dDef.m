%% MESH DEFORMATION AND VOLUME CHANGE VISUALIZATION
clear; clc; close all;

%% 1. FUNCTIONS AND CONFIGURATION
% Support Function Distance 
SF_R = 3; 

% Wendland C2 function (Basis Function)
Phi = @(r) (max(0, 1-abs(r))).^4 .* (4*abs(r) + 1);    

% Norm Function - Euclidean distance
Norm = @(a, b) sqrt(sum((a - b).^2, 3));

% Quadrilateral Area Function (Shoelace Formula)
Quad_Area = @(X, Y) 0.5 * abs(...
    (X(1:end-1,1:end-1).*Y(2:end,1:end-1) - Y(1:end-1,1:end-1).*X(2:end,1:end-1)) + ...
    (X(2:end,1:end-1).*Y(2:end,2:end)     - Y(2:end,1:end-1).*X(2:end,2:end)) + ...
    (X(2:end,2:end).*Y(1:end-1,2:end)     - Y(2:end,2:end).*X(1:end-1,2:end)) + ...
    (X(1:end-1,2:end).*Y(1:end-1,1:end-1) - Y(1:end-1,2:end).*X(1:end-1,1:end-1)));

%% 2. LOADING MESH
% File Handling
MeshFile = fopen("05_meshes\NACA0012257x129.xyz", "r");
Header = fscanf(MeshFile, '%d', 2);
Ni = Header(1);
Nj = Header(2);
MeshData = textscan(MeshFile, '%f %f %f');
fclose(MeshFile);

% Set Aerofoil and Evaluation Points (Rest of Mesh)
Ex = [reshape(MeshData{1}, [Ni*Nj, 1]), reshape(MeshData{3}, [Ni*Nj, 1])];
Ax = [Ex(1:Ni,1), Ex(1:Ni,2)];

% Base Grid Setup and Initial Area Calculation
X_orig_grid = reshape(Ex(:,1), [Ni, Nj]);
Y_orig_grid = reshape(Ex(:,2), [Ni, Nj]);
Area_orig = Quad_Area(X_orig_grid, Y_orig_grid);

%% 3. CONTROL POINTS & INFLUENCE MATRIX
% Sample control points densely along the aerofoil to capture camber accurately
xN_idx = 1:3:Ni; 
xN = [Ex(xN_idx, 1), Ex(xN_idx, 2)];
[xN_N, xN_D] = size(xN);

% Calculate symmetric influence matrix for control points (Cx)
Cx_rows = reshape(xN, [xN_N, 1, xN_D]);
Cx_cols = reshape(xN, [1, xN_N, xN_D]);
Norm_x = Norm(Cx_rows, Cx_cols);
Cx = Phi(Norm_x/SF_R);

% Calculate evaluation matrix for the entire volume mesh (CEx)
[Ex_N, ~] = size(Ex);
CEx_rows = reshape(Ex, [Ex_N, 1, 2]);
CEx_cols = reshape(xN, [1, xN_N, 2]);
Norm_Ex = Norm(CEx_rows, CEx_cols);
CEx = Phi(Norm_Ex/SF_R);

%% 4. DEFORMATION (NACA Camber)
% Deform the control points using your specific function
RxN = NACACamber(xN, 0.09, 0.3); % Adjust camber parameters as needed
DxN = RxN - xN;

% Deform the entire pure aerofoil line for plotting purposes
DAx = NACACamber(Ax, 0.09, 0.3);

% Solve for RBF Weights
Gamma_x = Cx \ DxN(:,1);
Gamma_y = Cx \ DxN(:,2);

% Interpolate deformation onto the entire volume mesh
Sx_x_Mesh = CEx * Gamma_x;
Sx_y_Mesh = CEx * Gamma_y;
Nx_Mesh = [Sx_x_Mesh + Ex(:,1), Sx_y_Mesh + Ex(:,2)];
         
X_def_grid = reshape(Nx_Mesh(:,1), [Ni, Nj]);
Y_def_grid = reshape(Nx_Mesh(:,2), [Ni, Nj]);

% Calculate cell volume change percentage
Area_def = Quad_Area(X_def_grid, Y_def_grid);
C_cell_ratio = ((Area_def - Area_orig) ./ Area_orig) * 100;

%% 5. PLOTTING OUTPUT
fig_main = figure('Name', 'O-Mesh Deformation', 'Units', 'centimeters', 'Position', [5 5, 30, 12]);
tiledlayout(1, 2, "TileSpacing", "compact"); 

% --- TILE 1: Undeformed Mesh ---
ax1 = nexttile; hold on; axis equal;
set(ax1, 'Color', 'k'); % Set black background and white axes

% Plot full mesh with grey lines
mesh(ax1, X_orig_grid, Y_orig_grid, zeros(Ni, Nj), 'EdgeColor', [0.5 0.5 0.5], 'FaceColor', 'none', 'LineWidth', 0.1);

% Plot aerofoil outline in bold red
plot(ax1, Ax(:,1), Ax(:,2), "m-", "LineWidth", 1.5); 

xlim(ax1, [-1, 2]); ylim(ax1, [-1.5, 1.5]); 
xlabel(ax1, "X Coordinate");
ylabel(ax1, "Y Coordinate");
view(ax1, 2); 
set(ax1, "Layer", "Top");

% --- TILE 2: Deformed Mesh (Colored by Volume Ratio) ---
ax2 = nexttile; hold on; axis equal;
set(ax2, 'Color', 'k');
% Surface plot mapping the volume change to color
surf(ax2, X_def_grid, Y_def_grid, zeros(Ni, Nj), C_cell_ratio, ...
    'EdgeColor', 'none', 'FaceColor', 'flat');
% Plot deformed aerofoil outline in bold red
plot(ax2, DAx(:,1), DAx(:,2), "m-", "LineWidth", 1.5); 

% Symmetric color limits centered on 0 for a diverging colormap
max_vol_change = max(abs(C_cell_ratio(:)));
if max_vol_change == 0
    max_vol_change = 1e-6; 
end
clim(ax2, [-max_vol_change, max_vol_change]);

% Add colorbar
c1 = colorbar(ax2);
c1.Label.String = '\Delta Cell Volume (%)';

% Using standard diverging colormap approach (jet/parula/turbo)
colormap(ax2, "jet");  

xlim(ax2, [-1, 2]); ylim(ax2, [-1.5, 1.5]);
xlabel(ax2, "X Coordinate");
ylabel(ax2, "Y Coordinate");
view(ax2, 2);
set(ax2, "Layer", "Top");

% --- APPLY GLOBAL STYLES ---
set(findall(fig_main, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig_main, '-property', 'FontWeight'), 'FontWeight', 'bold');
set(findall(fig_main, '-property', 'FontSize'), 'FontSize', 12);
set(findall(fig_main, 'Type', 'axes'), 'LineWidth', 2);

%% 6. SAVE FIGURES
if ~exist('04_figures', 'dir')
    mkdir('04_figures');
end

date_str = datestr(now, 'yyyymmdd');
save_filename = sprintf('F2_MeshVolumeDeformation_%s', date_str);

% Save as .fig
savefig(fig_main, fullfile('04_figures', [save_filename, '.fig']));
% Save as high-quality PDF
%exportgraphics(fig_main, fullfile('04_figures', [save_filename, '.pdf']), 'ContentType', 'vector');

fprintf('Plot successfully saved to /04_figures/%s\n', save_filename);