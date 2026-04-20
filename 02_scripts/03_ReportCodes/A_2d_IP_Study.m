%% VARIABLE NAMING CONVENTION                                              
%  Nx       Known (Control) Points                                         
%  Nx_idx   Nx indices in Ax                                               
%  Ax       Aerofoil Mesh Points                                           
%  RAx      Applied Deformation                                            
%  DAx      Actual Aerofoil Position After Rx                              
%  Mx       Volume Mesh Points                                             

clear; clc; close all;

%% LOAD MESH                                                               
% Load NACA Mesh                                                           
MeshFile = fopen("05_meshes\NACA0012257x129.xyz", "r");                    
Header = fscanf(MeshFile, '%d', 2);                                        
Ni = Header(1);
Nj = Header(2);
MeshData = textscan(MeshFile, '%f %f %f');
fclose(MeshFile);

% Set Aerofoil and Evaluation Points (Rest of Mesh)
Mx = [reshape(MeshData{1}, [Ni*Nj, 1]), reshape(MeshData{3}, [Ni*Nj, 1])];
Ax = [Mx(1:Ni,1), Mx(1:Ni,2)];

% Initial Control Points - LE and TE                                       
Nx_idx_base = [1; 129];
 
%% DEFORMATION
% True Deformation
% DAx = Rotate(Ax, 40, 0.5);
% DAx = Translate(DAx, [0, 0]);
DAx = NACACamber(Ax, 0.15, 0.3);
RAx = DAx - Ax;

%% PARAMETER SWEEP SETUP
SF_R = 3;
N = 254; % Fixed total target points
pct_vals = [0, 0.2, 0.4, 0.6, 0.8];
M = size(Ax, 1); 

% --- PREALLOCATE STORAGE FOR PLOTTING OUTSIDE THE LOOP ---
num_tests = length(pct_vals);
all_iterations = cell(num_tests, 1);
all_max_err = cell(num_tests, 1);
all_rmse = cell(num_tests, 1);
legend_labels = strings(num_tests, 1);

%% DATA GENERATION LOOP
fprintf('Running parameter sweep...\n');
for i = 1:num_tests
    pct = pct_vals(i);
    N_IP = round(N * pct);
    N_G = N - N_IP;
    
    % Reset the starting points
    Nx_idx = Nx_idx_base;
    
    % Initial Points (FPS)
    if N_IP > 0
        Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);
    end
    
    % Greedy Algorithm
    if N_G > 0
        [Nx_idx_final, max_err_hist, rmse_hist] = GreedyCholesky(Nx_idx, Ax, RAx, SF_R, 'K', N_G);
        
        % Shift the x-axis for error plots to start exactly at N_IP
        all_iterations{i} = N_IP + (1:length(max_err_hist)) - 1; 
        all_max_err{i} = max_err_hist;
        all_rmse{i} = rmse_hist;
    end
    
    % Format legend text for later
    legend_labels(i) = sprintf('N_{IP} = %d (%d%%)', N_IP, round(pct*100));
end
fprintf('Data generation complete. Plotting...\n');

%% PLOT 1: TRIPLE LAYOUT (SHAPE, MAX ERROR, RMSE)
fig_main = figure('Name', '2D Error Reduction', 'Position', [100, 100, 1400, 450]);
tiledlayout(1, 3, "TileSpacing", "compact");

% --- TILE 1: Aerofoil Shape ---
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');
% Plotting the outline (no internal mesh)
plot(ax1, Ax(:,1), Ax(:,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Original (NACA0012)');
plot(ax1, DAx(:,1), DAx(:,2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Target (Deformed)');
xlabel(ax1, 'X');
ylabel(ax1, 'Y');
title(ax1, 'Aerofoil Deformation Boundary');
legend(ax1, 'Location', 'best');

% --- TILE 2: Max Error ---
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');
set(ax2, 'YScale', 'log'); 
xlabel(ax2, 'Total Control Points Added (k)');
ylabel(ax2, 'Maximum Deformation Error');
title(ax2, sprintf('Max Error Convergence (N = %d)', N));

% --- TILE 3: RMSE ---
ax3 = nexttile;
hold(ax3, 'on'); grid(ax3, 'on');
set(ax3, 'YScale', 'log'); 
xlabel(ax3, 'Total Control Points Added (k)');
ylabel(ax3, 'Root Mean Square Error (RMSE)');
title(ax3, sprintf('RMSE Convergence (N = %d)', N));

% Extract default color order to keep lines matching across subplots
colors = colororder; 

% Plot the error data
for i = 1:num_tests
    c_idx = mod(i-1, size(colors,1)) + 1; % Cycle colors if > 7 tests
    c = colors(c_idx, :);
    
    if ~isempty(all_iterations{i})
        plot(ax2, all_iterations{i}, all_max_err{i}, '.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', c);
        plot(ax3, all_iterations{i}, all_rmse{i}, '.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', c);
    end
end

legend(ax2, legend_labels, 'Location', 'northeast');
legend(ax3, legend_labels, 'Location', 'northeast');
linkaxes([ax2, ax3], 'xy'); % Keep error axes synced

hold(ax1, 'off'); hold(ax2, 'off'); hold(ax3, 'off');

%% SAVE PLOTS
% Ensure output directory exists

date_str = datestr(now, 'yyyymmdd');
save_filename = sprintf('Test2_2DErrorReduction_%s', date_str);

% Save as .fig for future MATLAB editing and .pdf for report insertion
savefig(fig_main, fullfile('06_results', [save_filename, '.fig']));
%exportgraphics(fig_main, fullfile('06_results', [save_filename, '.pdf']), 'ContentType', 'vector');

fprintf('Plot saved to /06_results/%s.fig\n', save_filename);