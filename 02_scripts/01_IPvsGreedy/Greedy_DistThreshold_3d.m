%% VARIABLE NAMING CONVENTION                                              
%  Nx       Known (Control) Points                                         
%  Nx_idx   Nx indices in Ax                                               
%  Ax       Aerofoil Mesh Points                                           
%  RAx      Applied Deformation                                            
%  DAx      Actual Aerofoil Position After Rx                              

%% LOAD 3D MESH       
%
clear; clc; 

filename = '\05_meshes\01_3d_meshes\surfacepoints140K.plt';
fid = fopen(filename, 'r');
headerLines = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ~isempty(regexp(line, '^\s*[-0-9.eE]+', 'once'))
        break; 
    end
    headerLines = headerLines + 1;
end
fclose(fid);

fprintf('Loading numeric data...\n');
RawData = readmatrix(filename, 'NumHeaderLines', headerLines, 'FileType', 'text');
Ax = RawData(:, 1:3); % 3D Coordinates [X, Y, Z]
M = size(Ax, 1);

% Base Initial Point - Seed with the Root Trailing Edge (Index 1)                                       
Nx_idx_base = [1];
 
%% DEFORMATION
DAx = Ax;
% Applying a dummy 3D twist along the span for testing purposes.
% Replace this with your actual 3D deformation mapping function.
DAx(:,3) = DAx(:,3) + 0.005 * Ax(:,2).^2; 
RAx = DAx - Ax;

%% PARAMETER SWEEP: MAX ERROR THRESHOLD GREEDY
SF_R = 100; % Updated support radius
max_err_threshold = 1e-5; % Target tolerance
N_IP_vals = [100]; 

figure;
tiledlayout(1, 3, "TileSpacing", "compact");

% Setup Max Error Axis (Left)
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on');
set(ax1, 'YScale', 'log'); 
xlabel(ax1, 'Total Control Points Added (N)');
ylabel(ax1, 'Maximum Deformation Error');
title(ax1, 'Max Error Convergence');

yline(ax1, max_err_threshold, 'k--', 'Target Tolerance', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');

% Setup RMSE Axis (Middle)
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');
set(ax2, 'YScale', 'log'); 
xlabel(ax2, 'Total Control Points Added (N)');
ylabel(ax2, 'Root Mean Square Error (RMSE)');
title(ax2, 'RMSE Tracking');

% Setup Computational Cost Axis (Right)
ax3 = nexttile;
hold(ax3, 'on'); grid(ax3, 'on');
set(ax3, 'YScale', 'log'); 
xlabel(ax3, 'Total Control Points Added (N)');
ylabel(ax3, 'Cumulative Operations');
title(ax3, 'Theoretical Cost');

legend_labels = strings(length(N_IP_vals), 1);

for i = 1:length(N_IP_vals)
    N_IP = N_IP_vals(i);
    Nx_idx = Nx_idx_base;
    
    if N_IP > 0
        Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);
    end
    
    [Nx_idx_final, max_err_history, rmse_history] = GreedyCholesky(Nx_idx, Ax, RAx, SF_R, "err_tol", max_err_threshold);
    
    N_G = length(max_err_history);
    N_Total = length(Nx_idx_base) + N_IP + N_G - 1; % Account for seed point
    
    iterations = (N_Total - N_G + 1) : N_Total; 
    
    cost_full = zeros(1, N_Total);
    current_cost = 0;
    
    for k = 1:N_Total
        if k <= (N_IP + length(Nx_idx_base) - 1)
            current_cost = current_cost + M;
            if k == (N_IP + length(Nx_idx_base) - 1)
                current_cost = current_cost + (k^3);
            end
        else
            current_cost = current_cost + (k^2 + (M - k) * k);
        end
        cost_full(k) = current_cost;
    end
    
    p1 = plot(ax1, iterations, max_err_history, '.-', 'LineWidth', 1.5, 'MarkerSize', 10);
    line_color = p1.Color; 
    
    plot(ax2, iterations, rmse_history, '.-', 'LineWidth', 1.5, 'MarkerSize', 10, 'Color', line_color);
    plot(ax3, 1:N_Total, cost_full, '-', 'LineWidth', 2, 'Color', line_color);
    
    legend_labels(i) = sprintf('N_{IP} = %d (Final N = %d)', N_IP, N_Total);
end

legend(ax1, legend_labels, 'Location', 'northeast');
legend(ax2, legend_labels, 'Location', 'northeast');
legend(ax3, legend_labels, 'Location', 'northwest');
linkaxes([ax1, ax2], 'x');
hold(ax1, 'off');
hold(ax2, 'off');
hold(ax3, 'off');

%% SMOOTH OPTIMIZATION CURVE: N_IP vs TOTAL OPERATIONS
%
disp('first loop complete');
N_IP_vals_smooth = 1; 
total_operations = zeros(1, length(N_IP_vals_smooth));
final_N_totals = zeros(1, length(N_IP_vals_smooth));

for i = 1:length(N_IP_vals_smooth)
    if mod(i, 10) == 0
        disp(i);
    end
    N_IP = N_IP_vals_smooth(i);
    Nx_idx = Nx_idx_base;
    
    if N_IP > 0
        Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);
    end
    
    [~, max_err_history, ~] = GreedyCholesky(Nx_idx, Ax, RAx, SF_R, "err_tol", max_err_threshold);
    
    N_G = length(max_err_history);
    final_N_totals(i) = N_IP + N_G;
    
    [~, ~, ~, ~, total_new] = Costing(M, N_IP, N_G);
    total_operations(i) = total_new;
end


% Find the optimal point (minimum computational cost)
[min_cost, min_idx] = min(total_operations);
optimal_N_IP = N_IP_vals_smooth(min_idx);
optimal_N_total = final_N_totals(min_idx);
%

%% PLOT THE SMOOTH CURVE (Stacked Operational and Node Efficiency)
figure;
tiledlayout(2, 1, "TileSpacing", "compact");

% --- TOP TILE: Number of Operations ---
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, N_IP_vals_smooth, total_operations, 'b-', 'LineWidth', 2);
ylabel(ax1, 'Total Operations Required');
title(ax1, sprintf('Optimization Trade-off to Achieve Max Error < %g', max_err_threshold));
set(ax1, 'YScale', 'log'); 

% Highlight the optimal point on the operations curve
plot(ax1, optimal_N_IP, min_cost, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(ax1, optimal_N_IP + 1, min_cost, sprintf('Optimal: N_{IP} = %d', optimal_N_IP), 'VerticalAlignment', 'bottom');
hold(ax1, 'off');

% --- BOTTOM TILE: Total Number of Nodes ---
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, N_IP_vals_smooth, final_N_totals, 'k-', 'LineWidth', 2);
xlabel(ax2, 'Number of Initial Points (N_{IP})');
ylabel(ax2, 'Total Final Control Points (N_{Total})');

% Highlight the resulting total node count at the optimal operational point
plot(ax2, optimal_N_IP, optimal_N_total, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(ax2, optimal_N_IP + 1, optimal_N_total, sprintf('N_{Total} = %d', optimal_N_total), 'VerticalAlignment', 'top');
hold(ax2, 'off');

% Link the X-axes so that if you zoom in on one, the other matches perfectly
linkaxes([ax1, ax2], 'x');
%

%% FINAL VISUALIZATION: X-Y PROJECTION OF CONTROL POINTS
fprintf('Generating X-Y projection for N_IP = %d...\n', optimal_N_IP);

% Evaluate point selection at the optimal value
Nx_idx_fps = Nx_idx_base;
Nx_idx_fps = IP_Distance(Nx_idx_fps, Ax, optimal_N_IP);

[Nx_idx_final, ~, ~] = GreedyCholesky(Nx_idx_fps, Ax, RAx, SF_R, "err_tol", max_err_threshold);

% Extract spatial coordinates for plotting
fps_points = Ax(Nx_idx_fps, :);

% Find greedy points by taking the difference between the final list and the FPS list
greedy_idx = setdiff(Nx_idx_final, Nx_idx_fps, 'stable');
greedy_points = Ax(greedy_idx, :);

figure;
hold on; grid on; axis equal;

% Plot the background surface mesh in light grey
scatter(Ax(:,1), Ax(:,2), 2, [0.8 0.8 0.8], 'filled', 'DisplayName', 'Surface Mesh');

% Plot the Initial (FPS) Points in Red
scatter(fps_points(:,1), fps_points(:,2), 10, 'r', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Initial Points (FPS)');

% Plot the Greedy Selected Points in Blue
if ~isempty(greedy_points)
    scatter(greedy_points(:,1), greedy_points(:,2), 10, 'b', 'filled', 'MarkerEdgeColor', 'k', 'DisplayName', 'Greedy Points');
end

xlabel('X Coordinate (Chordwise)');
ylabel('Y Coordinate (Spanwise)');
title(sprintf('Control Point Distribution (X-Y Projection) | N_{IP} = %d', optimal_N_IP));
legend('Location', 'northeastoutside');

hold off;

file_name = sprintf('GreedyOptimizationData_SF%d_Tol%g.mat', SF_R, max_err_threshold);

% Combine the directory and filename robustly
%
save_filepath = fullfile('06_results', file_name);

fprintf('Saving plotting variables to %s...\n', save_filepath);

% Save only the variables necessary to recreate the plots
save(save_filepath, ...
    'Ax', 'DAx', 'RAx', 'M', 'SF_R', 'max_err_threshold', ...
    'N_IP_vals_smooth', 'total_operations', 'final_N_totals', ...
    'optimal_N_IP', 'optimal_N_total', 'min_cost', ...
    'Nx_idx_base', 'Nx_idx_fps', 'Nx_idx_final', ...
    'fps_points', 'greedy_idx', 'greedy_points');

fprintf('Data successfully saved.\n');

PlotDeformation3D(Ax, DAx);