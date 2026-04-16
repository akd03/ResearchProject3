%% VARIABLE NAMING CONVENTION                                              
%  Nx       Known (Control) Points                                         
%  Nx_idx   Nx indices in Ax                                               
%  Ax       Aerofoil Mesh Points                                           
%  RAx      Applied Deformation                                            
%  DAx      Actual Aerofoil Position After Rx                              
%  Mx       Volume Mesh Points                                             

%% LOAD MESH                                                               
MeshFile = fopen("05_meshes\NACA0012257x129.xyz", "r");                    
Header = fscanf(MeshFile, '%d', 2);                                        
Ni = Header(1);
Nj = Header(2);
MeshData = textscan(MeshFile, '%f %f %f');
fclose(MeshFile);

% Set Aerofoil and Evaluation Points
Mx = [reshape(MeshData{1}, [Ni*Nj, 1]), reshape(MeshData{3}, [Ni*Nj, 1])];
Ax = [Mx(1:Ni,1), Mx(1:Ni,2)];

% Initial Control Points - LE and TE                                       
Nx_idx_base = [1; 129];
 
%% DEFORMATION
DAx = NACACamber(Ax, 0.02, 0.5);
RAx = DAx - Ax;

%% PARAMETER SWEEP: MAX ERROR THRESHOLD GREEDY
SF_R = 3;
max_err_threshold = 1e-6; % Suggested starting threshold for Max Error
N_IP_vals = [0, 5, 10, 15, 20, 30, 50]; 
M = size(Ax, 1); 

figure;
tiledlayout(1, 3, "TileSpacing", "compact");

% Setup Max Error Axis (Left - Now contains the threshold trigger)
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on');
set(ax1, 'YScale', 'log'); 
xlabel(ax1, 'Total Control Points Added (N)');
ylabel(ax1, 'Maximum Deformation Error');
title(ax1, 'Max Error Convergence');

% Draw horizontal threshold line on the Max Error plot
yline(ax1, max_err_threshold, 'k--', 'Target Tolerance', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'right');

% Setup RMSE Axis (Middle - Tracked for observation)
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
    
    % Initial Points (FPS)
    if N_IP > 0
        Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);
    end
    
    % Greedy Algorithm (Max Error Threshold Variant)
    [Nx_idx_final, max_err_history, rmse_history] = GreedyCholesky(Nx_idx, Ax, RAx, SF_R, "err_tol", max_err_threshold);
    
    N_G = length(max_err_history);
    N_Total = N_IP + N_G;
    
    iterations = (N_IP + 1) : N_Total; 
    
    % Calculate cumulative cost dynamically
    cost_full = zeros(1, N_Total);
    current_cost = 0;
    
    for k = 1:N_Total
        if k <= N_IP
            current_cost = current_cost + M;
            if k == N_IP
                current_cost = current_cost + (N_IP^3);
            end
        else
            current_cost = current_cost + (k^2 + (M - k) * k);
        end
        cost_full(k) = current_cost;
    end
    
    % Plot and capture line color
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

%% SAVINGS ANALYSIS 
figure;
hold on; grid on;

final_costs = zeros(1, length(N_IP_vals));
for i = 1:length(N_IP_vals)
    % Extract the final cost reached for this run
    lines = findobj(ax3, 'Type', 'Line');
    % Account for reverse stacking of lines in axes children
    line_idx = length(N_IP_vals) - i + 1; 
    y_data = get(lines(line_idx), 'YData');
    final_costs(i) = y_data(end); 
end

baseline_cost = final_costs(1);
pct_saved = ((baseline_cost - final_costs) ./ baseline_cost) * 100;

plot(N_IP_vals, pct_saved, 'k-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Number of Initial Points (N_{IP})');
ylabel('Total Operations Saved (%)');
title('Computational Savings vs. N_{IP} to Reach Max Error Threshold');
xticks(N_IP_vals);

hold off;

%% SMOOTH OPTIMIZATION CURVE: N_IP vs TOTAL OPERATIONS

max_err_threshold = 1e-3;
N_IP_vals = 0:1:50; % High-resolution sweep
total_operations = zeros(1, length(N_IP_vals));
final_N_totals = zeros(1, length(N_IP_vals));

% Run silently to gather the data
for i = 1:length(N_IP_vals)
    N_IP = N_IP_vals(i);
    Nx_idx = Nx_idx_base;
    
    if N_IP > 0
        Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);
    end
    
    % Run the greedy algorithm
    [~, max_err_history, ~] = GreedyCholesky(Nx_idx, Ax, RAx, SF_R, "err_tol", max_err_threshold);
    
    N_G = length(max_err_history);
    final_N_totals(i) = N_IP + N_G;
    
    % Interface with your Costing function to get the final total operations
    [~, ~, ~, ~, total_new] = Costing(M, N_IP, N_G);
    total_operations(i) = total_new;
end

%% PLOT THE SMOOTH CURVE
figure;
hold on; grid on;

% Plot using a solid line without markers for a smooth appearance
plot(N_IP_vals, total_operations, 'b-', 'LineWidth', 2);

% Formatting
xlabel('Number of Initial Points (N_{IP})');
ylabel('Total Operations Required');
title(sprintf('Computational Cost to Achieve Max Error < %g', max_err_threshold));
set(gca, 'YScale', 'log'); % Often better for viewing large operational differences

% Highlight the optimal point (the minimum cost)
[min_cost, min_idx] = min(total_operations);
optimal_N_IP = N_IP_vals(min_idx);

plot(optimal_N_IP, min_cost, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
text(optimal_N_IP + 1, min_cost, sprintf('Optimal: N_{IP} = %d', optimal_N_IP), 'VerticalAlignment', 'bottom');

hold off;