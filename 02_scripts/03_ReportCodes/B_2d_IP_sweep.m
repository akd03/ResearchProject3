%% VARIABLE NAMING CONVENTION                                              
%  Nx       Known (Control) Points                                         
%  Nx_idx   Nx indices in Ax                                               
%  Ax       Aerofoil Mesh Points                                           
%  RAx      Applied Deformation                                            
%  DAx      Actual Aerofoil Position After Rx                              
%  Mx       Volume Mesh Points                                             

clear; clc; close all;

%% LOAD MESH                                                               
MeshFile = fopen("05_meshes\NACA0012257x129.xyz", "r");                    
Header = fscanf(MeshFile, '%d', 2);                                        
Ni = Header(1);
Nj = Header(2);
MeshData = textscan(MeshFile, '%f %f %f');
fclose(MeshFile);

Mx = [reshape(MeshData{1}, [Ni*Nj, 1]), reshape(MeshData{3}, [Ni*Nj, 1])];
Ax = [Mx(1:Ni,1), Mx(1:Ni,2)];
Nx_idx_base = [1; 129];
 
%% DEFORMATION
DAx = NACACamber(Ax, 0.09, 0.5);
RAx = DAx - Ax;

%% PARAMETER SWEEP SETUP
SF_R = 3;
M = size(Ax, 1); 

% Setup for Sweep 1 (Fixed N, vary N_IP)
N_vals = [25, 50, 100];
num_N_tests = length(N_vals);
s1_N_IP = cell(num_N_tests, 1);
s1_final_err = cell(num_N_tests, 1);
s1_total_cost = cell(num_N_tests, 1);

% Setup for Sweep 2 (Fixed Tolerance, vary N_IP)
err_threshold = 1e-6;
s2_N_IP_vals = 0:1:100; % Sweeping from 0 to 100 in steps of 5
num_tol_tests = length(s2_N_IP_vals);
s2_total_nodes = zeros(num_tol_tests, 1);
s2_total_cost = zeros(num_tol_tests, 1);

%% DATA GENERATION: SWEEP 1 (Fixed N)
fprintf('Running Sweep 1: Fixed N, Varying N_IP...\n');
for i = 1:num_N_tests
    N = N_vals(i);
    N_IP_range = 0:(N - 1);
    
    current_N_IP = zeros(length(N_IP_range), 1);
    current_err = zeros(length(N_IP_range), 1);
    current_cost = zeros(length(N_IP_range), 1);
    
    for j = 1:length(N_IP_range)
        N_IP = N_IP_range(j);
        N_G = N - N_IP;
        Nx_idx = Nx_idx_base;
        
        if N_IP > 0
            Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);
        end
        
        if N_G > 0
            [~, max_err_hist, ~] = GreedyCholesky(Nx_idx, Ax, RAx, SF_R, 'K', N_G);
            current_err(j) = max_err_hist(end);
        else
            % If N_G is 0, we need to evaluate the error of just the IP points manually
            % Assuming GreedyCholesky handles K=0 by evaluating initial error, 
            % otherwise we set to NaN or manually evaluate.
            [~, max_err_hist, ~] = GreedyCholesky(Nx_idx, Ax, RAx, SF_R, 'K', 1);
            current_err(j) = max_err_hist(1); 
        end
        
        current_N_IP(j) = N_IP;
        current_cost(j) = CalculateTotalCost(M, N, N_IP);
    end
    
    s1_N_IP{i} = current_N_IP;
    s1_final_err{i} = current_err;
    s1_total_cost{i} = current_cost;
end

%% DATA GENERATION: SWEEP 2 (Fixed Error Tolerance)
fprintf('Running Sweep 2: Fixed Error Tolerance, Varying N_IP...\n');
for i = 1:num_tol_tests
    N_IP = s2_N_IP_vals(i);
    Nx_idx = Nx_idx_base;
    
    if N_IP > 0
        Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);
    end
    
    % Run until error tolerance is met (K = Inf)
    [Nx_idx_final, max_err_hist, ~] = GreedyCholesky(Nx_idx, Ax, RAx, SF_R, 'err_tol', err_threshold, 'K', Inf);
    
    % The number of greedy points added is the length of the error history
    N_G = length(max_err_hist);
    N_total = N_IP + N_G;
    
    s2_total_nodes(i) = N_total;
    s2_total_cost(i) = CalculateTotalCost(M, N_total, N_IP);
end
fprintf('Data generation complete.\n');

%% PLOTTING SWEEP 1: FIXED N
fprintf('Plotting Sweep 1...\n');
fig_sweep1 = figure('Name', 'Sweep 1: Fixed Target Nodes', 'Position', [100, 100, 1000, 450]);
tiledlayout(1, 2, "TileSpacing", "compact");

% --- TILE 1: Final Max Error ---
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on');
set(ax1, 'YScale', 'log');
xlabel(ax1, 'Number of Initial Points (N_{IP})');
ylabel(ax1, 'Final Maximum Deformation Error');
title(ax1, 'Effect of Initial Points on Final Error');

% --- TILE 2: Total Operations ---
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');
set(ax2, 'YScale', 'log');
xlabel(ax2, 'Number of Initial Points (N_{IP})');
ylabel(ax2, 'Total Operations (FLOPs)');
title(ax2, 'Theoretical Computational Cost');

colors = colororder;
legend_labels_s1 = strings(num_N_tests, 1);

for i = 1:num_N_tests
    c_idx = mod(i-1, size(colors,1)) + 1; 
    c = colors(c_idx, :);
    
    plot(ax1, s1_N_IP{i}, s1_final_err{i}, '.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', c);
    plot(ax2, s1_N_IP{i}, s1_total_cost{i}, '.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', c);
    legend_labels_s1(i) = sprintf('Total N = %d', N_vals(i));
end

legend(ax1, legend_labels_s1, 'Location', 'best');
legend(ax2, legend_labels_s1, 'Location', 'best');
hold(ax1, 'off'); hold(ax2, 'off');

%% PLOTTING SWEEP 2: FIXED ERROR TOLERANCE
fprintf('Plotting Sweep 2...\n');
fig_sweep2 = figure('Name', 'Sweep 2: Fixed Error Tolerance', 'Position', [150, 150, 1000, 450]);
tiledlayout(1, 2, "TileSpacing", "compact");

% --- TILE 1: Required Control Nodes ---
ax3 = nexttile;
hold(ax3, 'on'); grid(ax3, 'on');
xlabel(ax3, 'Number of Initial Points (N_{IP})');
ylabel(ax3, 'Total Control Nodes Required (N_{total})');
title(ax3, sprintf('Nodes Required to Reach E_{tol} = %g', err_threshold));

% --- TILE 2: Total Operations ---
ax4 = nexttile;
hold(ax4, 'on'); grid(ax4, 'on');
set(ax4, 'YScale', 'log');
xlabel(ax4, 'Number of Initial Points (N_{IP})');
ylabel(ax4, 'Total Operations (FLOPs)');
title(ax4, 'Theoretical Computational Cost');

plot(ax3, s2_N_IP_vals, s2_total_nodes, 'k.-', 'LineWidth', 1.5, 'MarkerSize', 8);
plot(ax4, s2_N_IP_vals, s2_total_cost, 'k.-', 'LineWidth', 1.5, 'MarkerSize', 8);

hold(ax3, 'off'); hold(ax4, 'off');

%% SAVE PLOTS
if ~exist('06_results', 'dir')
    mkdir('06_results');
end

date_str = datestr(now, 'yyyymmdd');
save_filename_s1 = sprintf('Test3_SweepFixedN_%s', date_str);
save_filename_s2 = sprintf('Test4_SweepFixedTol_%s', date_str);

savefig(fig_sweep1, fullfile('06_results', [save_filename_s1, '.fig']));
savefig(fig_sweep2, fullfile('06_results', [save_filename_s2, '.fig']));

fprintf('Plots saved to /06_results/\n');

%% LOCAL FUNCTIONS
function total_cost = CalculateTotalCost(M, N, N_IP)
    % Calculates the total integrated cost at the final iteration
    total_cost = 0;
    
    % 1. FPS Cost
    total_cost = total_cost + (M * N_IP);
    
    % 2. Cholesky Setup Cost (Factoring the initial N_IP matrix)
    if N_IP > 0
        total_cost = total_cost + (N_IP^3);
    end
    
    % 3. Greedy Recurrence Loop Cost
    N_G = N - N_IP;
    if N_G > 0
        k = (N_IP + 1):N;
        total_cost = total_cost + sum(k.^2 + (M - k) .* k);
    end
end