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
DAx = NACACamber(Ax, 0.09, 0.5);
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
all_cost = cell(num_tests, 1);
all_k_full = cell(num_tests, 1);
legend_labels = strings(num_tests, 1);

%% DATA GENERATION LOOP (Errors & RMSE)
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
fprintf('Data generation complete.\n');

%% INDEPENDENT COST CALCULATION
% Calculated entirely outside the main loop for easy logic modification
fprintf('Calculating theoretical operational costs...\n');
for i = 1:num_tests
    pct = pct_vals(i);
    N_IP = round(N * pct);
    
    all_k_full{i} = 1:N;
    all_cost{i} = CalculateCumulativeCost(M, N, N_IP);
end

%% PLOT 1: QUADRUPLE LAYOUT (SHAPE, MAX ERROR, RMSE, COST)
fprintf('Plotting...\n');
% Increased figure width to accommodate 4 plots cleanly
fig_main = figure('Name', '2D Error & Cost Reduction', 'Position', [100, 100, 1800, 450]);
tiledlayout(1, 4, "TileSpacing", "compact");

% --- TILE 1: Aerofoil Shape ---
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on'); axis(ax1, 'equal');
plot(ax1, Ax(:,1), Ax(:,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Original (NACA0012)');
plot(ax1, DAx(:,1), DAx(:,2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Deformed (NACA9512)');
xlabel(ax1, 'X');
ylabel(ax1, 'Y');
title(ax1, 'Aerofoil Boundary');
legend(ax1, 'Location', 'best');

% --- TILE 2: Max Error ---
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');
set(ax2, 'YScale', 'log'); 
xlabel(ax2, 'Total Control Points Added (k)');
ylabel(ax2, 'Maximum Deformation Error');
title(ax2, sprintf('Max Error Convergence'));

% --- TILE 3: RMSE ---
ax3 = nexttile;
hold(ax3, 'on'); grid(ax3, 'on');
set(ax3, 'YScale', 'log'); 
xlabel(ax3, 'Total Control Points Added (k)');
ylabel(ax3, 'Root Mean Square Error (RMSE)');
title(ax3, sprintf('RMSE Convergence'));

% --- TILE 4: Theoretical Cost ---
ax4 = nexttile;
hold(ax4, 'on'); grid(ax4, 'on');
set(ax4, 'YScale', 'log'); 
xlabel(ax4, 'Total Control Points Added (k)');
ylabel(ax4, 'Cumulative Operations');
title(ax4, sprintf('Theoretical Cost'));

% Extract default color order to keep lines matching across subplots
colors = colororder; 

% Plot the data across all tiles
for i = 1:num_tests
    c_idx = mod(i-1, size(colors,1)) + 1; % Cycle colors
    c = colors(c_idx, :);
    
    % Plot Errors (Starts at N_IP)
    if ~isempty(all_iterations{i})
        plot(ax2, all_iterations{i}, all_max_err{i}, '.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', c);
        plot(ax3, all_iterations{i}, all_rmse{i}, '.-', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', c);
    end
    
    % Plot Cost (Full range 1 to N)
    plot(ax4, all_k_full{i}, all_cost{i}, '-', 'LineWidth', 1.5, 'Color', c);
end

legend(ax2, legend_labels, 'Location', 'northeast');

% Keep error axes synced (Cost is excluded as magnitudes differ)
linkaxes([ax2, ax3], 'xy'); 

hold(ax1, 'off'); hold(ax2, 'off'); hold(ax3, 'off'); hold(ax4, 'off');

%% SAVE PLOTS


date_str = datestr(now, 'yyyymmdd');
save_filename = sprintf('Test2_2DErrorAndCost_%s', date_str);

savefig(fig_main, fullfile('06_results', [save_filename, '.fig']));
%exportgraphics(fig_main, fullfile('06_results', [save_filename, '.pdf']), 'ContentType', 'vector');

fprintf('Plot saved to /06_results/%s.fig\n', save_filename);

%% LOCAL FUNCTIONS
function cost_full = CalculateCumulativeCost(M, N, N_IP)
    % Calculates the theoretical cumulative operations at each step k
    cost_full = zeros(1, N);
    current_cost = 0;
    
    for k = 1:N
        if k <= N_IP
            % Linear FPS Cost
            current_cost = current_cost + M;
            % Add the one-off Cholesky setup cost once FPS is finished
            if k == N_IP
                current_cost = current_cost + (N_IP^3);
            end
        else
            % Cubic Greedy Cost (Cholesky Update/Solve + Field Eval)
            current_cost = current_cost + (k^2 + (M - k) * k);
        end
        cost_full(k) = current_cost;
    end
end