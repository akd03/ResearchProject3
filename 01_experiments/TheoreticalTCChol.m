%% THEORETICAL TIME COMPLEXITY ANALYSIS
% Fixed parameters
M = 1000000;     % Number of nodes in the candidate mesh
N = 1000;        % Fixed total number of control points (N_IP + N_G)

% Independent variable array (Number of Initial Points)
N_IP_vals = 0:N;

% Preallocate cost arrays
C_fps = zeros(1, length(N_IP_vals));
C_greedy_old = zeros(1, length(N_IP_vals));
C_greedy_new = zeros(1, length(N_IP_vals));

%% CALCULATE OPERATIONS
for i = 1:length(N_IP_vals)
    n_ip = N_IP_vals(i);
    n_g = N - n_ip;
    
    % 1. FPS Cost: O(M * N_IP)
    C_fps(i) = M * n_ip;
    
    % 2. Greedy Costs
    if n_g > 0
        k = (n_ip + 1):N;
        
        % Original Method: k^3 (Matrix Solve) + M*k (Full Field Eval)
        C_greedy_old(i) = sum(k.^3 + M .* k);
        
        % New Method: k^2 (Cholesky Update) + (M-k)*k (Subset Field Eval)
        C_greedy_new(i) = sum(k.^2 + (M - k) .* k);
    else
        C_greedy_old(i) = 0; 
        C_greedy_new(i) = 0; 
    end
end

% Total Costs
C_total_old = C_fps + C_greedy_old;
C_total_new = C_fps + C_greedy_new;

%% PLOT COMPLEXITY TREND
figure;
hold on; grid on;

% Plot Original Method (Dashed lines)
plot(N_IP_vals, C_total_old, 'y-', 'LineWidth', 2, 'DisplayName', 'Total Cost (Original)');
plot(N_IP_vals, C_greedy_old, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Greedy Cost (Original)');

% Plot New Method (Solid lines)
plot(N_IP_vals, C_total_new, 'r-', 'LineWidth', 1, 'DisplayName', 'Total Cost (Cholesky)');
plot(N_IP_vals, C_greedy_new, 'm-', 'LineWidth', 1, 'DisplayName', 'Greedy Cost (Cholesky)');

% Plot FPS Cost (Shared - Red line)
plot(N_IP_vals, C_fps, 'g-', 'LineWidth', 1.5, 'DisplayName', 'FPS Cost');

% Formatting
xlabel('Number of Initial Points (N_{IP})');
ylabel('Theoretical Operations (Cost)');
title(sprintf('Computational Cost Comparison (M = %d, N = %d)', M, N));
legend('Location', 'southwest');

% Use a log scale for the Y-axis
set(gca, 'YScale', 'log');
xlim([0, N]);
hold off;

%% PERFORMANCE COMPARISON & WORKLOAD DISTRIBUTION
figure;
tiledlayout(1, 2, "TileSpacing", "compact");

% --- Tile 1: Absolute Difference in Operations ---
ax1 = nexttile;
hold(ax1, 'on'); grid(ax1, 'on');

% Calculate operations saved
ops_difference = C_total_old - C_total_new;

% Plot the absolute difference
plot(ax1, N_IP_vals, ops_difference, 'k-', 'LineWidth', 2);

xlabel(ax1, 'Number of Initial Points (N_{IP})');
ylabel(ax1, 'Operations Saved (Original - Cholesky)');
title(ax1, 'Absolute Computational Savings');

% Use a log scale because the savings drop from billions to zero
set(ax1, 'YScale', 'log'); 
xlim(ax1, [0, N]);
hold(ax1, 'off');

% --- Tile 2: Proportion of Operations (Cholesky Method) ---
ax2 = nexttile;
hold(ax2, 'on'); grid(ax2, 'on');

% Calculate proportions (percentages) for the new method
prop_fps = (C_fps ./ C_total_new) * 100;
prop_greedy = (C_greedy_new ./ C_total_new) * 100;

% Handle the NaN case at N_IP = N where greedy cost is 0, but division by 0 might occur if total was 0
prop_fps(isnan(prop_fps)) = 100; 
prop_greedy(isnan(prop_greedy)) = 0;

% Create a stacked area chart to show the 100% distribution
area_data = [prop_fps', prop_greedy'];
a = area(ax2, N_IP_vals, area_data, 'LineStyle', 'none');

% Match your requested color scheme aesthetics 
a(1).FaceColor = 'g'; % FPS Cost 
a(2).FaceColor = 'm'; % Greedy Cost 
a(1).FaceAlpha = 0.6;
a(2).FaceAlpha = 0.6;

xlabel(ax2, 'Number of Initial Points (N_{IP})');
ylabel(ax2, 'Proportion of Total Operations (%)');
title(ax2, 'Workload Distribution (Cholesky Method)');
legend(ax2, 'FPS Phase', 'Greedy Phase', 'Location', 'southwest');

xlim(ax2, [0, N]);
ylim(ax2, [0, 100]);
hold(ax2, 'off');