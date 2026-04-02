%% THEORETICAL TIME COMPLEXITY ANALYSIS (PARAMETER SWEEP)

% Fixed candidate mesh size
M = 100000;     

% Array of N values to test (e.g., 10%, 50%, and 100% of M)
N_vals = [10000, 50000, 100000];

for j = 1:length(N_vals)
    N = N_vals(j);
    
    % Generate 100 evenly spaced evaluation points to optimize script runtime
    step_size = max(1, round(N / 100));
    N_IP_vals = 0:step_size:N;
    
    % Ensure exactly N is included at the end of the array
    if N_IP_vals(end) ~= N
        N_IP_vals = [N_IP_vals, N];
    end
    
    % Preallocate cost arrays
    C_fps = zeros(1, length(N_IP_vals));
    C_greedy_old = zeros(1, length(N_IP_vals));
    C_greedy_new = zeros(1, length(N_IP_vals));
    C_total_old = zeros(1, length(N_IP_vals));
    C_total_new = zeros(1, length(N_IP_vals));
    
    %% CALCULATE OPERATIONS
    for i = 1:length(N_IP_vals)
        n_ip = N_IP_vals(i);
        n_g = N - n_ip;
        
        % Call local Costing function
        [C_fps(i), C_greedy_old(i), C_greedy_new(i), C_total_old(i), C_total_new(i)] = Costing(M, n_ip, n_g);
    end

    %% PLOT COMPLEXITY TREND
    figure;
    hold on; grid on;
    
    % Plot Original Method (Dashed lines)
    plot(N_IP_vals, C_total_old, 'y-', 'LineWidth', 2, 'DisplayName', 'Total Cost (Original)');
    plot(N_IP_vals, C_greedy_old, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Greedy Cost (Original)');
    
    % Plot New Method (Solid lines)
    plot(N_IP_vals, C_total_new, 'r-', 'LineWidth', 1, 'DisplayName', 'Total Cost (Cholesky)');
    plot(N_IP_vals, C_greedy_new, 'm-', 'LineWidth', 1, 'DisplayName', 'Greedy Cost (Cholesky)');
    
    % Plot FPS Cost (Shared - Green line)
    plot(N_IP_vals, C_fps, 'g-', 'LineWidth', 1.5, 'DisplayName', 'FPS Cost');
    
    % Formatting
    xlabel('Number of Initial Points (N_{IP})');
    ylabel('Theoretical Operations (Cost)');
    title(sprintf('Cost Comparison (M = %d, N = %d)', M, N));
    legend('Location', 'southwest');
    set(gca, 'YScale', 'log');
    xlim([0, N]);
    hold off;

    %% PERFORMANCE COMPARISON & WORKLOAD DISTRIBUTION
    figure;
    tiledlayout(1, 2, "TileSpacing", "compact");
    
    % --- Tile 1: Absolute Difference in Operations ---
    ax1 = nexttile;
    hold(ax1, 'on'); grid(ax1, 'on');
    
    % Calculate and plot operations saved
    ops_difference = C_total_old - C_total_new;
    plot(ax1, N_IP_vals, ops_difference, 'k-', 'LineWidth', 2);
    
    xlabel(ax1, 'Number of Initial Points (N_{IP})');
    ylabel(ax1, 'Operations Saved (Original - Cholesky)');
    title(ax1, sprintf('Absolute Savings (N = %d)', N));
    set(ax1, 'YScale', 'log'); 
    xlim(ax1, [0, N]);
    hold(ax1, 'off');
    
    % --- Tile 2: Proportion of Operations (Line Plots) ---
    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on');
    
    % Calculate proportions (percentages) for the new method
    prop_fps = (C_fps ./ C_total_new) * 100;
    prop_greedy = (C_greedy_new ./ C_total_new) * 100;
    
    % Handle the NaN case at N_IP = N
    prop_fps(isnan(prop_fps)) = 100; 
    prop_greedy(isnan(prop_greedy)) = 0;
    
    % Plot as standard lines
    plot(ax2, N_IP_vals, prop_fps, 'g-', 'LineWidth', 2, 'DisplayName', 'FPS Phase');
    plot(ax2, N_IP_vals, prop_greedy, 'm-', 'LineWidth', 2, 'DisplayName', 'Greedy Phase');
    
    xlabel(ax2, 'Number of Initial Points (N_{IP})');
    ylabel(ax2, 'Proportion of Total Operations (%)');
    title(ax2, 'Workload Distribution (Cholesky)');
    legend(ax2, 'Location', 'west');
    xlim(ax2, [0, N]);
    ylim(ax2, [0, 100]);
    hold(ax2, 'off');
end
