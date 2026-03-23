%% THEORETICAL TIME COMPLEXITY ANALYSIS
% Fixed parameters
M = 1000000;        % Number of nodes in the candidate mesh
N = 1000;        % Fixed total number of control points (N_IP + N_G)

% Independent variable array (Number of Initial Points)
N_IP_vals = 0:N;

% Preallocate cost arrays
C_fps = zeros(1, length(N_IP_vals));
C_greedy = zeros(1, length(N_IP_vals));
C_total = zeros(1, length(N_IP_vals));

%% CALCULATE OPERATIONS
for i = 1:length(N_IP_vals)
    n_ip = N_IP_vals(i);
    n_g = N - n_ip;
    
    % 1. FPS Cost: O(M * N_IP)
    C_fps(i) = M * n_ip;
    
    % 2. Greedy Cost: Summation from k = (N_IP + 1) to N
    % Operations per iteration: k^3 (Matrix Solve) + M*k (Field Eval)
    if n_g > 0
        k = (n_ip + 1):N;
        C_greedy(i) = sum(k.^3 + M .* k);
    else
        C_greedy(i) = 0; % Zero cost if no greedy points are used
    end
    
    % Total Cost
    C_total(i) = C_fps(i) + C_greedy(i);
end

%% PLOT COMPLEXITY TREND
figure;
hold on; grid on;

% Plot individual components and the total cost
plot(N_IP_vals, C_fps, 'r-', 'LineWidth', 2, 'DisplayName', 'FPS Cost');
plot(N_IP_vals, C_greedy, 'b-', 'LineWidth', 2, 'DisplayName', 'Greedy Cost');
plot(N_IP_vals, C_total, 'k--', 'LineWidth', 2, 'DisplayName', 'Total Cost');

% Formatting
xlabel('Number of Initial Points (N_{IP})');
ylabel('Theoretical Operations (Cost)');
title(sprintf('Computational Cost vs. Proportion of Initial Points (Fixed N = %d)', N));
legend('Location', 'northeast');

% Use a log scale for the Y-axis due to the massive quartic scale of the Greedy method
set(gca, 'YScale', 'log');

% Set X-axis limits
xlim([0, N]);
hold off;