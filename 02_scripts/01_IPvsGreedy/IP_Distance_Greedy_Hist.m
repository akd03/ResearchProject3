%% PARAMETER SWEEP & GREEDY ALGORITHM
SF_R = 3;
N_vals = [80];
pct_vals = [0, 0.1, 0.20, 0.3, 0.40, 0.5, 0.60, 0.7, 0.80, 0.9];
num_base = length(Nx_idx_base); % Store the number of base points (2)

for N = N_vals
    
    % --- SETUP ERROR CONVERGENCE FIGURE ---
    figure;
    tiledlayout(1, 2, "TileSpacing", "compact");
    
    ax1 = nexttile;
    hold(ax1, 'on'); grid(ax1, 'on');
    set(ax1, 'YScale', 'log'); 
    xlabel(ax1, 'Total Control Points Added (N_{IP} + N_G)');
    ylabel(ax1, 'Maximum Deformation Error');
    title(ax1, sprintf('Max Error Convergence (N = %d)', N));
    
    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on');
    set(ax2, 'YScale', 'log'); 
    xlabel(ax2, 'Total Control Points Added (N_{IP} + N_G)');
    ylabel(ax2, 'Root Mean Square Error (RMSE)');
    title(ax2, sprintf('RMSE Convergence (N = %d)', N));
    
    legend_labels = strings(length(pct_vals), 1);
    
    % --> INITIALIZE SEPARATE TRACKING ARRAYS
    all_IP_indices = []; 
    all_G_indices = [];
    
    for i = 1:length(pct_vals)
        pct = pct_vals(i);
        N_IP = round(N * pct);
        N_G = N - N_IP;
        
        Nx_idx = Nx_idx_base;
        
        if N_IP > 0
            Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);
        end
        
        if N_G > 0
            [Nx_idx_final, max_err_history, rmse_history] = Greedy(Nx_idx, Ax, RAx, SF_R, "K", N_G);
            
            % --> SLICE AND STORE THE INDICES
            % Initial points start immediately after the base points
            idx_IP = Nx_idx_final(num_base + 1 : num_base + N_IP);
            % Greedy points start immediately after the initial points
            idx_G  = Nx_idx_final(num_base + N_IP + 1 : end);
            
            all_IP_indices = [all_IP_indices; idx_IP];
            all_G_indices = [all_G_indices; idx_G];
            
            iterations = N_IP + (1:length(max_err_history)) - 1; 
            plot(ax1, iterations, max_err_history, '.-', 'LineWidth', 1.5, 'MarkerSize', 12);
            plot(ax2, iterations, rmse_history, '.-', 'LineWidth', 1.5, 'MarkerSize', 12);
        end
        
        legend_labels(i) = sprintf('N_{IP} = %d (%d%%)', N_IP, round(pct*100));
    end
    
    legend(ax1, legend_labels, 'Location', 'best');
    legend(ax2, legend_labels, 'Location', 'best');
    hold(ax1, 'off');
    hold(ax2, 'off');
    
    % --- GENERATE STACKED SPATIAL DISTRIBUTION BAR CHART ---
    figure;
    hold on; grid on;

    % Create bin edges so exactly 1 bin exists for each integer index (1 to length(Ax))
    edges = 0.5 : 1 : length(Ax) + 0.5;
    
    % Count occurrences for exactly 257 bins
    count_IP = histcounts(all_IP_indices, edges);
    count_G = histcounts(all_G_indices, edges);
    
    % Create the stacked bar chart. Width '1' removes gaps between bars.
    b = bar(1:length(Ax), [count_IP', count_G'], 1, 'stacked', 'EdgeColor', 'none');
    
    % Apply strict color coding
    b(1).FaceColor = 'r'; % Initial Points in Red
    b(2).FaceColor = 'b'; % Greedy Points in Blue

    % Add vertical reference lines for the geometric landmarks
    xline(1, 'k-', 'Lower TE', 'LabelVerticalAlignment', 'top', 'LineWidth', 1.5);
    xline(129, 'k-', 'Leading Edge', 'LabelVerticalAlignment', 'top', 'LineWidth', 1.5);
    xline(length(Ax), 'k-', 'Upper TE', 'LabelVerticalAlignment', 'top', 'LineWidth', 1.5);

    % Format the plot
    xlabel('Aerofoil Index');
    ylabel('Frequency');
    xlim([1, length(Ax)]);
    legend('Initial Points (Distance)', 'Greedy Points', 'Location', 'northwest');

    % Make the y-axis tighter
    ylim_val = max(count_IP + count_G) * 1.1;
    if ylim_val > 0
        ylim([0, ylim_val]);
    end

    hold off;
    
end