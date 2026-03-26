%% VARIABLE NAMING CONVENTION
%  Nx       Known (Control) Points
%  Nx_idx   Nx indices in Ax
%  Ax       Aerofoil Mesh Points
%  RAx      Applied Deformation 
%  DAx      Actual Aerofoil Position After Rx
%  Mx       Volume Mesh Points

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
%DAx = Rotate(Ax, 40, 0.5);
%DAx = Translate(DAx, [0, 0]);
DAx = NACACamber(Ax, 0.02, 0.5);
RAx = DAx - Ax;

%% PARAMETER SWEEP & GREEDY ALGORITHM
SF_R = 3;
N_vals = [255];
pct_vals = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
M = size(Ax, 1); % Total number of candidate nodes for complexity math

for N = N_vals
    figure;
    % Expand layout to 3 tiles
    tiledlayout(1, 3, "TileSpacing", "compact");
    
    % Setup Max Error Axis (Left Tile)
    ax1 = nexttile;
    hold(ax1, 'on'); grid(ax1, 'on');
    set(ax1, 'YScale', 'log'); 
    xlabel(ax1, 'Total Control Points Added (k)');
    ylabel(ax1, 'Maximum Deformation Error');
    title(ax1, sprintf('Max Error Convergence (N = %d)', N));
    
    % Setup RMSE Axis (Middle Tile)
    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on');
    set(ax2, 'YScale', 'log'); 
    xlabel(ax2, 'Total Control Points Added (k)');
    ylabel(ax2, 'Root Mean Square Error (RMSE)');
    title(ax2, sprintf('RMSE Convergence (N = %d)', N));

    % Setup Computational Cost Axis (Right Tile)
    ax3 = nexttile;
    hold(ax3, 'on'); grid(ax3, 'on');
    set(ax3, 'YScale', 'log'); 
    xlabel(ax3, 'Total Control Points Added (k)');
    ylabel(ax3, 'Cumulative Operations');
    title(ax3, sprintf('Theoretical Cost (N = %d)', N));
    
    legend_labels = strings(length(pct_vals), 1);
    
    for i = 1:length(pct_vals)
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
            [Nx_idx_final, max_err_history, rmse_history] = Greedy(Nx_idx, Ax, RAx, SF_R, "K", N_G);
            
            % Shift the x-axis for error plots to start exactly at N_IP
            iterations = N_IP + (1:length(max_err_history)) - 1; 
            
            % Calculate cumulative cost for the entire 1 to N sequence
            k_full = 1:N;
            cost_full = zeros(1, N);
            current_cost = 0;
            
            for k = 1:N
                if k <= N_IP
                    % Linear FPS Cost
                    current_cost = current_cost + M;
                else
                    % Quartic Greedy Cost (Matrix Solve + Field Eval)
                    current_cost = current_cost + (k^3 + M * k);
                end
                cost_full(k) = current_cost;
            end
            
            % Plot to ax1 and extract the automatically assigned color
            p1 = plot(ax1, iterations, max_err_history, '.-', 'LineWidth', 1.5, 'MarkerSize', 12);
            line_color = p1.Color; 
            
            % Force ax2 and ax3 to use the exact same color
            plot(ax2, iterations, rmse_history, '.-', 'LineWidth', 1.5, 'MarkerSize', 12, 'Color', line_color);
            plot(ax3, k_full, cost_full, '-', 'LineWidth', 2, 'Color', line_color);
        end
        
        % Format legend text
        legend_labels(i) = sprintf('N_{IP} = %d (%d%%)', N_IP, round(pct*100));
    end
    
    % Apply legends to all three tiles
    legend(ax1, legend_labels, 'Location', 'northeast');
    legend(ax2, legend_labels, 'Location', 'northeast');
    legend(ax3, legend_labels, 'Location', 'northwest');
    
    hold(ax1, 'off');
    hold(ax2, 'off');
    hold(ax3, 'off');
end



%% ANIMATION OF POINT ADDITION
%num_base = length(Nx_idx_base);
%PlotPointSeq(Ax, Nx_idx_final, num_base, N_IP, 0.1);