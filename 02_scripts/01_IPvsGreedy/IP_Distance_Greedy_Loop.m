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
DAx = NACACamber(Ax, 0.02, 0.4);
RAx = DAx - Ax;

%% PARAMETER SWEEP & GREEDY ALGORITHM
SF_R = 3;
%N_vals = [30, 40, 50, 60, 100];
N_vals = [255];
pct_vals = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
%N_vals = [60];
%pct_vals = [0.5];

for N = N_vals
    figure;
    tiledlayout(1, 2, "TileSpacing", "compact");
    
    % Setup Max Error Axis (Left Tile)
    ax1 = nexttile;
    hold(ax1, 'on'); grid(ax1, 'on');
    set(ax1, 'YScale', 'log'); 
    xlabel(ax1, 'Total Control Points Added (N_{IP} + N_G)');
    ylabel(ax1, 'Maximum Deformation Error');
    title(ax1, sprintf('Max Error Convergence (N = %d)', N));

    % Setup RMSE Axis (Right Tile)
    ax2 = nexttile;
    hold(ax2, 'on'); grid(ax2, 'on');
    set(ax2, 'YScale', 'log'); 
    xlabel(ax2, 'Total Control Points Added (N_{IP} + N_G)');
    ylabel(ax2, 'Root Mean Square Error (RMSE)');
    title(ax2, sprintf('RMSE Convergence (N = %d)', N));
    
    legend_labels = strings(length(pct_vals), 1);
    
    for i = 1:length(pct_vals)
        pct = pct_vals(i);
        N_IP = round(N * pct);
        N_G = N - N_IP;
        
        % Reset the starting points
        Nx_idx = Nx_idx_base;
        
        % Initial Points
        if N_IP > 0
            Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);
        end
        
        % Greedy Algorithm
        if N_G > 0
            % Capture all three outputs from the updated Greedy function
            [Nx_idx_final, max_err_history, rmse_history] = Greedy(Nx_idx, Ax, RAx, SF_R, "K", N_G);
            
            % Shift the x-axis so the line starts exactly at N_IP
            iterations = N_IP + (1:length(max_err_history)) - 1; 
            
            % Plot to the respective axes
            plot(ax1, iterations, max_err_history, '.-', 'LineWidth', 1.5, 'MarkerSize', 12);
            plot(ax2, iterations, rmse_history, '.-', 'LineWidth', 1.5, 'MarkerSize', 12);
        end
        
        % Format legend text
        legend_labels(i) = sprintf('N_{IP} = %d (%d%%)', N_IP, round(pct*100));
    end
    
    % Apply legends to both tiles
    legend(ax1, legend_labels, 'Location', 'best');
    legend(ax2, legend_labels, 'Location', 'best');
    
    hold(ax1, 'off');
    hold(ax2, 'off');
end



%% ANIMATION OF POINT ADDITION
%num_base = length(Nx_idx_base);
%PlotPointSeq(Ax, Nx_idx_final, num_base, N_IP, 0.1);