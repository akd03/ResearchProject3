%% VARIABLE NAMING CONVENTION
%  Nx       Known (Control) Points
%  Nx_idx   Nx indices in Ax
%  Ax       Aerofoil Mesh Points
%  RAx      Applied Deformation 
%  DAx      Actual Aerofoil Position After Rx
%  Mx       Volume Mesh Points

%% LOAD MESH
%
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
%}


 

%% DEFORMATION
% True Deformation
DAx = Rotate(Ax, 30, 0);
DAx = Translate(DAx, [0, 1]);
RAx = DAx - Ax;

%% PARAMETER SWEEP & GREEDY ALGORITHM
SF_R = 3;
%N_vals = [10, 20, 30, 40, 50, 60, 100];
%pct_vals = [0, 0.20, 0.40, 0.60, 0.80];

N_vals = [60];
pct_vals = [0.5];

for N = N_vals
    figure;
    hold on; grid on;
    set(gca, 'YScale', 'log'); % Log scale for better error visualization
    
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
            [Nx_idx_final, max_err_history] = Greedy(Nx_idx, Ax, RAx, SF_R, "K", N_G);
            
            % 3. Plotting the Error History
            % Shift the x-axis so the line starts exactly at N_IP
            % The -1 ensures the first plotted point corresponds to the initial error
            iterations = N_IP + (1:length(max_err_history)) - 1; 
            plot(iterations, max_err_history, '.-', 'LineWidth', 1.5, 'MarkerSize', 12);
        end
        
        % Format legend text
        legend_labels(i) = sprintf('N_{IP} = %d (%d%%)', N_IP, round(pct*100));
    end
    
    % Finalize Figure Formatting
    xlabel('Total Control Points Added (N_{IP} + N_G)');
    ylabel('Maximum Deformation Error');
    title(sprintf('Greedy Convergence for N = %d', N));
    legend(legend_labels, 'Location', 'best');
    hold off;
end



%% ANIMATION OF POINT ADDITION
figure;
hold on; 
axis equal;
title('Control Point Selection Sequence');
xlabel('X'); 
ylabel('Y');

% Plot the base aerofoil mesh in light gray for context
plot(Ax(:,1), Ax(:,2), '.-', 'Color', [0.7 0.7 0.7]); 

% Configurable delay between points in seconds
anim_delay = 0.2; 

num_base = length(Nx_idx_base);
num_total = length(Nx_idx_final);

% Iterate through the chronological list of indices
for k = 1:num_total
    
    % Get the geometric coordinates for the current chronological point
    idx = Nx_idx_final(k);
    pt_x = Ax(idx, 1);
    pt_y = Ax(idx, 2);
    
    % Determine color based on when it was added to the list
    if k <= num_base
        pt_color = 'g'; % First 2 points: Base (LE, TE)
    elseif k <= num_base + N_IP
        pt_color = 'r'; % Next N_IP points: Initial Selection
    else
        pt_color = 'b'; % Remaining points: Greedy Selection
    end
    
    % Plot the new point as a black 'x'
    h_new = plot(pt_x, pt_y, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    drawnow;
    
    % Pause for visibility
    pause(anim_delay);
    
    % Replace the 'x' with the permanent colored 'o'
    delete(h_new);
    plot(pt_x, pt_y, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', pt_color, 'MarkerSize', 8);
    drawnow;
end
hold off;