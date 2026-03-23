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
Nx_idx = [1, 129];
%}
%% CARTESIAN MESH GENERATION
%{
x_min = -5;
x_max = 5;
y_min = -5;
y_max = 5;
AerofoilData = readmatrix("05_meshes\RAE2822.csv");
x_coords = AerofoilData(:, 1);
y_upper  = AerofoilData(:, 2);                  
y_lower  = AerofoilData(:, 3);
Ax_up = [flipud(x_coords), flipud(y_upper)];
Ax_lo = [x_coords, y_lower];
Ax = [Ax_up; Ax_lo];

Mx_x = x_min:0.1:x_max;
Mx_y = y_min:0.1:y_max;
[Mx_X, Mx_Y] = meshgrid(Mx_x, Mx_y);
Mx = [Mx_X(:), Mx_Y(:)];
Nx_idx = [1, 65];
%}

%% DEFORMATION
% True Deformation
%DAx = Rotate(Ax, 30, 0);
%DAx = Translate(DAx, [0, 1]);
DAx = NACACamber(Ax, 0.05, 0.4);
RAx = DAx - Ax;

%% INITIAL POINT SELECTION
N_IP = 0;
Nx_idx = IP_Distance(Nx_idx, Ax, N_IP);

%% GREEDY ALGORITHM
N_G = 20;
SF_R = 10;
[Nx_idx, max_err_history] = Greedy(Nx_idx, Ax, RAx, SF_R, "K", N_G);

%% FINAL MESH DEFORMATION
Nx = Ax(Nx_idx, :); 
RNx = RAx(Nx_idx, :);
% Solve for weights
CNx = Phi_WC2(Norm(Nx, Nx) / SF_R);              
G_x = CNx \ RNx(:, 1);             
G_y = CNx \ RNx(:, 2);             
% Interpolate Deformation
CMx = Phi_WC2(Norm(Mx, Nx) / SF_R);
SMx_x = CMx * G_x;                  
SMx_y = CMx * G_y;           
DMx = [Mx(:,1)+SMx_x, Mx(:,2)+SMx_y];

%% RESULTS/PLOTTING
%% RESULTS/PLOTTING
figure; 
tiledlayout(2, 2, "TileSpacing", "tight"); 

% Top Left: Initial Geometry and Control Points
nexttile(1); hold on; axis equal;
%plot(Mx(:,1), Mx(:,2), "b.");
plot(Ax(:,1), Ax(:,2), "r.-");
plot(Ax(Nx_idx(1:2), 1), Ax(Nx_idx(1:2),2), "go", "MarkerFaceColor", "g", "MarkerSize", 8);
plot(Ax(Nx_idx(3:2+N_IP), 1), Ax(Nx_idx(3:2+N_IP),2), "bo", "MarkerFaceColor", "b", "MarkerSize", 8);
plot(Ax(Nx_idx(3+N_IP:end), 1), Ax(Nx_idx(3+N_IP:end),2), "ro", "MarkerFaceColor", "r", "MarkerSize", 8);
title('Initial Aerofoil');

% Bottom Left: Deformed Geometry and Control Points
nexttile(3); hold on; axis equal;
xline(0, "k");
yline(0, "k");
plot(DAx(:,1), DAx(:,2), "b-"); 
%plot(DMx(:,1), DMx(:,2), "r.");
%plot(DAx(Nx_idx,1), DAx(Nx_idx,2), "ro", "MarkerFaceColor", "r", "MarkerSize", 8);
title('Deformed Aerofoil');

% Right (Spanning two rows): Max Error History
nexttile(2, [2, 1]); hold on; grid on;
% Plot the error history starting from iteration 1
iterations = 1:length(max_err_history);
plot(iterations, max_err_history, "k.-", "LineWidth", 1.5, "MarkerSize", 12);
xlabel('Greedy Points Added');
ylabel('Maximum Deformation Error');
title('Greedy Algorithm Convergence');

%PlotOMesh(Mx, DMx, Ni, Nj);