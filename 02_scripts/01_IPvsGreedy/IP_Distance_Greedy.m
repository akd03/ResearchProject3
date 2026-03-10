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
Nx_idx = [1, 129];



%% DEFORMATION
% True Deformation
DAx_1 = Rotate(Ax, 120, 0.5);
DAx = Translate(DAx_1, [0, 0]);
RAx = DAx - Ax;

%% INITIAL POINT SELECTION
Nx_idx = IP_Distance(Nx_idx, Ax, 15);

%% GREEDY ALGORITHM
SF_R = 10;
[Nx_idx, max_err_history] = Greedy(Nx_idx, Ax, RAx, SF_R, "K", 20);
Nx = Ax(Nx_idx, :); 
RNx = RAx(Nx_idx, :);

% 2. Build the Square Training Matrix (Nc x Nc)
% This determines how the control points interact with each other
CNx = Phi_WC2(Norm(Nx, Nx) / SF_R);              

% 3. Solve for the 30 Weights (Nc x 1)
G_x = CNx \ RNx(:, 1);             
G_y = CNx \ RNx(:, 2);             

% 4. Build the Rectangular Evaluation Matrix (Nm x Nc)
% Notice the order: Norm(Mx, Nx) yields a 30000 x 30 matrix
CMx = Phi_WC2(Norm(Mx, Nx) / SF_R);

% 5. Interpolate the deformation for the entire mesh (Nm x 1)
RMx_x = CMx * G_x;                  
RMx_y = CMx * G_y;           

% Final Deformation Array
DMx = [Mx(:,1)+RMx_x, Mx(:,2)+RMx_y];

%% RESULTS/PLOTTING
figure; 
tiledlayout(1, 2, "TileSpacing", "tight"); 
nexttile; hold on; axis equal;
%plot(Mx(:,1), Mx(:,2), "b.");
plot(Ax(:,1), Ax(:,2), "r.-");
plot(Ax(Nx_idx(1:2), 1), Ax(Nx_idx(1:2),2), "go", "MarkerFaceColor", "g", "MarkerSize", 8);
plot(Ax(Nx_idx(3:end), 1), Ax(Nx_idx(3:end),2), "bo", "MarkerFaceColor", "b", "MarkerSize", 8);
nexttile; hold on; axis equal;
xline(0, "k");
yline(0, "k");
plot(DAx(:,1), DAx(:,2), "b.-"); 
%plot(DMx(:,1), DMx(:,2), "r.");
plot(DAx(Nx_idx,1), DAx(Nx_idx,2), "ro", "MarkerFaceColor", "r", "MarkerSize", 8);

PlotDeformationField(Mx, DMx, 100);
PlotOMesh(Mx, DMx, Ni, Nj);