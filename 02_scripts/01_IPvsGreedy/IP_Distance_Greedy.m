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
DAx = Rotate(Ax, 10, 1);




%% INITIAL POINT SELECTION
Nx_idx = IP_Distance(Nx_idx, Ax, 15);

figure; 
tiledlayout(1, 2, "TileSpacing", "tight"); 
nexttile; hold on; axis equal;
plot(Ax(:,1), Ax(:,2), "r.-");
plot(Ax(Nx_idx(1:2), 1), Ax(Nx_idx(1:2),2), "go", "MarkerFaceColor", "g", "MarkerSize", 8);
plot(Ax(Nx_idx(3:end), 1), Ax(Nx_idx(3:end),2), "bo", "MarkerFaceColor", "b", "MarkerSize", 8);
nexttile; hold on; axis equal;
xline(0, "k");
yline(0, "k");
plot(DAx(:,1), DAx(:,2), "b.-"); 

%% GREEDY ALGORITHM




%xN = Ax(xN_idx, :);





%% RESULTS/PLOTTING