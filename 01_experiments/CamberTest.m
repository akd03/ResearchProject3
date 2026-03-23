%% LOAD TRUE NACA MESH

MeshFile = fopen("05_meshes\NACA0012257x129.xyz", "r");
Header = fscanf(MeshFile, '%d', 2);
Ni = Header(1);
Nj = Header(2);
MeshData = textscan(MeshFile, '%f %f %f');
fclose(MeshFile);

% Extract Aerofoil Points
Mx = [MeshData{1}, MeshData{3}];
Ax = Mx(1:Ni, 1:2);


%% Deforming Line Along Arc
Lx = [linspace(0, 1, 11)', zeros(11, 1)]; 
R = 3;
xc = 0.25;
theta_Lx = (Lx(end,1) - Lx(1,1)) / R;

t_Lx = (Lx(:,1) - xc) ./ (Lx(end,1)-Lx(1,1));
DLx = [R * sin(theta_Lx * t_Lx) + xc, R * cos(theta_Lx * t_Lx) - R];


figure; hold on; axis equal;
plot(Lx(:,1), Lx(:,2), "bo-"); 
plot(DLx(:,1), DLx(:,2), "ro-"); 