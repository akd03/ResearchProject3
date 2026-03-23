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
Lx = [linspace(0, 1, 41)', zeros(41, 1)]; 
R = 3;
xc = 0.25;
theta_Lx = (Lx(end,1) - Lx(1,1)) / R;

t_Lx = (Lx(:,1) - xc) ./ (Lx(end,1)-Lx(1,1));
DLx = [R * sin(theta_Lx * t_Lx) + xc, R * cos(theta_Lx * t_Lx) - R];


figure; hold on; axis equal;
plot(Lx(:,1), Lx(:,2), "bo-"); 
plot(DLx(:,1), DLx(:,2), "ro-"); 


%% DEFINE NACA 2412 CAMBER FUNCTIONS
m = 0.05;
p = 0.7;

% Camber line equation (Piecewise)
yc_func = @(x) (x < p) .* (m/p^2 .* (2*p.*x - x.^2)) + ...
               (x >= p) .* (m/(1-p)^2 .* ((1-2*p) + 2*p.*x - x.^2));

% Gradient equation (Piecewise)
dyc_func = @(x) (x < p) .* (2*m/p^2 .* (p - x)) + ...
                (x >= p) .* (2*m/(1-p)^2 .* (p - x));

%% DEFORM MESH
DAx_NACA = GeneralCamber(Ax, yc_func, dyc_func);

% Plot to verify
figure; hold on; axis equal; grid on;
plot(Ax(:,1), Ax(:,2), 'b-', 'DisplayName', 'Original 0012');
plot(DAx_NACA(:,1), DAx_NACA(:,2), 'r-', 'DisplayName', 'Cambered 2412');
plot(Lx(:,1), Lx(:,2), "b-");
plot(Lx(:,1), yc_func(Lx(:,1)), "r-");
legend;