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
DAx_NACA = NACACamber(Ax, 0.02, 0.5);

% Plot to verify
figure; hold on; axis equal; grid on;
plot(Ax(:,1), Ax(:,2), 'b-', 'DisplayName', 'Original 0012');
plot(DAx_NACA(:,1), DAx_NACA(:,2), 'r-', 'DisplayName', 'Cambered 2412');
legend;



%% CIRCULAR CAMBER (Fixed LE/TE at y=0)
R_circ = 2.5; % Target bend radius

% 1. Calculate the center of the circle (h, k)
h = 0.5;
% k must be negative so the arc bows upwards
k = -sqrt(R_circ^2 - h^2); 

% 2. Define the exact circle equation and its derivative
yc_circ = @(x) k + sqrt(R_circ^2 - (x - h).^2);
dyc_circ = @(x) (h - x) ./ sqrt(R_circ^2 - (x - h).^2);

% 3. Apply the deformation using your generalized function
DAx_Circ = GeneralCamber(Ax, yc_circ, dyc_circ);

%% PLOT CIRCULAR CAMBER
figure; 
hold on; axis equal; grid on;

% Plot Aerofoils
plot(Ax(:,1), Ax(:,2), 'b-', 'DisplayName', 'Original 0012');
plot(DAx_Circ(:,1), DAx_Circ(:,2), 'm-', 'DisplayName', sprintf('Circular Camber (R=%g)', R_circ));

% Plot Camber Lines
plot(Lx(:,1), Lx(:,2), "b-", 'HandleVisibility', 'off'); % Straight chord
plot(Lx(:,1), yc_circ(Lx(:,1)), "m--", 'LineWidth', 1.5, 'DisplayName', 'Circular Chord Line');

xlabel('X'); ylabel('Y');
title('Circular Camber Deformation (Fixed LE/TE)');
legend('Location', 'best');
hold off;