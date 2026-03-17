%% NACA 0012 ANALYTICAL FUNCTIONS
NACA = @(x) 5 * 0.12 * (0.2969.*x.^0.5 - 0.1260.*x - 0.3516.*x.^2 + 0.2843.*x.^3 - 0.1036.*x.^4);
dNACAdx = @(x) 5 * 0.12 * (0.5*0.2969.*x.^-0.5 - 0.1260 - 2*0.3516.*x + 3*0.2843.*x.^2 - 4*0.1036.*x.^3);
d2NACAdx2 = @(x) 5 * 0.12 * (-0.25*0.2969.*x.^-1.5 - 2*0.3516 + 6*0.2843.*x - 12*0.1036.*x.^2);

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

%% NUMERICAL CURVATURE (CDS)
% Calculate the full loop to maintain smooth transitions at the LE
[grad1, grad2] = CentralDifferenceDeriv(Ax);

% Numerical kappa using the parametric equation
kappa_num = abs(grad1(:,1) .* grad2(:,2) - grad1(:,2) .* grad2(:,1)) ./ ...
            ((grad1(:,1).^2 + grad1(:,2).^2).^1.5);

%% ANALYTICAL CURVATURE (TOP SURFACE EXACT MATCH)
% Extract just the upper surface (Leading Edge to Trailing Edge)
% In your O-mesh, LE is index 129, upper surface TE is index 257.
idx_top = 129:257;
x_top = Ax(idx_top, 1);

% Evaluate analytical equations at the EXACT x-coordinates of the mesh
y_prime = dNACAdx(x_top);
y_double_prime = d2NACAdx2(x_top);

% Analytical kappa using the explicit equation
kappa_ana = abs(y_double_prime) ./ (1 + y_prime.^2).^1.5;

%% PLOT COMPARISON
figure;
tiledlayout(1, 2, "TileSpacing", "compact");

ax1 = nexttile;
hold on; grid on;
set(ax1, 'YScale', 'log'); 

% Plot against physical x-coordinates
plot(x_top, kappa_ana, "r-", "LineWidth", 2, "DisplayName", "Analytical Curvature");
plot(x_top, kappa_num(idx_top), "b--", "LineWidth", 2, "DisplayName", "Numerical Curvature (CDS)");

xlabel('Chord Position (x)');
ylabel('Curvature (\kappa)');
title('Curvature Comparison: Top Surface of NACA 0012');
legend;

% Limit the y-axis to see the main chord body clearly 
% (Otherwise the massive leading edge spike dominates the visual scale)
ylim([0, 100]);
xlim([0, 1]);

%% RADIUS OF CURVATURE COMPARISON

% 1. Calculate Numerical Radius
r_num = 1 ./ kappa_num(idx_top);

% 2. Calculate Analytical Radius
r_ana = 1 ./ kappa_ana;

% 3. Override the Leading Edge Singularity (x = 0)
% The mathematical limit for a NACA 4-digit LE radius is 1.1019 * t^2
% Assuming idx_top(1) is exactly the leading edge where x = 0
t = 0.12;
r_ana(1) = 1.1019 * t^2; 

%% PLOT RADIUS COMPARISON
ax2 = nexttile;
hold on; grid on;


% Plot against physical x-coordinates
plot(x_top, r_ana, "r-", "LineWidth", 1.5, "DisplayName", "Analytical Radius");
plot(x_top, r_num, "b-", "LineWidth", 1.5, "DisplayName", "Numerical Radius (CDS)");

xlabel('Chord Position (x)');
ylabel('Radius of Curvature (r)');
title('Radius of Curvature Comparison: Top Surface of NACA 0012');
legend('Location', 'northwest');

% Limit the y-axis. The radius goes to infinity as the aerofoil 
% flattens out towards the trailing edge. Limiting to y=2 keeps 
% the leading edge and mid-chord features readable.
ylim([0, 8.5]);
xlim([0, 1]);