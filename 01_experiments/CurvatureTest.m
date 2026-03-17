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

%% RELATIVE ERROR CALCULATION
% Calculate the normalized error in curvature
rel_err_kappa = abs(kappa_num(idx_top) - kappa_ana) ./ kappa_ana;

%% PLOT COMPARISON
figure;
tiledlayout(1, 3, "TileSpacing", "compact");

% --- TILE 1: CURVATURE ---
ax1 = nexttile;
hold on; grid on;
set(ax1, 'YScale', 'log'); 
plot(x_top, kappa_ana, "r-", "LineWidth", 1.5, "DisplayName", "Analytical Curvature");
plot(x_top, kappa_num(idx_top), "b-", "LineWidth", 1.5, "DisplayName", "Numerical Curvature (CDS)");
xlabel('Chord (x)');
ylabel('Curvature (\kappa)');
title('Curvature Comparison: Top Surface');
legend;
ylim([0.1, 100]); % Adjusted lower limit slightly for log scale
xlim([0, 1]);

%% RADIUS OF CURVATURE COMPARISON
% 1. Calculate Numerical Radius
r_num = 1 ./ kappa_num(idx_top);
% 2. Calculate Analytical Radius
r_ana = 1 ./ kappa_ana;
% 3. Override the Leading Edge Singularity (x = 0)
t = 0.12;
r_ana(1) = 1.1019 * t^2; 

% --- TILE 2: RADIUS ---
ax2 = nexttile;
hold on; grid on;
plot(x_top, r_ana, "r-", "LineWidth", 1.5, "DisplayName", "Analytical Radius");
plot(x_top, r_num, "b-", "LineWidth", 1.5, "DisplayName", "Numerical Radius (CDS)");
xlabel('Chord (x)');
ylabel('Radius of Curvature (r)');
title('Radius of Curvature Comparison');
legend('Location', 'northwest');
ylim([0, 8.5]);
xlim([0, 1]);

% --- TILE 3: NORMALIZED ERROR ---
ax3 = nexttile;
hold on; grid on;
plot(x_top, rel_err_kappa*100, "k-", "LineWidth", 1.5, "DisplayName", "Relative Error (\kappa)");
xlabel('Chord (x)');
ylabel('Percentage Error (%)');
title('Normalized Curvature Error');
legend('Location', 'best');
xlim([0, 1]);

%% FULL AEROFOIL NUMERICAL CURVATURE PLOT
figure;
hold on; grid on;

% Define the full index range
idx_full = 1:length(kappa_num);

% Plot the full numerical curvature array
plot(idx_full, kappa_num, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Numerical Curvature');

% Add vertical reference lines for geometric landmarks
% Assuming a 257-point O-mesh starting at the lower TE and wrapping to the upper TE
xline(1, 'k--', 'Lower TE', 'LabelVerticalAlignment', 'bottom');
xline(129, 'r--', 'Leading Edge', 'LabelVerticalAlignment', 'bottom');
xline(length(kappa_num), 'k--', 'Upper TE', 'LabelVerticalAlignment', 'bottom');

% Format the plot
set(gca, 'YScale', 'log'); 
xlabel('Node Index (Lower TE \rightarrow LE \rightarrow Upper TE)');
ylabel('Curvature (\kappa)');
title('Numerical Curvature Across Entire O-Mesh');
xlim([1, length(kappa_num)]);
% Set a reasonable lower limit to avoid graphing floating-point noise near zero
ylim([1e-2, 100]); 
hold off;