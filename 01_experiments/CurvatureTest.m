NACA = @(x) 5 * 0.12 * (0.2969 * x.^0.5 ...
                      - 0.1260 * x ...
                      - 0.3516 * x.^2 ...
                      + 0.2843 * x.^3 ...
                      - 0.1036 * x.^4);
dNACAdx = @(x) 5 * 0.12 * (0.5 * 0.2969 * x.^-0.5 ...
                           - 0.1260 ...
                           - 2 * 0.3516 * x ...
                           + 3 * 0.2843 * x.^2 ...
                           - 4 * 0.1036 * x.^3);


%% LOAD TRUE NACA MESH
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

%% CALC CONTINUOUS NACA
% Define x as a column vector
x = (0:0.001:1)'; 
N12 = NACA(x);
dN12 = dNACAdx(x);

% 1. Lower Surface (Trailing Edge to Leading Edge)
% Reverse the arrays so x goes from 1 down to 0, and make y negative
x_lower = flipud(x);
y_lower = -flipud(N12);

% 2. Upper Surface (Leading Edge to Trailing Edge)
% Keep arrays going from 0 up to 1, and keep y positive
x_upper = x;
y_upper = N12;

% 3. Combine into a continuous clockwise list
% We start at x_upper(2) to prevent duplicating the exact Leading Edge point (0,0)
Ax_exact = [x_lower, y_lower; x_upper(2:end), y_upper(2:end)];

% PLOT
figure; 
hold on; 
axis equal;
% Plot your full O-mesh boundary in blue
plot(Ax(:, 1), Ax(:, 2), "b.-", "DisplayName", "O-Mesh");
% Plot the exact continuous equation in red
plot(Ax_exact(:, 1), Ax_exact(:, 2), "r-", "LineWidth", 1.5, "DisplayName", "Exact NACA 0012");
legend;