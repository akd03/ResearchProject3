function PlotDeformation3D(Ax, DAx)
%PLOTDEFORMATION3D Visualizes the original and deformed aerodynamic surfaces.
arguments (Input)
    Ax (:,3) double     
    DAx (:,3) double    
end

%% 1. INITIALIZATION AND DISPLACEMENT
M = size(Ax, 1);
Disp = sqrt(sum((DAx - Ax).^2, 2));

figure;
hold on; grid on; axis equal;

%% 2. DYNAMIC MESH PARSING
y_diffs = diff(Ax(:,2));
reset_idx = find(y_diffs < -0.1, 1); 

if isempty(reset_idx)
    blocks = {1:M};
else
    blocks = {1:reset_idx, (reset_idx+1):M};
end

h_orig = [];
h_def = [];

% Preallocate structures to hold the boundary edges for stitching
boundaries = struct('X', [], 'Y', [], 'Z', [], 'C', []);

%% 3. RESHAPE AND PLOT
for b = 1:length(blocks)
    idx = blocks{b};
    if isempty(idx), continue; end
    
    X_b = Ax(idx, 1); Y_b = Ax(idx, 2); Z_b = Ax(idx, 3);
    DX_b = DAx(idx, 1); DY_b = DAx(idx, 2); DZ_b = DAx(idx, 3);
    C_b  = Disp(idx); 
    
    y_b_diffs = diff(Y_b);
    Ni = find(abs(y_b_diffs) > 1e-4, 1);
    if isempty(Ni), Ni = length(idx); end
    Nj = length(idx) / Ni;
    
    if mod(length(idx), Ni) == 0
        DX_grid = reshape(DX_b, Ni, Nj);
        DY_grid = reshape(DY_b, Ni, Nj);
        DZ_grid = reshape(DZ_b, Ni, Nj);
        C_grid  = reshape(C_b, Ni, Nj);
        
        % Plot Original Mesh (Semi-transparent)
        p1 = surf(reshape(X_b, Ni, Nj), reshape(Y_b, Ni, Nj), reshape(Z_b, Ni, Nj), ...
            'FaceColor', [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
        
        % Plot Deformed Mesh
        p2 = surf(DX_grid, DY_grid, DZ_grid, C_grid, 'EdgeColor', 'k', 'FaceAlpha', 1.0);
        
        if b == 1
            h_orig = p1; h_def = p2;
        end
        
        % Store the first and last rows of this block to stitch the LE and TE later
        boundaries(b).X_first = DX_grid(1, :);   boundaries(b).X_last = DX_grid(end, :);
        boundaries(b).Y_first = DY_grid(1, :);   boundaries(b).Y_last = DY_grid(end, :);
        boundaries(b).Z_first = DZ_grid(1, :);   boundaries(b).Z_last = DZ_grid(end, :);
        boundaries(b).C_first = C_grid(1, :);    boundaries(b).C_last = C_grid(end, :);
    else
        scatter3(DX_b, DY_b, DZ_b, 5, C_b, 'filled');
    end
end

%% 4. STITCH THE LEADING AND TRAILING EDGES
% Only execute if exactly two structured blocks (Upper and Lower) were found
if length(blocks) == 2 && ~isempty(boundaries(1).X_first) && ~isempty(boundaries(2).X_first)
    
    % Dynamically find the shared number of spanwise sections to prevent index errors
    % if one block has a wingtip cap and the other does not.
    Nj_shared = min(length(boundaries(1).X_first), length(boundaries(2).X_first));
    
    % STITCH 1: Lower block end to Upper block start (Usually the Leading Edge)
    X_stitch_1 = [boundaries(1).X_last(1:Nj_shared); boundaries(2).X_first(1:Nj_shared)];
    Y_stitch_1 = [boundaries(1).Y_last(1:Nj_shared); boundaries(2).Y_first(1:Nj_shared)];
    Z_stitch_1 = [boundaries(1).Z_last(1:Nj_shared); boundaries(2).Z_first(1:Nj_shared)];
    C_stitch_1 = [boundaries(1).C_last(1:Nj_shared); boundaries(2).C_first(1:Nj_shared)];
    
    surf(X_stitch_1, Y_stitch_1, Z_stitch_1, C_stitch_1, 'EdgeColor', 'k', 'FaceAlpha', 1.0);
    
    % STITCH 2: Upper block end to Lower block start (Usually the Trailing Edge)
    X_stitch_2 = [boundaries(2).X_last(1:Nj_shared); boundaries(1).X_first(1:Nj_shared)];
    Y_stitch_2 = [boundaries(2).Y_last(1:Nj_shared); boundaries(1).Y_first(1:Nj_shared)];
    Z_stitch_2 = [boundaries(2).Z_last(1:Nj_shared); boundaries(1).Z_first(1:Nj_shared)];
    C_stitch_2 = [boundaries(2).C_last(1:Nj_shared); boundaries(1).C_first(1:Nj_shared)];
    
    surf(X_stitch_2, Y_stitch_2, Z_stitch_2, C_stitch_2, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
end

%% 5. LIGHTING AND FORMATTING
lighting gouraud;
camlight('headlight');
colormap('turbo');
cb = colorbar; cb.Label.String = 'Displacement Magnitude';

xlabel('X (Chordwise)'); ylabel('Y (Spanwise)'); zlabel('Z (Thickness)');
title('3D Aerofoil Surface Deformation');

if ~isempty(h_orig) && ~isempty(h_def)
    legend([h_orig, h_def], {'Original Mesh', 'Deformed Mesh'}, 'Location', 'northeast');
end
view(30, 30);
hold off;
end