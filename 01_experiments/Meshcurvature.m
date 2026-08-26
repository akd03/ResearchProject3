clear; clc; 

filename = '\05_meshes\01_3d_meshes\surfacepoints140k.plt';
fid = fopen(filename, 'r');
headerLines = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ~isempty(regexp(line, '^\s*[-0-9.eE]+', 'once'))
        break; 
    end
    headerLines = headerLines + 1;
end
fclose(fid);

fprintf('Loading numeric data...\n');
RawData = readmatrix(filename, 'NumHeaderLines', headerLines, 'FileType', 'text');

Points = RawData(:, 1:3); 
%% 1. MESH PREPROCESSING AND MERGING
% Load your raw 3881x3 points here
% Points = ... 

% Corrected Split: Bottom surface comes first (41 * 49 = 2009 points)
% Followed by the Top surface (39 * 48 = 1872 points)
bot_points = Points(1:2009, :);
top_points = Points(2010:end, :);

% Reshape into structured grids
X_bot_full = reshape(bot_points(:,1), 41, 49);
Y_bot_full = reshape(bot_points(:,2), 41, 49);
Z_bot_full = reshape(bot_points(:,3), 41, 49);

X_top = reshape(top_points(:,1), 39, 48);
Y_top = reshape(top_points(:,2), 39, 48);
Z_top = reshape(top_points(:,3), 39, 48);

% Truncate the extra 49th spanwise profile on the bottom to allow merging
X_bot = X_bot_full(:, 1:48);
Y_bot = Y_bot_full(:, 1:48);
Z_bot = Z_bot_full(:, 1:48);

% Merge into a continuous C-grid (TE Bottom -> LE -> TE Top)
% Note: flipud is used assuming the raw profiles go LE -> TE. 
% This connects them smoothly at the Leading Edge.
X_merged = [flipud(X_bot); X_top(2:end, :)]; 
Y_merged = [flipud(Y_bot); Y_top(2:end, :)];
Z_merged = [flipud(Z_bot); Z_top(2:end, :)];

%% 2. CALCULATE CURVATURE ON MERGED GRID
[k1_merged, k2_merged] = Curvature3D(X_merged, Y_merged, Z_merged);

% Calculate Maximum Absolute Principal Curvature
K_max_merged = max(abs(k1_merged), abs(k2_merged));

%% 3. SPLIT DATA FOR PLOTTING
% The merged grid has length 41 (bot) + 39 (top) - 1 (shared LE) = 79 points per profile
% Re-split them to plot top and bottom separately
K_max_bot = flipud(K_max_merged(1:41, :)); % Flip back to original orientation
K_max_top = K_max_merged(41:end, :);       % Index 41 is the shared LE

%% 4. PLOT X-Y COLORMAPS
figure('Name', 'Maximum Absolute Principal Curvature', 'Position', [100, 100, 1200, 500]);
tiledlayout(1, 2, 'TileSpacing', 'compact');

% Plot Top Surface
ax1 = nexttile;
surf(ax1, X_top, Y_top, Z_top, K_max_top, 'EdgeColor', 'none');
view(ax1, 2); % Top-down X-Y view
shading(ax1, 'interp');
colormap(ax1, 'jet');
colorbar(ax1);
title(ax1, 'Top Surface Curvature');
xlabel(ax1, 'Chordwise (X)');
ylabel(ax1, 'Spanwise (Y)');
axis(ax1, 'equal');
grid(ax1, 'off');

% Plot Bottom Surface
ax2 = nexttile;
surf(ax2, X_bot, Y_bot, Z_bot, K_max_bot, 'EdgeColor', 'none');
view(ax2, 2); % Top-down X-Y view
shading(ax2, 'interp');
colormap(ax2, 'jet');
colorbar(ax2);
title(ax2, 'Bottom Surface Curvature');
xlabel(ax2, 'Chordwise (X)');
ylabel(ax2, 'Spanwise (Y)');
axis(ax2, 'equal');
grid(ax2, 'off');

% Link axes for easy zooming/panning
linkaxes([ax1, ax2], 'xy');


