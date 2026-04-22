%% MESH PROFILE VISUALIZER (X-Z PROJECTION)
%% 140k, 41 points on lower surface, 49 profiles. 
%% MESH PROFILE VISUALIZER: TOP AND BOTTOM SURFACE TEST
clear; clc; close all;

filename = '\05_meshes\01_3d_meshes\surfacepoints8M.plt';

%% 1. READ HEADER AND EXTRACT DATA
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

fprintf('Loading numeric data (skipping %d header lines)...\n', headerLines);
RawData = readmatrix(filename, 'NumHeaderLines', headerLines, 'FileType', 'text');
X = RawData(:, 1);
Z = RawData(:, 3); 
M = length(X);

%% 2. KNOWN LOWER SURFACE PARAMETERS
n_lower = 41;          
num_profiles = 49;     
total_lower_points = n_lower * num_profiles;

fprintf('Lower surface assumes: %d profiles of %d points (%d points total)\n', num_profiles, n_lower, total_lower_points);
fprintf('Points remaining for top surface: %d\n\n', M - total_lower_points);

%{ 
%% 3. TOP SURFACE GUESS
% ---> CHANGE THIS NUMBER to test the top surface <---
n_top = 39;  % Guessing points per top profile

%% 4. PLOT PROFILES IN X-Z PLANE
figure;
hold on; grid on; 

% --- Plot Lower Surface (Greyed out for reference) ---
for i = 48:49%1:49
    idx_start = (i-1)*n_lower + 1;
    % Safety bound to prevent indexing errors if the 50 profile guess is too large
    idx_end = min(i*n_lower, M); 
    
    plot(X(idx_start:idx_end), Z(idx_start:idx_end), 'Color', [0.8 0.8 0.8], 'LineWidth', 1);
end

% --- Plot Top Surface (Colored to check stacking and closure) ---
colors = lines(7); 
start_upper_idx = total_lower_points; 

for i = 48%1:num_profiles-1
    idx_start = start_upper_idx + (i-1)*n_top + 1;
    idx_end = start_upper_idx + i*n_top;
    
    % Prevent indexing past the end of the file
    if idx_end > M
        fprintf('WARNING: Reached end of file at upper profile %d! Stopping plot early.\n', i);
        idx_end = M;
        plot(X(idx_start:end), Z(idx_start:end), 'r.-', 'LineWidth', 2);
        break;
    end
    
    c = colors(mod(i-1, 7) + 1, :);
    
    % Plot the segment
    plot(X(idx_start:idx_end), Z(idx_start:idx_end), '.-', 'Color', c, 'LineWidth', 1.2, 'MarkerSize', 8);
    
    % Mark the exact FIRST point of each top profile (Look to see if this sits on the LE/TE)
    plot(X(idx_start), Z(idx_start), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4);
end

% Check how many points are left over
points_used = total_lower_points + (num_profiles * n_top);
if points_used > M
    fprintf('Calculation exceeded total points (%d > %d). Lower your guess.\n', points_used, M);
elseif points_used < M
    fprintf('WARNING: %d points left over at the end.\n', M - points_used);
else
    fprintf('SUCCESS: Top surface divided perfectly leaving 0 points!\n');
end

xlabel('X Coordinate');
ylabel('Z Coordinate');
title(sprintf('Top Surface Check (n_{top} = %d points)', n_top));
hold off;

%}
Ax = RawData(:, 1:3); 
DAx = Ax;
DAx(:,3) = DAx(:,3) + 0.05 * Ax(:,2); 
RAx = DAx - Ax;

PlotDeformation3D(Ax, DAx);