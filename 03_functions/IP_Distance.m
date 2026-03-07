function [xN_IP] = IP_Distance(xN, Ax, Norm_Ax, K);
%IP_DISTANCE Select K initial points based on euclidean distance.
%   Inputs: xN, Ax, Norm_Ax, Ax, K
%   xN may contain just LE and TE to begin with. 
%   K must be smaller than Ax
arguments (Input)
    xN                   % Current List of Initial Points
    Ax                   % Complete List of Aerofoil Points
    Norm_Ax              % Norms Between All Aerofoil Points
    K                    % Number of Initial Points to Select
end

arguments (Output)
    xN_IP                % New List of Initial Points
end

%% FUNCTION BODY
num_existing = size(xN, 1);

% 1. Map existing xN coordinates to their row indices in Ax
% Preallocate the array for existing points + K new points
selected_idx = zeros(num_existing + K, 1);
for i = 1:num_existing
    [~, idx] = min(sum((Ax - xN(i, :)).^2, 2));
    selected_idx(i) = idx;
end

% 2. Initialize the running minimum distance array
% Extract columns for all currently selected points and find the minimum row-wise
min_dist = min(Norm_Ax(:, selected_idx(1:num_existing)), [], 2);

% 3. Iterative Selection Loop
for i = 1:K
    % Find the point that maximizes the minimum distance
    [~, next_idx] = max(min_dist);
    
    % Store the new index
    current_step = num_existing + i;
    selected_idx(current_step) = next_idx;
    
    % Update the running minimum distance directly from the Norm_Ax matrix
    min_dist = min(min_dist, Norm_Ax(:, next_idx));
end

% 4. Extract final coordinates
xN_IP = Ax(selected_idx, :);


end