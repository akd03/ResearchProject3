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
% Setup 
N = size(xN, 1);
idx_N = zeros(N + K, 1);
for i = 1:N
    [~, idx] = min(sum((Ax - xN(i, :)).^2, 2));
    idx_N(i) = idx;
end

% 2. Initialize the running minimum distance array
% Extract columns for all currently selected points and find the minimum row-wise
MinDist_Ax = min(Norm_Ax(:, idx_N(1:N)), [], 2);

% 3. Iterative Selection Loop
for i = 1:K
    % Find the point that maximizes the minimum distance
    [~, next_idx] = max(MinDist_Ax);
    
    % Store the new index
    current_step = N + i;
    idx_N(current_step) = next_idx;
    
    % Update the running minimum distance directly from the Norm_Ax matrix
    MinDist_Ax = min(MinDist_Ax, Norm_Ax(:, next_idx));
end

% 4. Extract final coordinates
xN_IP = Ax(idx_N, :);


end