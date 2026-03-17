function PlotPointSeq(Ax, Nx_idx_final, num_base, N_IP, anim_delay)
%PLOTPOINTSEQ Animates the chronological selection of control points.
% Inputs:
%   Ax           : N x 2 matrix of Aerofoil Mesh Points
%   Nx_idx_final : Array of selected control point indices
%   num_base     : Number of base points initially defined (e.g., LE, TE)
%   N_IP         : Number of Initial Points selected by distance/curvature
%   anim_delay   : Configurable pause between point additions in seconds

arguments
    Ax           (:,2) double
    Nx_idx_final (:,1) double
    num_base     (1,1) double
    N_IP         (1,1) double
    anim_delay   (1,1) double = 0.2
end

%% FUNCTION BODY
num_total = length(Nx_idx_final);

figure;
hold on; 
axis equal;
title('Control Point Selection Sequence');
xlabel('X'); 
ylabel('Y');

% Plot the base aerofoil mesh in light gray for context
plot(Ax(:,1), Ax(:,2), '.-', 'Color', [0.7 0.7 0.7]); 

% Iterate through the chronological list of indices
for k = 1:num_total
    
    % Get the geometric coordinates for the current chronological point
    idx = Nx_idx_final(k);
    pt_x = Ax(idx, 1);
    pt_y = Ax(idx, 2);
    
    % Determine color based on when it was added to the list
    if k <= num_base
        pt_color = 'g'; % First points: Base (LE, TE)
    elseif k <= num_base + N_IP
        pt_color = 'r'; % Next N_IP points: Initial Selection
    else
        pt_color = 'b'; % Remaining points: Greedy Selection
    end
    
    % Plot the new point as a black 'x'
    h_new = plot(pt_x, pt_y, 'kx', 'MarkerSize', 10, 'LineWidth', 2);
    drawnow;
    
    % Pause for visibility
    pause(anim_delay);
    
    % Replace the 'x' with the permanent colored 'o'
    delete(h_new);
    plot(pt_x, pt_y, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', pt_color, 'MarkerSize', 8);
    drawnow;
end
hold off;

end