function DAx = BendTwist3D(Ax, S_root, S_tip, Phi_tip_deg, Theta_tip_deg)
%BENDTWIST3D Applies inextensible bending and twist using a pre-defined elastic axis.
%
% Inputs:
%   Ax            - Surface mesh coordinates [M x 3]
%   S_root        - Quarter chord root position [1 x 3]
%   S_tip         - Quarter chord tip position [1 x 3]
%   Phi_tip_deg   - Maximum bending slope at the tip (degrees)
%   Theta_tip_deg - Maximum twist angle at the tip (degrees)
%
% Output:
%   DAx           - Deformed surface mesh coordinates [M x 3]

    M = size(Ax, 1);
    
    % Ensure spar inputs are row vectors for proper broadcasting
    if iscolumn(S_root); S_root = S_root'; end
    if iscolumn(S_tip); S_tip = S_tip'; end

    %% 1. DEFINE THE STRAIGHT SPAR VECTOR
    v_spar = S_tip - S_root;
    L_total = norm(v_spar);
    t_const = v_spar / L_total;
    
    % Define constant Normal and Chordwise vectors for the straight spar
    n_temp = cross([-1, 0, 0], t_const);
    n_const = n_temp / norm(n_temp);
    c_const = cross(n_const, t_const);
    
    %% 2. ORTHOGONAL PROJECTION OF SURFACE NODES
    p_rel = Ax - S_root;
    
    % Dot product to find physical distance along the spar
    d_mapped = p_rel * t_const'; 
    s_mapped = d_mapped / L_total;
    
    % Calculate exact mapped point on the straight spar
    S_mapped = S_root + d_mapped * t_const;
    
    % Relative position vector from spar to surface node
    r = Ax - S_mapped;
    
    %% 3. INEXTENSIBLE BENDING INTEGRATION
    % We integrate the bent tangent along a high-res 1D array
    num_int_pts = 1000;
    s_int = linspace(0, 1, num_int_pts)';
    ds = s_int(2) - s_int(1);
    
    Phi_tip_rad = deg2rad(Phi_tip_deg);
    phi_int = Phi_tip_rad .* s_int;
    
    cos_phi = cos(phi_int);
    sin_phi = sin(phi_int);
    
    % Rotate the constant tangent by the bending angle
    t_bent_int = t_const .* cos_phi + cross(repmat(c_const, num_int_pts, 1), repmat(t_const, num_int_pts, 1), 2) .* sin_phi + ...
                 c_const .* (sum(c_const .* t_const, 2) .* (1 - cos_phi));
    t_bent_int = t_bent_int ./ vecnorm(t_bent_int, 2, 2);
    
    % Integrate to find new spatial coordinates of the bent spar
    S_bent_int = zeros(num_int_pts, 3);
    S_bent_int(1, :) = S_root;
    for k = 2:num_int_pts
        S_bent_int(k, :) = S_bent_int(k-1, :) + ...
            (t_bent_int(k, :) + t_bent_int(k-1, :)) / 2 * (ds * L_total);
    end
    
%% 4. INTERPOLATE KINEMATICS FOR SURFACE NODES
    % Lookup the exact bent coordinates and tangents for the mapped points
    S_bent_mapped = interp1(s_int, S_bent_int, s_mapped, 'linear', 'extrap');
    t_bent_mapped = interp1(s_int, t_bent_int, s_mapped, 'linear', 'extrap');
    t_bent_mapped = t_bent_mapped ./ vecnorm(t_bent_mapped, 2, 2);
    
    % --- RIGID INBOARD CLAMP ---
    % Override the extrapolation for any points inboard of the root (s < 0)
    inboard_mask = s_mapped < 0;
    
    S_bent_mapped(inboard_mask, :) = S_mapped(inboard_mask, :);
    t_bent_mapped(inboard_mask, :) = repmat(t_const, sum(inboard_mask), 1);
    % ---------------------------
    
    %% 5. DOUBLE ROTATION
    % Clamp the spanwise parameter so inboard points experience 0 bend/twist
    s_def = max(s_mapped, 0);
    
    phi_mapped = Phi_tip_rad .* s_def;
    Theta_tip_rad = deg2rad(Theta_tip_deg);
    theta_mapped = Theta_tip_rad .* s_def;
    
    c_mapped = repmat(c_const, M, 1);
    
    % Rotation 1: Bend
    cos_phi_m = cos(phi_mapped);
    sin_phi_m = sin(phi_mapped);
    r_bent = r .* cos_phi_m + cross(c_mapped, r, 2) .* sin_phi_m + ...
             c_mapped .* (sum(c_mapped .* r, 2) .* (1 - cos_phi_m));
             
    % Rotation 2: Twist
    cos_theta = cos(theta_mapped);
    sin_theta = sin(theta_mapped);
    r_final = r_bent .* cos_theta + cross(t_bent_mapped, r_bent, 2) .* sin_theta + ...
              t_bent_mapped .* (sum(t_bent_mapped .* r_bent, 2) .* (1 - cos_theta));
              
    %% 6. FINAL COORDINATES
    DAx = S_bent_mapped + r_final;
end