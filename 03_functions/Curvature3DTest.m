function [k1, k2] = Curvature3DTest(X, Y, Z)
    % Obtains the dimensions of the structured mesh
    [Ni, Nj] = size(X);
    
    % Initialize curvature matrices
    k1 = zeros(Ni, Nj);
    k2 = zeros(Ni, Nj);
    
    % Loop through the internal nodes to extract 3x3 neighborhoods
    for i = 2:(Ni-1)
        for j = 2:(Nj-1)
            % 1. Extract the 3x3 local coordinate neighborhood
            x_patch = X(i-1:i+1, j-1:j+1);
            y_patch = Y(i-1:i+1, j-1:j+1);
            z_patch = Z(i-1:i+1, j-1:j+1);
            
            % Reshape into a 9x3 matrix of points
            P = [x_patch(:), y_patch(:), z_patch(:)];
            
            % 2. Mean-shift the neighborhood to the central evaluation point
            p_center = [X(i,j), Y(i,j), Z(i,j)];
            P_shifted = P - p_center;
            
            % 3. PCA Alignment to define the local reference frame
            % The third principal component corresponds to the surface normal
            [coeff, ~, ~] = pca(P_shifted);
            
            % Enforce a right-handed coordinate system
            if det(coeff) < 0
                coeff(:, 2) = -coeff(:, 2);
            end
            
            % Ensure the normal points upwards relative to the global Z-axis
            % This prevents the principal curvature signs from arbitrarily flipping
            if coeff(3, 3) < 0
                coeff(:, 3) = -coeff(:, 3);
                coeff(:, 2) = -coeff(:, 2); % Maintain right-handedness
            end
            
            % Rotate the neighborhood points into the local tangent space
            P_local = P_shifted * coeff;
            x_loc = P_local(:, 1);
            y_loc = P_local(:, 2);
            z_loc = P_local(:, 3);
            
            % 4. Construct the overdetermined design matrix A (9x5)
            A_mat = [0.5 * x_loc.^2, x_loc .* y_loc, 0.5 * y_loc.^2, x_loc, y_loc];
            
            % 5. Solve the linear least-squares system
            c = A_mat \ z_loc;
            
            % Extract the second fundamental form coefficients
            A_coeff = c(1);
            B_coeff = c(2);
            C_coeff = c(3);
            
            % 6. Calculate Principal Curvatures
            T1 = (A_coeff + C_coeff) / 2;
            T2 = sqrt(T1^2 - A_coeff*C_coeff + B_coeff^2);
            
            k1(i, j) = T1 + T2;
            k2(i, j) = T1 - T2;
        end
    end
    
    % 7. Handle Boundaries (Simple replication from the nearest internal node)
    % Top and bottom rows
    k1(1, :) = k1(2, :);
    k1(end, :) = k1(end-1, :);
    k2(1, :) = k2(2, :);
    k2(end, :) = k2(end-1, :);
    
    % Left and right columns
    k1(:, 1) = k1(:, 2);
    k1(:, end) = k1(:, end-1);
    k2(:, 1) = k2(:, 2);
    k2(:, end) = k2(:, end-1);
    
    % Overwrite the four extreme corners
    k1(1, 1) = k1(2, 2); 
    k1(1, end) = k1(2, end-1);
    k1(end, 1) = k1(end-1, 2); 
    k1(end, end) = k1(end-1, end-1);
    
    k2(1, 1) = k2(2, 2); 
    k2(1, end) = k2(2, end-1);
    k2(end, 1) = k2(end-1, 2); 
    k2(end, end) = k2(end-1, end-1);
end