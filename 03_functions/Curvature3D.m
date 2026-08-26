function [k1, k2] = Curvature3D(X, Y, Z)
    [Ni, Nj] = size(X);
    k1 = zeros(Ni, Nj);
    k2 = zeros(Ni, Nj);
    
    for i = 2:(Ni-1)
        for j = 2:(Nj-1)
            % 1. Extract the 3x3 local coordinate neighborhood
            P = [reshape(X(i-1:i+1, j-1:j+1), 9, 1), ...
                 reshape(Y(i-1:i+1, j-1:j+1), 9, 1), ...
                 reshape(Z(i-1:i+1, j-1:j+1), 9, 1)];
            
            % 2. Mean-shift
            P_shifted = P - [X(i,j), Y(i,j), Z(i,j)];
            
            % 3. PCA Alignment
            [coeff, ~, ~] = pca(P_shifted);
            
            if det(coeff) < 0
                coeff(:, 2) = -coeff(:, 2);
            end
            
            % 4. Rotate to local tangent space
            P_local = P_shifted * coeff;
            x_loc = P_local(:, 1);
            y_loc = P_local(:, 2);
            z_loc = P_local(:, 3);
            
            % 5. Fit quadric surface
            A_mat = [0.5 * x_loc.^2, x_loc .* y_loc, 0.5 * y_loc.^2, x_loc, y_loc];
            c = A_mat \ z_loc;
            
            % 6. Principal Curvatures
            A_coeff = c(1); B_coeff = c(2); C_coeff = c(3);
            T1 = (A_coeff + C_coeff) / 2;
            T2 = sqrt(T1^2 - A_coeff*C_coeff + B_coeff^2);
            
            k1(i, j) = T1 + T2;
            k2(i, j) = T1 - T2;
        end
    end
    
    % 7. Extrapolate boundaries
    k1(1, :) = k1(2, :); k1(end, :) = k1(end-1, :);
    k2(1, :) = k2(2, :); k2(end, :) = k2(end-1, :);
    k1(:, 1) = k1(:, 2); k1(:, end) = k1(:, end-1);
    k2(:, 1) = k2(:, 2); k2(:, end) = k2(:, end-1);
    
    % Corners
    k1(1, 1) = k1(2, 2); k1(1, end) = k1(2, end-1);
    k1(end, 1) = k1(end-1, 2); k1(end, end) = k1(end-1, end-1);
    k2(1, 1) = k2(2, 2); k2(1, end) = k2(2, end-1);
    k2(end, 1) = k2(end-1, 2); k2(end, end) = k2(end-1, end-1);
end