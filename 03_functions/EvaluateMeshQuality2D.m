function [CellMetrics, MeshSummary] = EvaluateMeshQuality2D(Mx_def, Mx_base, Ni, Nj)
% EVALUATEMESHQUALITY2D Calculates Knupp and Orthogonality metrics for a 2D structured mesh.
%
% INPUTS:
%   Mx_def  - [M x 2] Array of deformed mesh coordinates
%   Mx_base - [M x 2] Array of initial/undeformed mesh coordinates (for reference sizes)
%   Ni, Nj  - Mesh dimensions
%
% OUTPUTS:
%   CellMetrics - Struct containing [Ni-1 x Nj-1] matrices for each metric
%   MeshSummary - Struct containing the global average for each metric

    % Define grid dimensions for cells
    num_cells_i = Ni - 1;
    num_cells_j = Nj - 1;
    
    % Generate cell indices
    [I, J] = ndgrid(1:num_cells_i, 1:num_cells_j);
    
    % Map (I,J) back to the 1D array indices for the 4 corners (CCW winding)
    % n0: Bottom-Left, n1: Bottom-Right, n2: Top-Right, n3: Top-Left
    idx0 = (J - 1) * Ni + I;
    idx1 = (J - 1) * Ni + (I + 1);
    idx2 = J * Ni + (I + 1);
    idx3 = J * Ni + I;
    
    % --- 1. CALCULATE REFERENCE SIZES (w) FROM BASE MESH ---
    % Used for the relative size metric calculation
    dx01_b = Mx_base(idx1, 1) - Mx_base(idx0, 1);
    dy01_b = Mx_base(idx1, 2) - Mx_base(idx0, 2);
    dx03_b = Mx_base(idx3, 1) - Mx_base(idx0, 1);
    dy03_b = Mx_base(idx3, 2) - Mx_base(idx0, 2);
    dx23_b = Mx_base(idx3, 1) - Mx_base(idx2, 1);
    dy23_b = Mx_base(idx3, 2) - Mx_base(idx2, 2);
    dx21_b = Mx_base(idx1, 1) - Mx_base(idx2, 1);
    dy21_b = Mx_base(idx1, 2) - Mx_base(idx2, 2);
    
    % Local areas at corner 0 and 2 for base mesh
    alpha0_b = dx01_b .* dy03_b - dy01_b .* dx03_b;
    alpha2_b = dx23_b .* dy21_b - dy23_b .* dx21_b;
    w_ref = 0.5 * (alpha0_b + alpha2_b); % Total reference area

    % --- 2. CALCULATE JACOBIAN COMPONENTS FOR DEFORMED MESH ---
    % Corner 0 (k=0)
    dx01 = Mx_def(idx1, 1) - Mx_def(idx0, 1); dy01 = Mx_def(idx1, 2) - Mx_def(idx0, 2);
    dx03 = Mx_def(idx3, 1) - Mx_def(idx0, 1); dy03 = Mx_def(idx3, 2) - Mx_def(idx0, 2);
    a0 = dx01 .* dy03 - dy01 .* dx03;
    L11_0 = dx01.^2 + dy01.^2; 
    L22_0 = dx03.^2 + dy03.^2;
    dot_0 = dx01 .* dx03 + dy01 .* dy03;

    % Corner 1 (k=1)
    dx12 = Mx_def(idx2, 1) - Mx_def(idx1, 1); dy12 = Mx_def(idx2, 2) - Mx_def(idx1, 2);
    dx10 = -dx01;                             dy10 = -dy01;
    a1 = dx12 .* dy10 - dy12 .* dx10;
    L11_1 = dx12.^2 + dy12.^2; 
    L22_1 = dx10.^2 + dy10.^2;
    dot_1 = dx12 .* dx10 + dy12 .* dy10;

    % Corner 2 (k=2)
    dx23 = Mx_def(idx3, 1) - Mx_def(idx2, 1); dy23 = Mx_def(idx3, 2) - Mx_def(idx2, 2);
    dx21 = -dx12;                             dy21 = -dy12;
    a2 = dx23 .* dy21 - dy23 .* dx21;
    L11_2 = dx23.^2 + dy23.^2; 
    L22_2 = dx21.^2 + dy21.^2;
    dot_2 = dx23 .* dx21 + dy23 .* dy21;

    % Corner 3 (k=3)
    dx30 = -dx03;                             dy30 = -dy03;
    dx32 = -dx23;                             dy32 = -dy23;
    a3 = dx30 .* dy32 - dy30 .* dx32;
    L11_3 = dx30.^2 + dy30.^2; 
    L22_3 = dx32.^2 + dy32.^2;
    dot_3 = dx30 .* dx32 + dy30 .* dy32;
    
    % --- 3. VALIDITY MASK ---
    % Elements with any local negative area are invalid and set to 0
    valid_mask = (a0 > 0) & (a1 > 0) & (a2 > 0) & (a3 > 0);

    % --- 4. KNUPP METRICS CALCULATIONS ---
    % Shape Metric
    sum_shape = (L11_0 + L22_0)./a0 + (L11_1 + L22_1)./a1 + ...
                (L11_2 + L22_2)./a2 + (L11_3 + L22_3)./a3;
    f_shape = 8 ./ sum_shape;
    f_shape(~valid_mask) = 0;

    % Skew Metric
    sum_skew = sqrt(L11_0 .* L22_0)./a0 + sqrt(L11_1 .* L22_1)./a1 + ...
               sqrt(L11_2 .* L22_2)./a2 + sqrt(L11_3 .* L22_3)./a3;
    f_skew = 4 ./ sum_skew;
    f_skew(~valid_mask) = 0;

    % Relative Size Metric
    area_def = 0.5 * (a0 + a2); 
    tau = area_def ./ w_ref;
    f_size = min(tau, 1./tau);
    f_size(area_def <= 0) = 0;

    % Combined Metrics
    f_size_shape = f_size .* f_shape;
    f_size_skew = f_size .* f_skew;

    % --- 5. ANGULAR ORTHOGONALITY METRIC ---
    % Calculates the maximum deviation from 90 degrees (pi/2) across all 4 corners
    theta0 = acos(dot_0 ./ sqrt(L11_0 .* L22_0));
    theta1 = acos(dot_1 ./ sqrt(L11_1 .* L22_1));
    theta2 = acos(dot_2 ./ sqrt(L11_2 .* L22_2));
    theta3 = acos(dot_3 ./ sqrt(L11_3 .* L22_3));
    
    max_dev = max(cat(3, abs(theta0 - pi/2), abs(theta1 - pi/2), ...
                         abs(theta2 - pi/2), abs(theta3 - pi/2)), [], 3);
                     
    % Normalize to [0, 1] where 1 is perfectly orthogonal (0 deviation)
    f_ortho = 1 - (max_dev ./ (pi/2));
    f_ortho(~valid_mask) = 0;

    % --- 6. ASSEMBLE OUTPUTS ---
    CellMetrics.Shape = f_shape;
    CellMetrics.Skew = f_skew;
    CellMetrics.Size = f_size;
    CellMetrics.SizeShape = f_size_shape;
    CellMetrics.SizeSkew = f_size_skew;
    CellMetrics.Orthogonality = f_ortho;
    
    % Calculate global averages for the mesh (excluding strictly zero/invalid cells)
    MeshSummary.AvgShape = mean(f_shape(valid_mask));
    MeshSummary.AvgSkew = mean(f_skew(valid_mask));
    MeshSummary.AvgSizeShape = mean(f_size_shape(valid_mask));
    MeshSummary.AvgSizeSkew = mean(f_size_skew(valid_mask));
    MeshSummary.AvgOrthogonality = mean(f_ortho(valid_mask));
endfunction [outputArg1,outputArg2] = aa(inputArg1,inputArg2)
%AA Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    inputArg1
    inputArg2
end

arguments (Output)
    outputArg1
    outputArg2
end

outputArg1 = inputArg1;
outputArg2 = inputArg2;
end