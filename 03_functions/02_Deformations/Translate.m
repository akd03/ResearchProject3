function [DAx] = Translate(Ax, Def)
%TRANSLATE Translate Aerofoil Mesh
%   Inputs: Ax, Def

arguments (Input)
    Ax           % Aerofoil surface mesh
    Def          % n-d Translation vector (same dimesion as Ax coords)
end

arguments (Output)
    DAx          % Deformed Aerofoil surface mesh
end

%% FUNCTION BODY
DAx = Ax + Def;
end