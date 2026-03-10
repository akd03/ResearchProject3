function [DAx] = Rotate(Ax, Ang, Pos)
%ROTATE Rotate Aerofoil Mesh
%   Inputs: Ax, Ang, Pos

arguments (Input)
    Ax           % Aerofoil surface mesh
    Ang          % Angle (degrees) CW +ve about left-hand leading edge
    Pos          % Chord Position to rotate about (0-1)
end

arguments (Output)
    DAx          % Deformed Aerofoil surface mesh
end

%% FUNCTION BODY
DAx = Translate(Ax, [-Pos, 0]);
theta = deg2rad(Ang);
R = [cos(theta) -sin(theta); 
     sin(theta) cos(theta)];
DAx = DAx * R;
DAx = Translate(DAx, [Pos, 0]);
end