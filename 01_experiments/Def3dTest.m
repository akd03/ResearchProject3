%% VARIABLE NAMING CONVENTION                                              
%  Nx       Known (Control) Points                                         
%  Nx_idx   Nx indices in Ax                                               
%  Ax       Aerofoil Mesh Points                                           
%  RAx      Applied Deformation                                            
%  DAx      Actual Aerofoil Position After Rx                              

%% LOAD 3D MESH       
%
clear; clc; 

filename = '\05_meshes\01_3d_meshes\surfacepoints140k.plt';
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

fprintf('Loading numeric data...\n');
RawData = readmatrix(filename, 'NumHeaderLines', headerLines, 'FileType', 'text');

Ax = RawData(:, 1:3); 
DAx = BendTwist3D(Ax, [24.9,0,-3.28], [47, 35.4, 1.55], 0, -20);


PlotDeformation3D(Ax, DAx);
